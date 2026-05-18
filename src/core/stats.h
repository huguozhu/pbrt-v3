
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_STATS_H
#define PBRT_CORE_STATS_H

// core/stats.h*
// 统计系统模块：提供pbrt渲染器的性能统计和性能分析基础设施。
// 包含计数器、内存计数器、整数/浮点数分布统计、百分比/比率统计，
// 以及基于profiling category的性能分析支持。
#include "pbrt.h"
#include <map>
#include <cfloat>
#include <chrono>
#include <string>
#include <functional>
#include <mutex>

namespace pbrt {

// Statistics Declarations
// 统计类型声明：统计系统的核心类型前向声明和注册器

// StatRegisterer: 统计回调注册器。
// 通过静态初始化机制自动注册各统计变量的报告回调函数，
// 在渲染结束时由StatsAccumulator统一调用以汇总统计结果。
class StatsAccumulator;
class StatRegisterer {
  public:
    // StatRegisterer Public Methods
    // StatRegisterer 公有方法
    // 构造函数：注册一个统计报告回调函数
    StatRegisterer(std::function<void(StatsAccumulator &)> func) {
        static std::mutex mutex;
        std::lock_guard<std::mutex> lock(mutex);
        if (!funcs)
            funcs = new std::vector<std::function<void(StatsAccumulator &)>>;
        funcs->push_back(func);
    }
    // CallCallbacks: 调用所有已注册的回调函数，将统计数据汇总到累加器中
    static void CallCallbacks(StatsAccumulator &accum);

  private:
    // StatRegisterer Private Data
    // StatRegisterer 私有数据
    // funcs: 存储所有统计报告回调函数的静态向量指针
    static std::vector<std::function<void(StatsAccumulator &)>> *funcs;
};

// PrintStats: 打印所有统计结果到指定文件流
void PrintStats(FILE *dest);
// ClearStats: 清除所有统计信息
void ClearStats();
// ReportThreadStats: 报告当前线程的统计信息
void ReportThreadStats();

// StatsAccumulator: 统计累加器。
// 用于汇总所有线程的统计数据，支持多种统计类型：
// 计数器、内存计数器、整数/浮点数分布、百分比、比率。
class StatsAccumulator {
  public:
    // StatsAccumulator Public Methods
    // StatsAccumulator 公有方法

    // ReportCounter: 报告计数器的累加值（如渲染调用次数）
    void ReportCounter(const std::string &name, int64_t val) {
        counters[name] += val;
    }
    // ReportMemoryCounter: 报告内存计数器的值（如纹理占用的内存量）
    void ReportMemoryCounter(const std::string &name, int64_t val) {
        memoryCounters[name] += val;
    }
    // ReportIntDistribution: 报告整数型分布统计（总和、计数、最小值、最大值）
    void ReportIntDistribution(const std::string &name, int64_t sum,
                               int64_t count, int64_t min, int64_t max) {
        intDistributionSums[name] += sum;
        intDistributionCounts[name] += count;
        if (intDistributionMins.find(name) == intDistributionMins.end())
            intDistributionMins[name] = min;
        else
            intDistributionMins[name] =
                std::min(intDistributionMins[name], min);
        if (intDistributionMaxs.find(name) == intDistributionMaxs.end())
            intDistributionMaxs[name] = max;
        else
            intDistributionMaxs[name] =
                std::max(intDistributionMaxs[name], max);
    }
    // ReportFloatDistribution: 报告浮点型分布统计（总和、计数、最小值、最大值）
    void ReportFloatDistribution(const std::string &name, double sum,
                                 int64_t count, double min, double max) {
        floatDistributionSums[name] += sum;
        floatDistributionCounts[name] += count;
        if (floatDistributionMins.find(name) == floatDistributionMins.end())
            floatDistributionMins[name] = min;
        else
            floatDistributionMins[name] =
                std::min(floatDistributionMins[name], min);
        if (floatDistributionMaxs.find(name) == floatDistributionMaxs.end())
            floatDistributionMaxs[name] = max;
        else
            floatDistributionMaxs[name] =
                std::max(floatDistributionMaxs[name], max);
    }
    // ReportPercentage: 报告百分比统计（分子/分母，如已渲染像素占比）
    void ReportPercentage(const std::string &name, int64_t num, int64_t denom) {
        percentages[name].first += num;
        percentages[name].second += denom;
    }
    // ReportRatio: 报告比率统计（分子/分母，如光线与交点的比率）
    void ReportRatio(const std::string &name, int64_t num, int64_t denom) {
        ratios[name].first += num;
        ratios[name].second += denom;
    }

    // Print: 将所有统计数据输出到指定文件流
    void Print(FILE *file);
    // Clear: 清除所有累积的统计数据
    void Clear();

  private:
    // StatsAccumulator Private Data
    // StatsAccumulator 私有数据：各统计类型的存储映射表
    std::map<std::string, int64_t> counters;
    std::map<std::string, int64_t> memoryCounters;
    std::map<std::string, int64_t> intDistributionSums;
    std::map<std::string, int64_t> intDistributionCounts;
    std::map<std::string, int64_t> intDistributionMins;
    std::map<std::string, int64_t> intDistributionMaxs;
    std::map<std::string, double> floatDistributionSums;
    std::map<std::string, int64_t> floatDistributionCounts;
    std::map<std::string, double> floatDistributionMins;
    std::map<std::string, double> floatDistributionMaxs;
    std::map<std::string, std::pair<int64_t, int64_t>> percentages;
    std::map<std::string, std::pair<int64_t, int64_t>> ratios;
};

// Prof: 性能分析分类枚举。
// 定义了pbrt渲染管线中各阶段的性能分析类别，
// 用于精确测量各模块的执行时间。
enum class Prof {
    SceneConstruction,
    AccelConstruction,
    TextureLoading,
    MIPMapCreation,

    IntegratorRender,
    SamplerIntegratorLi,
    SPPMCameraPass,
    SPPMGridConstruction,
    SPPMPhotonPass,
    SPPMStatsUpdate,
    BDPTGenerateSubpath,
    BDPTConnectSubpaths,
    LightDistribLookup,
    LightDistribSpinWait,
    LightDistribCreation,
    DirectLighting,
    BSDFEvaluation,
    BSDFSampling,
    BSDFPdf,
    BSSRDFEvaluation,
    BSSRDFSampling,
    PhaseFuncEvaluation,
    PhaseFuncSampling,
    AccelIntersect,
    AccelIntersectP,
    LightSample,
    LightPdf,
    MediumSample,
    MediumTr,
    TriIntersect,
    TriIntersectP,
    CurveIntersect,
    CurveIntersectP,
    ShapeIntersect,
    ShapeIntersectP,
    ComputeScatteringFuncs,
    GenerateCameraRay,
    MergeFilmTile,
    SplatFilm,
    AddFilmSample,
    StartPixel,
    GetSample,
    TexFiltTrilerp,
    TexFiltEWA,
    TexFiltPtex,
    NumProfCategories
};

static_assert((int)Prof::NumProfCategories <= 64,
              "No more than 64 profiling categories may be defined.");

inline uint64_t ProfToBits(Prof p) { return 1ull << (int)p; }

static const char *ProfNames[] = {
    "Scene parsing and creation",
    "Acceleration structure creation",
    "Texture loading",
    "MIP map generation",

    "Integrator::Render()",
    "SamplerIntegrator::Li()",
    "SPPM camera pass",
    "SPPM grid construction",
    "SPPM photon pass",
    "SPPM photon statistics update",
    "BDPT subpath generation",
    "BDPT subpath connections",
    "SpatialLightDistribution lookup",
    "SpatialLightDistribution spin wait",
    "SpatialLightDistribution creation",
    "Direct lighting",
    "BSDF::f()",
    "BSDF::Sample_f()",
    "BSDF::PDF()",
    "BSSRDF::f()",
    "BSSRDF::Sample_f()",
    "PhaseFunction::p()",
    "PhaseFunction::Sample_p()",
    "Accelerator::Intersect()",
    "Accelerator::IntersectP()",
    "Light::Sample_*()",
    "Light::Pdf()",
    "Medium::Sample()",
    "Medium::Tr()",
    "Triangle::Intersect()",
    "Triangle::IntersectP()",
    "Curve::Intersect()",
    "Curve::IntersectP()",
    "Other Shape::Intersect()",
    "Other Shape::IntersectP()",
    "Material::ComputeScatteringFunctions()",
    "Camera::GenerateRay[Differential]()",
    "Film::MergeTile()",
    "Film::AddSplat()",
    "Film::AddSample()",
    "Sampler::StartPixelSample()",
    "Sampler::GetSample[12]D()",
    "MIPMap::Lookup() (trilinear)",
    "MIPMap::Lookup() (EWA)",
    "Ptex lookup",
};

static_assert((int)Prof::NumProfCategories ==
                  sizeof(ProfNames) / sizeof(ProfNames[0]),
              "ProfNames[] array and Prof enumerant have different "
              "numbers of entries!");

extern PBRT_THREAD_LOCAL uint64_t ProfilerState;
inline uint64_t CurrentProfilerState() { return ProfilerState; }

// ProfilePhase: 性能分析阶段类（RAII风格）。
// 构造时设置指定类别的分析标志位，析构时自动恢复，
// 用于精确测量代码段的执行时间。支持嵌套使用。
class ProfilePhase {
  public:
    // ProfilePhase Public Methods
    // ProfilePhase 公有方法
    ProfilePhase(Prof p) {
        categoryBit = ProfToBits(p);
        reset = (ProfilerState & categoryBit) == 0;
        ProfilerState |= categoryBit;
    }
    ~ProfilePhase() {
        if (reset) ProfilerState &= ~categoryBit;
    }
    ProfilePhase(const ProfilePhase &) = delete;
    ProfilePhase &operator=(const ProfilePhase &) = delete;

  private:
    // ProfilePhase Private Data
    // ProfilePhase 私有数据
    bool reset;        // 是否需要恢复之前的分析状态（防止嵌套覆盖）
    uint64_t categoryBit;  // 当前阶段对应的性能分析类别位掩码
};

// InitProfiler: 初始化性能分析器
void InitProfiler();
// SuspendProfiler: 暂停性能分析
void SuspendProfiler();
// ResumeProfiler: 恢复性能分析
void ResumeProfiler();
// ProfilerWorkerThreadInit: 初始化工作线程的性能分析状态
void ProfilerWorkerThreadInit();
// ReportProfilerResults: 输出性能分析结果到指定文件流
void ReportProfilerResults(FILE *dest);
// ClearProfiler: 清除性能分析数据
void ClearProfiler();
// CleanupProfiler: 清理性能分析器资源
void CleanupProfiler();
void SuspendProfiler();
void ResumeProfiler();
void ProfilerWorkerThreadInit();
void ReportProfilerResults(FILE *dest);
void ClearProfiler();
void CleanupProfiler();

// Statistics Macros
// 统计宏：提供便捷的统计变量声明和注册方式

// STAT_COUNTER: 声明一个计数器统计变量
#define STAT_COUNTER(title, var)                           \
    static PBRT_THREAD_LOCAL int64_t var;                  \
    static void STATS_FUNC##var(StatsAccumulator &accum) { \
        accum.ReportCounter(title, var);                   \
        var = 0;                                           \
    }                                                      \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)
// STAT_MEMORY_COUNTER: 声明一个内存计数器统计变量（如纹理内存占用）
#define STAT_MEMORY_COUNTER(title, var)                    \
    static PBRT_THREAD_LOCAL int64_t var;                  \
    static void STATS_FUNC##var(StatsAccumulator &accum) { \
        accum.ReportMemoryCounter(title, var);             \
        var = 0;                                           \
    }                                                      \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)

#ifndef PBRT_HAVE_CONSTEXPR
#define STATS_INT64_T_MIN LLONG_MAX
#define STATS_INT64_T_MAX INT64_MIN
#define STATS_DBL_T_MIN DBL_MAX
#define STATS_DBL_T_MAX -DBL_MAX
#else
#define STATS_INT64_T_MIN std::numeric_limits<int64_t>::max()
#define STATS_INT64_T_MAX std::numeric_limits<int64_t>::lowest()
#define STATS_DBL_T_MIN std::numeric_limits<double>::max()
#define STATS_DBL_T_MAX std::numeric_limits<double>::lowest()
#endif

// STAT_INT_DISTRIBUTION: 声明一个整数型分布统计变量
#define STAT_INT_DISTRIBUTION(title, var)                                  \
    static PBRT_THREAD_LOCAL int64_t var##sum;                             \
    static PBRT_THREAD_LOCAL int64_t var##count;                           \
    static PBRT_THREAD_LOCAL int64_t var##min = (STATS_INT64_T_MIN);       \
    static PBRT_THREAD_LOCAL int64_t var##max = (STATS_INT64_T_MAX);       \
    static void STATS_FUNC##var(StatsAccumulator &accum) {                 \
        accum.ReportIntDistribution(title, var##sum, var##count, var##min, \
                                    var##max);                             \
        var##sum = 0;                                                      \
        var##count = 0;                                                    \
        var##min = std::numeric_limits<int64_t>::max();                    \
        var##max = std::numeric_limits<int64_t>::lowest();                 \
    }                                                                      \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)

// STAT_FLOAT_DISTRIBUTION: 声明一个浮点型分布统计变量
#define STAT_FLOAT_DISTRIBUTION(title, var)                                  \
    static PBRT_THREAD_LOCAL double var##sum;                                \
    static PBRT_THREAD_LOCAL int64_t var##count;                             \
    static PBRT_THREAD_LOCAL double var##min = (STATS_DBL_T_MIN);            \
    static PBRT_THREAD_LOCAL double var##max = (STATS_DBL_T_MAX);            \
    static void STATS_FUNC##var(StatsAccumulator &accum) {                   \
        accum.ReportFloatDistribution(title, var##sum, var##count, var##min, \
                                      var##max);                             \
        var##sum = 0;                                                        \
        var##count = 0;                                                      \
        var##min = std::numeric_limits<double>::max();                       \
        var##max = std::numeric_limits<double>::lowest();                    \
    }                                                                        \
    static StatRegisterer STATS_REG##var(STATS_FUNC##var)

// ReportValue: 向分布统计变量报告一个值（自动更新总和、计数、最小值和最大值）
#define ReportValue(var, value)                                   \
    do {                                                          \
        var##sum += value;                                        \
        var##count += 1;                                          \
        var##min = std::min(var##min, decltype(var##min)(value)); \
        var##max = std::max(var##max, decltype(var##min)(value)); \
    } while (0)

// STAT_PERCENT: 声明一个百分比统计变量（分子/分母）
#define STAT_PERCENT(title, numVar, denomVar)                 \
    static PBRT_THREAD_LOCAL int64_t numVar, denomVar;        \
    static void STATS_FUNC##numVar(StatsAccumulator &accum) { \
        accum.ReportPercentage(title, numVar, denomVar);      \
        numVar = denomVar = 0;                                \
    }                                                         \
    static StatRegisterer STATS_REG##numVar(STATS_FUNC##numVar)

// STAT_RATIO: 声明一个比率统计变量（分子/分母）
#define STAT_RATIO(title, numVar, denomVar)                   \
    static PBRT_THREAD_LOCAL int64_t numVar, denomVar;        \
    static void STATS_FUNC##numVar(StatsAccumulator &accum) { \
        accum.ReportRatio(title, numVar, denomVar);           \
        numVar = denomVar = 0;                                \
    }                                                         \
    static StatRegisterer STATS_REG##numVar(STATS_FUNC##numVar)

}  // namespace pbrt

#endif  // PBRT_CORE_STATS_H
