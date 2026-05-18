
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


// media/grid.cpp*
// 文件描述: 网格密度体积介质的实现。使用三维密度网格和delta-tracking
// 进行光线-体积交互采样。
#include "media/grid.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"
#include "interaction.h"

namespace pbrt {

STAT_RATIO("Media/Grid steps per Tr() call", nTrSteps, nTrCalls);

// GridDensityMedium Method Definitions
// 计算给定位置p处的密度值(使用三线性插值从体素网格中采样)
Float GridDensityMedium::Density(const Point3f &p) const {
    // 计算体素坐标和偏移量(将p从[0,1]映射到体素空间)
    Point3f pSamples(p.x * nx - .5f, p.y * ny - .5f, p.z * nz - .5f);
    Point3i pi = (Point3i)Floor(pSamples);
    Vector3f d = pSamples - (Point3f)pi;

    // 三线性插值计算局部密度值
    Float d00 = Lerp(d.x, D(pi), D(pi + Vector3i(1, 0, 0)));
    Float d10 = Lerp(d.x, D(pi + Vector3i(0, 1, 0)), D(pi + Vector3i(1, 1, 0)));
    Float d01 = Lerp(d.x, D(pi + Vector3i(0, 0, 1)), D(pi + Vector3i(1, 0, 1)));
    Float d11 = Lerp(d.x, D(pi + Vector3i(0, 1, 1)), D(pi + Vector3i(1, 1, 1)));
    Float d0 = Lerp(d.y, d00, d10);
    Float d1 = Lerp(d.y, d01, d11);
    return Lerp(d.z, d0, d1);
}

// 使用delta-tracking(Woodcock追踪)方法采样介质中的散射事件。
// 将光线从世界空间变换到介质空间，然后迭代步进寻找散射位置。
Spectrum GridDensityMedium::Sample(const Ray &rWorld, Sampler &sampler,
                                   MemoryArena &arena,
                                   MediumInteraction *mi) const {
    ProfilePhase _(Prof::MediumSample);
    Ray ray = WorldToMedium(
        Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
    // 计算光线与介质包围盒的交叠区间[tMin, tMax]
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

    // 执行delta-tracking迭代来采样介质交互
    // 基本思想: 使用一个常数"主要波长"(sigma_t/invMaxDensity)生成试探步长，
    // 然后以实际密度与最大密度之比作为接受概率进行随机拒绝。
    Float t = tMin;
    while (true) {
        t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
        if (t >= tMax) break;
        if (Density(ray(t)) * invMaxDensity > sampler.Get1D()) {
            // 填充介质交互信息并返回
            PhaseFunction *phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
            *mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this,
                                    phase);
            return sigma_s / sigma_t;
        }
    }
    return Spectrum(1.f);
}

// 使用比率追踪方法计算光线穿过网格介质的透射率
Spectrum GridDensityMedium::Tr(const Ray &rWorld, Sampler &sampler) const {
    ProfilePhase _(Prof::MediumTr);
    ++nTrCalls;

    Ray ray = WorldToMedium(
        Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
    // 计算光线与介质包围盒的交叠区间[tMin, tMax]
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

    // 执行比率追踪来估计透射率值
    Float Tr = 1, t = tMin;
    while (true) {
        ++nTrSteps;
        t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
        if (t >= tMax) break;
        Float density = Density(ray(t));
        Tr *= 1 - std::max((Float)0, density * invMaxDensity);
        // 当透射率变低时，开始应用俄罗斯轮盘赌来提前终止采样
        const Float rrThreshold = .1;
        if (Tr < rrThreshold) {
            Float q = std::max((Float).05, 1 - Tr);
            if (sampler.Get1D() < q) return 0;
            Tr /= 1 - q;
        }
    }
    return Spectrum(Tr);
}

}  // namespace pbrt
