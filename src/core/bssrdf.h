
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

#ifndef PBRT_CORE_BSSRDF_H
#define PBRT_CORE_BSSRDF_H

// core/bssrdf.h*
// BSSRDF: 次表面散射(Bidirectional Surface Scattering Reflectance Distribution Function)模型，
// 用于模拟光在半透明材质（如皮肤、玉石、牛奶）内部的次表面散射效果
#include "interaction.h"
#include "reflection.h"
#include "stats.h"

namespace pbrt {

// BSSRDF Utility Declarations
// 菲涅尔矩计算，用于次表面散射的漫反射透射率近似
Float FresnelMoment1(Float invEta);
Float FresnelMoment2(Float invEta);

// BSSRDF Declarations
// BSSRDF基类: 定义次表面散射的通用接口
class BSSRDF {
  public:
    // BSSRDF Public Methods
    BSSRDF(const SurfaceInteraction &po, Float eta) : po(po), eta(eta) {}
    virtual ~BSSRDF() {}

    // BSSRDF Interface
    // S: 计算次表面散射分布函数值，pi为出射点，wi为出射方向
    virtual Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi) = 0;
    // Sample_S: 对次表面散射进行重要性采样，返回采样结果和概率密度
    virtual Spectrum Sample_S(const Scene &scene, Float u1, const Point2f &u2,
                              MemoryArena &arena, SurfaceInteraction *si,
                              Float *pdf) const = 0;

  protected:
    // BSSRDF Protected Data
    const SurfaceInteraction &po;  // 入射点（光从介质进入表面的位置）
    Float eta;                     // 相对折射率
};

// SeparableBSSRDF: 可分离的BSSRDF，将次表面散射分解为空间函数Sp和方向函数Sw
class SeparableBSSRDF : public BSSRDF {

class SeparableBSSRDF : public BSSRDF {
    friend class SeparableBSSRDFAdapter;

  public:
    // SeparableBSSRDF Public Methods
    SeparableBSSRDF(const SurfaceInteraction &po, Float eta,
                    const Material *material, TransportMode mode)
        : BSSRDF(po, eta),
          ns(po.shading.n),
          ss(Normalize(po.shading.dpdu)),
          ts(Cross(ns, ss)),
          material(material),
          mode(mode) {}
    Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi) {
        ProfilePhase pp(Prof::BSSRDFEvaluation);
        Float Ft = FrDielectric(CosTheta(po.wo), 1, eta);
        return (1 - Ft) * Sp(pi) * Sw(wi);
    }
    Spectrum Sw(const Vector3f &w) const {
        Float c = 1 - 2 * FresnelMoment1(1 / eta);
        return (1 - FrDielectric(CosTheta(w), 1, eta)) / (c * Pi);
    }
    Spectrum Sp(const SurfaceInteraction &pi) const {
        return Sr(Distance(po.p, pi.p));
    }
    Spectrum Sample_S(const Scene &scene, Float u1, const Point2f &u2,
                      MemoryArena &arena, SurfaceInteraction *si,
                      Float *pdf) const;
    Spectrum Sample_Sp(const Scene &scene, Float u1, const Point2f &u2,
                       MemoryArena &arena, SurfaceInteraction *si,
                       Float *pdf) const;
    Float Pdf_Sp(const SurfaceInteraction &si) const;

    // SeparableBSSRDF Interface
    virtual Spectrum Sr(Float d) const = 0;
    virtual Float Sample_Sr(int ch, Float u) const = 0;
    virtual Float Pdf_Sr(int ch, Float r) const = 0;

  private:
    // SeparableBSSRDF Private Data
    const Normal3f ns;           // 法线方向（着色空间）
    const Vector3f ss, ts;       // 切线和副法线方向（着色空间）
    const Material *material;    // 关联的材质
    const TransportMode mode;    // 传输模式（辐射度或重要性）
};

// TabulatedBSSRDF: 使用预计算表格的BSSRDF模型，基于dipole或多极子模型近似
class TabulatedBSSRDF : public SeparableBSSRDF {
  public:
    // TabulatedBSSRDF Public Methods
    TabulatedBSSRDF(const SurfaceInteraction &po, const Material *material,
                    TransportMode mode, Float eta, const Spectrum &sigma_a,
                    const Spectrum &sigma_s, const BSSRDFTable &table)
        : SeparableBSSRDF(po, eta, material, mode), table(table) {
        sigma_t = sigma_a + sigma_s;
        for (int c = 0; c < Spectrum::nSamples; ++c)
            rho[c] = sigma_t[c] != 0 ? (sigma_s[c] / sigma_t[c]) : 0;
    }
    // Sr: 根据空间距离计算散射分布值
    Spectrum Sr(Float distance) const;
    // Pdf_Sr: 计算空间采样的概率密度
    Float Pdf_Sr(int ch, Float distance) const;
    // Sample_Sr: 对空间散射距离进行重要性采样
    Float Sample_Sr(int ch, Float sample) const;

  private:
    // TabulatedBSSRDF Private Data
    const BSSRDFTable &table;  // 预计算表格引用
    Spectrum sigma_t, rho;     // 消光系数和反照率
};

// BSSRDFTable: 存储BSSRDF的预计算表格数据，包含漫反射轮廓和CDF

struct BSSRDFTable {
    // BSSRDFTable Public Data
    const int nRhoSamples, nRadiusSamples;       // 反照率和半径采样维度
    std::unique_ptr<Float[]> rhoSamples;           // 反照率采样点
    std::unique_ptr<Float[]> radiusSamples;        // 半径采样点
    std::unique_ptr<Float[]> profile;              // BSSRDF轮廓数据
    std::unique_ptr<Float[]> rhoEff;               // 有效反照率
    std::unique_ptr<Float[]> profileCDF;           // 轮廓CDF用于重要性采样

    // BSSRDFTable Public Methods
    BSSRDFTable(int nRhoSamples, int nRadiusSamples);
    // EvalProfile: 在给定索引处评估BSSRDF轮廓值
    inline Float EvalProfile(int rhoIndex, int radiusIndex) const {
        return profile[rhoIndex * nRadiusSamples + radiusIndex];
    }
};

// SeparableBSSRDFAdapter: 将SeparableBSSRDF适配为BxDF接口，用于BSDF中处理次表面散射的方向部分
class SeparableBSSRDFAdapter : public BxDF {
  public:
    // SeparableBSSRDFAdapter Public Methods
    SeparableBSSRDFAdapter(const SeparableBSSRDF *bssrdf)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), bssrdf(bssrdf) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        Spectrum f = bssrdf->Sw(wi);
        // Update BSSRDF transmission term to account for adjoint light
        // transport
        if (bssrdf->mode == TransportMode::Radiance)
            f *= bssrdf->eta * bssrdf->eta;
        return f;
    }
    std::string ToString() const { return "[ SeparableBSSRDFAdapter ]"; }

  private:
    const SeparableBSSRDF *bssrdf;
};

// 光束扩散函数：单次散射和多次散射近似
Float BeamDiffusionSS(Float sigma_s, Float sigma_a, Float g, Float eta,
                      Float r);
Float BeamDiffusionMS(Float sigma_s, Float sigma_a, Float g, Float eta,
                      Float r);
// 根据漫反射参数计算BSSRDF表格
void ComputeBeamDiffusionBSSRDF(Float g, Float eta, BSSRDFTable *t);
// 从漫反射颜色反推材质吸收和散射系数
void SubsurfaceFromDiffuse(const BSSRDFTable &table, const Spectrum &rhoEff,
                           const Spectrum &mfp, Spectrum *sigma_a,
                           Spectrum *sigma_s);

}  // namespace pbrt

#endif  // PBRT_CORE_BSSRDF_H
