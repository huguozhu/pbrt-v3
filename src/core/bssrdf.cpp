
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


// core/bssrdf.cpp*
//
// 此文件实现了次表面散射（BSSRDF）相关的核心功能。
// BSSRDF（Bidirectional Surface Scattering Reflectance Distribution Function）
// 描述了光在介质表面下进行次表面散射后重新从表面射出的过程。
// 包含菲涅尔矩计算、光束扩散 BSSRDF 计算、以及查表 BSSRDF 的采样与求值。
//

#include "bssrdf.h"
#include "interpolation.h"
#include "parallel.h"
#include "scene.h"

namespace pbrt {

// BSSRDF Utility Functions

// 计算菲涅尔一阶矩 FresnelMoment1(eta)。
// 根据折射率 eta 是否小于 1 选择不同的多项式拟合系数。
Float FresnelMoment1(Float eta) {
    Float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta,
          eta5 = eta4 * eta;
    if (eta < 1)
        return 0.45966f - 1.73965f * eta + 3.37668f * eta2 - 3.904945 * eta3 +
               2.49277f * eta4 - 0.68441f * eta5;
    else
        return -4.61686f + 11.1136f * eta - 10.4646f * eta2 + 5.11455f * eta3 -
               1.27198f * eta4 + 0.12746f * eta5;
}

// 计算菲涅尔二阶矩 FresnelMoment2(eta)。
// 与一阶矩类似，根据 eta < 1 或 eta >= 1 使用不同的多项式拟合公式。
Float FresnelMoment2(Float eta) {
    Float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta,
          eta5 = eta4 * eta;
    if (eta < 1) {
        return 0.27614f - 0.87350f * eta + 1.12077f * eta2 - 0.65095f * eta3 +
               0.07883f * eta4 + 0.04860f * eta5;
    } else {
        Float r_eta = 1 / eta, r_eta2 = r_eta * r_eta, r_eta3 = r_eta2 * r_eta;
        return -547.033f + 45.3087f * r_eta3 - 218.725f * r_eta2 +
               458.843f * r_eta + 404.557f * eta - 189.519f * eta2 +
               54.9327f * eta3 - 9.00603f * eta4 + 0.63942f * eta5;
    }
}

// 使用偶极子近似（dipole approximation）计算多重散射 BSSRDF 扩散剖面。
// 通过对点光源深度进行蒙特卡洛采样，累加偶极子辐射度和矢量辐照度贡献。
// sigma_s: 散射系数, sigma_a: 吸收系数, g: 各向异性参数,
// eta: 相对折射率, r: 表面距离。
Float BeamDiffusionMS(Float sigma_s, Float sigma_a, Float g, Float eta,
                      Float r) {
    const int nSamples = 100;
    Float Ed = 0;
    // Precompute information for dipole integrand

    // Compute reduced scattering coefficients $\sigmaps, \sigmapt$ and albedo
    // $\rhop$
    Float sigmap_s = sigma_s * (1 - g);
    Float sigmap_t = sigma_a + sigmap_s;
    Float rhop = sigmap_s / sigmap_t;

    // Compute non-classical diffusion coefficient $D_\roman{G}$ using
    // Equation (15.24)
    Float D_g = (2 * sigma_a + sigmap_s) / (3 * sigmap_t * sigmap_t);

    // Compute effective transport coefficient $\sigmatr$ based on $D_\roman{G}$
    Float sigma_tr = std::sqrt(sigma_a / D_g);

    // Determine linear extrapolation distance $\depthextrapolation$ using
    // Equation (15.28)
    Float fm1 = FresnelMoment1(eta), fm2 = FresnelMoment2(eta);
    Float ze = -2 * D_g * (1 + 3 * fm2) / (1 - 2 * fm1);

    // Determine exitance scale factors using Equations (15.31) and (15.32)
    Float cPhi = .25f * (1 - 2 * fm1), cE = .5f * (1 - 3 * fm2);
    for (int i = 0; i < nSamples; ++i) {
        // Sample real point source depth $\depthreal$
        Float zr = -std::log(1 - (i + .5f) / nSamples) / sigmap_t;

        // Evaluate dipole integrand $E_{\roman{d}}$ at $\depthreal$ and add to
        // _Ed_
        Float zv = -zr + 2 * ze;
        Float dr = std::sqrt(r * r + zr * zr), dv = std::sqrt(r * r + zv * zv);

        // Compute dipole fluence rate $\dipole(r)$ using Equation (15.27)
        Float phiD = Inv4Pi / D_g * (std::exp(-sigma_tr * dr) / dr -
                                     std::exp(-sigma_tr * dv) / dv);

        // Compute dipole vector irradiance $-\N{}\cdot\dipoleE(r)$ using
        // Equation (15.27)
        Float EDn = Inv4Pi * (zr * (1 + sigma_tr * dr) *
                                  std::exp(-sigma_tr * dr) / (dr * dr * dr) -
                              zv * (1 + sigma_tr * dv) *
                                  std::exp(-sigma_tr * dv) / (dv * dv * dv));

        // Add contribution from dipole for depth $\depthreal$ to _Ed_
        Float E = phiD * cPhi + EDn * cE;
        Float kappa = 1 - std::exp(-2 * sigmap_t * (dr + zr));
        Ed += kappa * rhop * rhop * E;
    }
    return Ed / nSamples;
}

// 计算单次散射 BSSRDF 扩散剖面。
// 对临界角以下的深度进行采样，累加单次散射对出射辐射度的贡献。
Float BeamDiffusionSS(Float sigma_s, Float sigma_a, Float g, Float eta,
                      Float r) {
    // Compute material parameters and minimum $t$ below the critical angle
    Float sigma_t = sigma_a + sigma_s, rho = sigma_s / sigma_t;
    Float tCrit = r * std::sqrt(eta * eta - 1);
    Float Ess = 0;
    const int nSamples = 100;
    for (int i = 0; i < nSamples; ++i) {
        // Evaluate single scattering integrand and add to _Ess_
        Float ti = tCrit - std::log(1 - (i + .5f) / nSamples) / sigma_t;

        // Determine length $d$ of connecting segment and $\cos\theta_\roman{o}$
        Float d = std::sqrt(r * r + ti * ti);
        Float cosThetaO = ti / d;

        // Add contribution of single scattering at depth $t$
        Ess += rho * std::exp(-sigma_t * (d + tCrit)) / (d * d) *
               PhaseHG(cosThetaO, g) * (1 - FrDielectric(-cosThetaO, 1, eta)) *
               std::abs(cosThetaO);
    }
    return Ess / nSamples;
}

// 计算光束扩散 BSSRDF 表。预处理不同半径和反照率下的扩散剖面值，
// 并计算有效反照率和累积分布函数（CDF）用于重要性采样。
void ComputeBeamDiffusionBSSRDF(Float g, Float eta, BSSRDFTable *t) {
    // Choose radius values of the diffusion profile discretization
    t->radiusSamples[0] = 0;
    t->radiusSamples[1] = 2.5e-3f;
    for (int i = 2; i < t->nRadiusSamples; ++i)
        t->radiusSamples[i] = t->radiusSamples[i - 1] * 1.2f;

    // Choose albedo values of the diffusion profile discretization
    for (int i = 0; i < t->nRhoSamples; ++i)
        t->rhoSamples[i] =
            (1 - std::exp(-8 * i / (Float)(t->nRhoSamples - 1))) /
            (1 - std::exp(-8));
    ParallelFor([&](int i) {
        // Compute the diffusion profile for the _i_th albedo sample

        // Compute scattering profile for chosen albedo $\rho$
        for (int j = 0; j < t->nRadiusSamples; ++j) {
            Float rho = t->rhoSamples[i], r = t->radiusSamples[j];
            t->profile[i * t->nRadiusSamples + j] =
                2 * Pi * r * (BeamDiffusionSS(rho, 1 - rho, g, eta, r) +
                              BeamDiffusionMS(rho, 1 - rho, g, eta, r));
        }

        // Compute effective albedo $\rho_{\roman{eff}}$ and CDF for importance
        // sampling
        t->rhoEff[i] =
            IntegrateCatmullRom(t->nRadiusSamples, t->radiusSamples.get(),
                                &t->profile[i * t->nRadiusSamples],
                                &t->profileCDF[i * t->nRadiusSamples]);
    }, t->nRhoSamples);
}

// 根据有效反照率 rhoEff 和平均自由程 mfp，反推吸收系数 sigma_a 和散射系数 sigma_s。
// 通过 Catmull-Rom 插值反函数将漫反射参数转换为次表面散射参数。
void SubsurfaceFromDiffuse(const BSSRDFTable &t, const Spectrum &rhoEff,
                           const Spectrum &mfp, Spectrum *sigma_a,
                           Spectrum *sigma_s) {
    for (int c = 0; c < Spectrum::nSamples; ++c) {
        Float rho = InvertCatmullRom(t.nRhoSamples, t.rhoSamples.get(),
                                     t.rhoEff.get(), rhoEff[c]);
        (*sigma_s)[c] = rho / mfp[c];
        (*sigma_a)[c] = (1 - rho) / mfp[c];
    }
}

// BSSRDF Method Definitions

// BSSRDFTable 构造函数。分配用于存储反照率样本、半径样本、
// 扩散剖面、有效反照率和剖面 CDF 的数组内存。
BSSRDFTable::BSSRDFTable(int nRhoSamples, int nRadiusSamples)
    : nRhoSamples(nRhoSamples),
      nRadiusSamples(nRadiusSamples),
      rhoSamples(new Float[nRhoSamples]),
      radiusSamples(new Float[nRadiusSamples]),
      profile(new Float[nRadiusSamples * nRhoSamples]),
      rhoEff(new Float[nRhoSamples]),
      profileCDF(new Float[nRadiusSamples * nRhoSamples]) {}

// 计算查表 BSSRDF 在表面距离 r 处的值 Sr。
// 将距离转换为光学半径，通过 Catmull-Rom 张量样条插值从表中获取扩散剖面值，
// 然后变换回世界空间单位。
Spectrum TabulatedBSSRDF::Sr(Float r) const {
    Spectrum Sr(0.f);
    for (int ch = 0; ch < Spectrum::nSamples; ++ch) {
        // Convert $r$ into unitless optical radius $r_{\roman{optical}}$
        Float rOptical = r * sigma_t[ch];

        // Compute spline weights to interpolate BSSRDF on channel _ch_
        int rhoOffset, radiusOffset;
        Float rhoWeights[4], radiusWeights[4];
        if (!CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(),
                               rho[ch], &rhoOffset, rhoWeights) ||
            !CatmullRomWeights(table.nRadiusSamples, table.radiusSamples.get(),
                               rOptical, &radiusOffset, radiusWeights))
            continue;

        // Set BSSRDF value _Sr[ch]_ using tensor spline interpolation
        Float sr = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                Float weight = rhoWeights[i] * radiusWeights[j];
                if (weight != 0)
                    sr += weight *
                          table.EvalProfile(rhoOffset + i, radiusOffset + j);
            }
        }

        // Cancel marginal PDF factor from tabulated BSSRDF profile
        if (rOptical != 0) sr /= 2 * Pi * rOptical;
        Sr[ch] = sr;
    }
    // Transform BSSRDF value into world space units
    Sr *= sigma_t * sigma_t;
    return Sr.Clamp();
}

// 对可分离 BSSRDF 进行采样：先采样空间项 Sp，然后在采样点处初始化 BSDF，
// 使得后续可以在该点进行完整的着色计算。
Spectrum SeparableBSSRDF::Sample_S(const Scene &scene, Float u1,
                                   const Point2f &u2, MemoryArena &arena,
                                   SurfaceInteraction *si, Float *pdf) const {
    ProfilePhase pp(Prof::BSSRDFSampling);
    Spectrum Sp = Sample_Sp(scene, u1, u2, arena, si, pdf);
    if (!Sp.IsBlack()) {
        // Initialize material model at sampled surface interaction
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);
        si->bsdf->Add(ARENA_ALLOC(arena, SeparableBSSRDFAdapter)(this));
        si->wo = Vector3f(si->shading.n);
    }
    return Sp;
}

// 对可分离 BSSRDF 的空间项 Sp 进行采样。
// 随机选择一个投影轴，选择光谱通道，在极坐标下对 BSSRDF 剖面进行采样，
// 然后沿光线与场景求交，随机选择一个交点作为采样结果。
Spectrum SeparableBSSRDF::Sample_Sp(const Scene &scene, Float u1,
                                    const Point2f &u2, MemoryArena &arena,
                                    SurfaceInteraction *pi, Float *pdf) const {
    ProfilePhase pp(Prof::BSSRDFEvaluation);
    // Choose projection axis for BSSRDF sampling
    Vector3f vx, vy, vz;
    if (u1 < .5f) {
        vx = ss;
        vy = ts;
        vz = Vector3f(ns);
        u1 *= 2;
    } else if (u1 < .75f) {
        // Prepare for sampling rays with respect to _ss_
        vx = ts;
        vy = Vector3f(ns);
        vz = ss;
        u1 = (u1 - .5f) * 4;
    } else {
        // Prepare for sampling rays with respect to _ts_
        vx = Vector3f(ns);
        vy = ss;
        vz = ts;
        u1 = (u1 - .75f) * 4;
    }

    // Choose spectral channel for BSSRDF sampling
    int ch = Clamp((int)(u1 * Spectrum::nSamples), 0, Spectrum::nSamples - 1);
    u1 = u1 * Spectrum::nSamples - ch;

    // Sample BSSRDF profile in polar coordinates
    Float r = Sample_Sr(ch, u2[0]);
    if (r < 0) return Spectrum(0.f);
    Float phi = 2 * Pi * u2[1];

    // Compute BSSRDF profile bounds and intersection height
    Float rMax = Sample_Sr(ch, 0.999f);
    if (r >= rMax) return Spectrum(0.f);
    Float l = 2 * std::sqrt(rMax * rMax - r * r);

    // Compute BSSRDF sampling ray segment
    Interaction base;
    base.p =
        po.p + r * (vx * std::cos(phi) + vy * std::sin(phi)) - l * vz * 0.5f;
    base.time = po.time;
    Point3f pTarget = base.p + l * vz;

    // Intersect BSSRDF sampling ray against the scene geometry

    // Declare _IntersectionChain_ and linked list
    struct IntersectionChain {
        SurfaceInteraction si;
        IntersectionChain *next = nullptr;
    };
    IntersectionChain *chain = ARENA_ALLOC(arena, IntersectionChain)();

    // Accumulate chain of intersections along ray
    IntersectionChain *ptr = chain;
    int nFound = 0;
    while (true) {
        Ray r = base.SpawnRayTo(pTarget);
        if (r.d == Vector3f(0, 0, 0) || !scene.Intersect(r, &ptr->si))
            break;

        base = ptr->si;
        // Append admissible intersection to _IntersectionChain_
        if (ptr->si.primitive->GetMaterial() == this->material) {
            IntersectionChain *next = ARENA_ALLOC(arena, IntersectionChain)();
            ptr->next = next;
            ptr = next;
            nFound++;
        }
    }

    // Randomly choose one of several intersections during BSSRDF sampling
    if (nFound == 0) return Spectrum(0.0f);
    int selected = Clamp((int)(u1 * nFound), 0, nFound - 1);
    while (selected-- > 0) chain = chain->next;
    *pi = chain->si;

    // Compute sample PDF and return the spatial BSSRDF term $\Sp$
    *pdf = this->Pdf_Sp(*pi) / nFound;
    return this->Sp(*pi);
}

// 计算空间项 Sp 的采样概率密度函数（PDF）。
// 将采样点坐标转换到局部坐标系，计算三个投影轴上的剖面半径，
// 累加所有轴和光谱通道的概率密度。
Float SeparableBSSRDF::Pdf_Sp(const SurfaceInteraction &pi) const {
    // Express $\pti-\pto$ and $\bold{n}_i$ with respect to local coordinates at
    // $\pto$
    Vector3f d = po.p - pi.p;
    Vector3f dLocal(Dot(ss, d), Dot(ts, d), Dot(ns, d));
    Normal3f nLocal(Dot(ss, pi.n), Dot(ts, pi.n), Dot(ns, pi.n));

    // Compute BSSRDF profile radius under projection along each axis
    Float rProj[3] = {std::sqrt(dLocal.y * dLocal.y + dLocal.z * dLocal.z),
                      std::sqrt(dLocal.z * dLocal.z + dLocal.x * dLocal.x),
                      std::sqrt(dLocal.x * dLocal.x + dLocal.y * dLocal.y)};

    // Return combined probability from all BSSRDF sampling strategies
    Float pdf = 0, axisProb[3] = {.25f, .25f, .5f};
    Float chProb = 1 / (Float)Spectrum::nSamples;
    for (int axis = 0; axis < 3; ++axis)
        for (int ch = 0; ch < Spectrum::nSamples; ++ch)
            pdf += Pdf_Sr(ch, rProj[axis]) * std::abs(nLocal[axis]) * chProb *
                   axisProb[axis];
    return pdf;
}

// 对查表 BSSRDF 的径向剖面进行采样。
// 使用 Catmull-Rom 二维采样从预计算表中采样半径值，并转换回世界空间单位。
// 返回 -1 表示采样失败。
Float TabulatedBSSRDF::Sample_Sr(int ch, Float u) const {
    if (sigma_t[ch] == 0) return -1;
    return SampleCatmullRom2D(table.nRhoSamples, table.nRadiusSamples,
                              table.rhoSamples.get(), table.radiusSamples.get(),
                              table.profile.get(), table.profileCDF.get(),
                              rho[ch], u) /
           sigma_t[ch];
}

// 计算查表 BSSRDF 径向剖面在半径 r 处的概率密度。
// 使用 Catmull-Rom 样条插值从表中获取剖面值，再除以其有效反照率进行归一化。
Float TabulatedBSSRDF::Pdf_Sr(int ch, Float r) const {
    // Convert $r$ into unitless optical radius $r_{\roman{optical}}$
    Float rOptical = r * sigma_t[ch];

    // Compute spline weights to interpolate BSSRDF density on channel _ch_
    int rhoOffset, radiusOffset;
    Float rhoWeights[4], radiusWeights[4];
    if (!CatmullRomWeights(table.nRhoSamples, table.rhoSamples.get(), rho[ch],
                           &rhoOffset, rhoWeights) ||
        !CatmullRomWeights(table.nRadiusSamples, table.radiusSamples.get(),
                           rOptical, &radiusOffset, radiusWeights))
        return 0.f;

    // Return BSSRDF profile density for channel _ch_
    Float sr = 0, rhoEff = 0;
    for (int i = 0; i < 4; ++i) {
        if (rhoWeights[i] == 0) continue;
        rhoEff += table.rhoEff[rhoOffset + i] * rhoWeights[i];
        for (int j = 0; j < 4; ++j) {
            if (radiusWeights[j] == 0) continue;
            sr += table.EvalProfile(rhoOffset + i, radiusOffset + j) *
                  rhoWeights[i] * radiusWeights[j];
        }
    }

    // Cancel marginal PDF factor from tabulated BSSRDF profile
    if (rOptical != 0) sr /= 2 * Pi * rOptical;
    return std::max((Float)0, sr * sigma_t[ch] * sigma_t[ch] / rhoEff);
}

}  // namespace pbrt
