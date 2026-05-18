
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

// core/microfacet.cpp*
//
// 此文件实现了微表面分布模型（Microfacet Distribution），包括 Beckmann 和
// Trowbridge-Reitz（GGX）两种分布。微表面理论将表面建模为大量微小镜面微面
// 的集合，每个微面遵循几何光学。该模块提供了法线分布函数 D()、遮蔽函数
// Lambda() 以及微面法线采样等功能，用于基于物理的渲染中的表面反射建模。
//
#include "microfacet.h"
#include "reflection.h"

namespace pbrt {

// Microfacet Utility Functions
// 微表面工具函数：Beckmann 分布的单变量采样（第11类）
// 根据入射角 cosThetaI 和随机数 U1/U2 计算微面法线的斜率分量
static void BeckmannSample11(Float cosThetaI, Float U1, Float U2,
                             Float *slope_x, Float *slope_y) {
    /* Special case (normal incidence) */
    if (cosThetaI > .9999) {
        Float r = std::sqrt(-std::log(1.0f - U1));
        Float sinPhi = std::sin(2 * Pi * U2);
        Float cosPhi = std::cos(2 * Pi * U2);
        *slope_x = r * cosPhi;
        *slope_y = r * sinPhi;
        return;
    }

    /* The original inversion routine from the paper contained
       discontinuities, which causes issues for QMC integration
       and techniques like Kelemen-style MLT. The following code
       performs a numerical inversion with better behavior */
    Float sinThetaI =
        std::sqrt(std::max((Float)0, (Float)1 - cosThetaI * cosThetaI));
    Float tanThetaI = sinThetaI / cosThetaI;
    Float cotThetaI = 1 / tanThetaI;

    /* Search interval -- everything is parameterized
       in the Erf() domain */
    Float a = -1, c = Erf(cotThetaI);
    Float sample_x = std::max(U1, (Float)1e-6f);

    /* Start with a good initial guess */
    // Float b = (1-sample_x) * a + sample_x * c;

    /* We can do better (inverse of an approximation computed in
     * Mathematica) */
    Float thetaI = std::acos(cosThetaI);
    Float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
    Float b = c - (1 + c) * std::pow(1 - sample_x, fit);

    /* Normalization factor for the CDF */
    static const Float SQRT_PI_INV = 1.f / std::sqrt(Pi);
    Float normalization =
        1 /
        (1 + c + SQRT_PI_INV * tanThetaI * std::exp(-cotThetaI * cotThetaI));

    int it = 0;
    while (++it < 10) {
        /* Bisection criterion -- the oddly-looking
           Boolean expression are intentional to check
           for NaNs at little additional cost */
        if (!(b >= a && b <= c)) b = 0.5f * (a + c);

        /* Evaluate the CDF and its derivative
           (i.e. the density function) */
        Float invErf = ErfInv(b);
        Float value =
            normalization *
                (1 + b + SQRT_PI_INV * tanThetaI * std::exp(-invErf * invErf)) -
            sample_x;
        Float derivative = normalization * (1 - invErf * tanThetaI);

        if (std::abs(value) < 1e-5f) break;

        /* Update bisection intervals */
        if (value > 0)
            c = b;
        else
            a = b;

        b -= value / derivative;
    }

    /* Now convert back into a slope value */
    *slope_x = ErfInv(b);

    /* Simulate Y component */
    *slope_y = ErfInv(2.0f * std::max(U2, (Float)1e-6f) - 1.0f);

    CHECK(!std::isinf(*slope_x));
    CHECK(!std::isnan(*slope_x));
    CHECK(!std::isinf(*slope_y));
    CHECK(!std::isnan(*slope_y));
}

// Beckmann 分布的全三维采样：给定入射方向 wi 和粗糙度参数 alpha_x/alpha_y，
// 按照可见法线分布采样微面法线方向。
static Vector3f BeckmannSample(const Vector3f &wi, Float alpha_x, Float alpha_y,
                               Float U1, Float U2) {
    // 1. stretch wi
    Vector3f wiStretched =
        Normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    Float slope_x, slope_y;
    BeckmannSample11(CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);

    // 3. rotate
    Float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
    slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. compute normal
    return Normalize(Vector3f(-slope_x, -slope_y, 1.f));
}

// MicrofacetDistribution Method Definitions
// 微表面分布基类析构函数
MicrofacetDistribution::~MicrofacetDistribution() {}

// Beckmann 分布的法线分布函数 D(wh)
// 计算给定微面法线 wh 下的 Beckmann 分布概率密度值
Float BeckmannDistribution::D(const Vector3f &wh) const {
    Float tan2Theta = Tan2Theta(wh);
    if (std::isinf(tan2Theta)) return 0.;
    Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
    return std::exp(-tan2Theta * (Cos2Phi(wh) / (alphax * alphax) +
                                  Sin2Phi(wh) / (alphay * alphay))) /
           (Pi * alphax * alphay * cos4Theta);
}

// Trowbridge-Reitz (GGX) 分布的法线分布函数 D(wh)
// 计算给定微面法线 wh 下的 GGX 分布概率密度值
Float TrowbridgeReitzDistribution::D(const Vector3f &wh) const {
    Float tan2Theta = Tan2Theta(wh);
    if (std::isinf(tan2Theta)) return 0.;
    const Float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
    Float e =
        (Cos2Phi(wh) / (alphax * alphax) + Sin2Phi(wh) / (alphay * alphay)) *
        tan2Theta;
    return 1 / (Pi * alphax * alphay * cos4Theta * (1 + e) * (1 + e));
}

// Beckmann 分布的遮蔽函数 Lambda(w)
// 计算给定方向 w 下未被遮蔽的微面比例，用于 Smith 遮蔽模型
Float BeckmannDistribution::Lambda(const Vector3f &w) const {
    Float absTanTheta = std::abs(TanTheta(w));
    if (std::isinf(absTanTheta)) return 0.;
    // Compute _alpha_ for direction _w_
    Float alpha =
        std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
    Float a = 1 / (alpha * absTanTheta);
    if (a >= 1.6f) return 0;
    return (1 - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a);
}

// Trowbridge-Reitz (GGX) 分布的遮蔽函数 Lambda(w)
// 计算给定方向 w 下未被遮蔽的微面比例，用于 Smith 遮蔽模型
Float TrowbridgeReitzDistribution::Lambda(const Vector3f &w) const {
    Float absTanTheta = std::abs(TanTheta(w));
    if (std::isinf(absTanTheta)) return 0.;
    // Compute _alpha_ for direction _w_
    Float alpha =
        std::sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
    Float alpha2Tan2Theta = (alpha * absTanTheta) * (alpha * absTanTheta);
    return (-1 + std::sqrt(1.f + alpha2Tan2Theta)) / 2;
}

// 将 Beckmann 分布参数转换为字符串表示
std::string BeckmannDistribution::ToString() const {
    return StringPrintf("[ BeckmannDistribution alphax: %f alphay: %f ]",
                        alphax, alphay);
}

// 将 Trowbridge-Reitz 分布参数转换为字符串表示
std::string TrowbridgeReitzDistribution::ToString() const {
    return StringPrintf("[ TrowbridgeReitzDistribution alphax: %f alphay: %f ]",
                        alphax, alphay);
}

// Beckmann 分布的微面法线采样函数
// 根据入射方向 wo 和随机样本 u，采样微面法线方向 wh
Vector3f BeckmannDistribution::Sample_wh(const Vector3f &wo,
                                         const Point2f &u) const {
    if (!sampleVisibleArea) {
        // Sample full distribution of normals for Beckmann distribution
        // 采样完整的法线分布（不可见区域采样）

        // Compute $\tan^2 \theta$ and $\phi$ for Beckmann distribution sample
        // 计算 Beckmann 分布采样的 tan^2(theta) 和 phi
        Float tan2Theta, phi;
        if (alphax == alphay) {
            // 各向同性情况：直接通过对数映射采样
            Float logSample = std::log(1 - u[0]);
            DCHECK(!std::isinf(logSample));
            tan2Theta = -alphax * alphax * logSample;
            phi = u[1] * 2 * Pi;
        } else {
            // Compute _tan2Theta_ and _phi_ for anisotropic Beckmann
            // distribution
            // 各向异性情况：根据粗糙度比率计算 phi 和 tan2Theta
            Float logSample = std::log(1 - u[0]);
            DCHECK(!std::isinf(logSample));
            phi = std::atan(alphay / alphax *
                            std::tan(2 * Pi * u[1] + 0.5f * Pi));
            if (u[1] > 0.5f) phi += Pi;
            Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
            Float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
            tan2Theta = -logSample /
                        (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
        }

        // Map sampled Beckmann angles to normal direction _wh_
        // 将采样的 Beckmann 角度映射为法线方向 wh
        Float cosTheta = 1 / std::sqrt(1 + tan2Theta);
        Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
        Vector3f wh = SphericalDirection(sinTheta, cosTheta, phi);
        if (!SameHemisphere(wo, wh)) wh = -wh;
        return wh;
    } else {
        // Sample visible area of normals for Beckmann distribution
        // 采样可见法线区域（只对可见的微面进行采样，减少方差）
        Vector3f wh;
        bool flip = wo.z < 0;
        wh = BeckmannSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
        if (flip) wh = -wh;
        return wh;
    }
}

// Trowbridge-Reitz 分布的单变量采样（第11类）
// 根据入射角余弦值和随机数 U1/U2 计算微面法线的斜率分量
static void TrowbridgeReitzSample11(Float cosTheta, Float U1, Float U2,
                                    Float *slope_x, Float *slope_y) {
    // special case (normal incidence)
    // 特殊情况：垂直入射
    if (cosTheta > .9999) {
        Float r = sqrt(U1 / (1 - U1));
        Float phi = 6.28318530718 * U2;
        *slope_x = r * cos(phi);
        *slope_y = r * sin(phi);
        return;
    }

    Float sinTheta =
        std::sqrt(std::max((Float)0, (Float)1 - cosTheta * cosTheta));
    Float tanTheta = sinTheta / cosTheta;
    Float a = 1 / tanTheta;
    Float G1 = 2 / (1 + std::sqrt(1.f + 1.f / (a * a)));

    // sample slope_x
    // 采样 x 方向斜率：基于 Smith 遮蔽函数的逆变换采样
    Float A = 2 * U1 / G1 - 1;
    Float tmp = 1.f / (A * A - 1.f);
    if (tmp > 1e10) tmp = 1e10;
    Float B = tanTheta;
    Float D = std::sqrt(
        std::max(Float(B * B * tmp * tmp - (A * A - B * B) * tmp), Float(0)));
    Float slope_x_1 = B * tmp - D;
    Float slope_x_2 = B * tmp + D;
    *slope_x = (A < 0 || slope_x_2 > 1.f / tanTheta) ? slope_x_1 : slope_x_2;

    // sample slope_y
    // 采样 y 方向斜率：使用有理函数拟合分布
    Float S;
    if (U2 > 0.5f) {
        S = 1.f;
        U2 = 2.f * (U2 - .5f);
    } else {
        S = -1.f;
        U2 = 2.f * (.5f - U2);
    }
    Float z =
        (U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
        (U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
    *slope_y = S * z * std::sqrt(1.f + *slope_x * *slope_x);

    CHECK(!std::isinf(*slope_y));
    CHECK(!std::isnan(*slope_y));
}

// Trowbridge-Reitz 分布的全三维采样：给定入射方向 wi 和粗糙度参数，
// 按照可见法线分布采样微面法线方向。
static Vector3f TrowbridgeReitzSample(const Vector3f &wi, Float alpha_x,
                                      Float alpha_y, Float U1, Float U2) {
    // 1. stretch wi
    // 第一步：将入射方向按粗糙度拉伸
    Vector3f wiStretched =
        Normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    Float slope_x, slope_y;
    TrowbridgeReitzSample11(CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);

    // 3. rotate
    Float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
    slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. compute normal
    return Normalize(Vector3f(-slope_x, -slope_y, 1.));
}

// Trowbridge-Reitz 分布的微面法线采样函数
// 根据入射方向 wo 和随机样本 u，采样微面法线方向 wh
Vector3f TrowbridgeReitzDistribution::Sample_wh(const Vector3f &wo,
                                                const Point2f &u) const {
    Vector3f wh;
    if (!sampleVisibleArea) {
        // 采样完整的法线分布（不可见区域采样）
        Float cosTheta = 0, phi = (2 * Pi) * u[1];
        if (alphax == alphay) {
            Float tanTheta2 = alphax * alphax * u[0] / (1.0f - u[0]);
            cosTheta = 1 / std::sqrt(1 + tanTheta2);
        } else {
            phi =
                std::atan(alphay / alphax * std::tan(2 * Pi * u[1] + .5f * Pi));
            if (u[1] > .5f) phi += Pi;
            Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
            const Float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
            const Float alpha2 =
                1 / (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
            Float tanTheta2 = alpha2 * u[0] / (1 - u[0]);
            cosTheta = 1 / std::sqrt(1 + tanTheta2);
        }
        Float sinTheta =
            std::sqrt(std::max((Float)0., (Float)1. - cosTheta * cosTheta));
        wh = SphericalDirection(sinTheta, cosTheta, phi);
        if (!SameHemisphere(wo, wh)) wh = -wh;
    } else {
        // 采样可见法线区域
        bool flip = wo.z < 0;
        wh = TrowbridgeReitzSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
        if (flip) wh = -wh;
    }
    return wh;
}

// 微表面分布的概率密度函数 Pdf
// 根据是否使用可见法线采样，返回对应的概率密度值
Float MicrofacetDistribution::Pdf(const Vector3f &wo,
                                  const Vector3f &wh) const {
    if (sampleVisibleArea)
        // 可见法线采样：包含遮蔽函数和投影项
        return D(wh) * G1(wo) * AbsDot(wo, wh) / AbsCosTheta(wo);
    else
        // 不可见法线采样：仅基于法线分布和余弦项
        return D(wh) * AbsCosTheta(wh);
}

}  // namespace pbrt
