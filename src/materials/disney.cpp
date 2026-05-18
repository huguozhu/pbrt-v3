
/*
    pbrt source code is Copyright(c) 1998-2017
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

/*

Implementation of the Disney BSDF with Subsurface Scattering, as described in:
http://blog.selfshadow.com/publications/s2015-shading-course/burley/s2015_pbs_disney_bsdf_notes.pdf.

That model is based on the Disney BRDF, described in:
https://disney-animation.s3.amazonaws.com/uploads/production/publication_asset/48/asset/s2012_pbs_disney_brdf_notes_v3.pdf

Many thanks for Brent Burley and Karl Li for answering many questions about
the details of the implementation.

The initial implementation of the BRDF was adapted from
https://github.com/wdas/brdf/blob/master/src/brdfs/disney.brdf, which is
licensed under a slightly-modified Apache 2.0 license.

*/

// materials/disney.cpp*
// 文件描述: 迪士尼材质BSDF/BSSRDF的实现。基于Burley等人的迪士尼BRDF/BSDF模型，
// 包含迪士尼漫反射、次表面散射、清漆、金属、头发光泽等多个BRDF lobe。
#include "materials/disney.h"
#include "bssrdf.h"
#include "interaction.h"
#include "paramset.h"
#include "reflection.h"
#include "stats.h"
#include "stringprint.h"
#include "texture.h"
#include "rng.h"

namespace pbrt {

// 辅助函数: 计算平方
inline Float sqr(Float x) { return x * x; }

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
// Schlick菲涅耳近似: R = R(0) + (1 - R(0)) (1 - cos theta)^5
// 其中R(0)是法线入射时的反射率
inline Float SchlickWeight(Float cosTheta) {
    Float m = Clamp(1 - cosTheta, 0, 1);
    return (m * m) * (m * m) * m;
}

inline Float FrSchlick(Float R0, Float cosTheta) {
    return Lerp(SchlickWeight(cosTheta), R0, 1);
}

inline Spectrum FrSchlick(const Spectrum &R0, Float cosTheta) {
    return Lerp(SchlickWeight(cosTheta), R0, Spectrum(1.));
}

// 对于介电质，法线入射反射率R(0) = (eta-1)^2/(eta+1)^2(假设从空气入射)
inline Float SchlickR0FromEta(Float eta) { return sqr(eta - 1) / sqr(eta + 1); }

///////////////////////////////////////////////////////////////////////////
// DisneyDiffuse - 迪士尼漫反射BRDF
// 基于Burley 2015论文，在漫反射中加入菲涅耳效应，
// 使反射在掠射角时逐渐减弱。

class DisneyDiffuse : public BxDF {
  public:
    DisneyDiffuse(const Spectrum &R)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
    Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
    std::string ToString() const;

  private:
    Spectrum R;
};

// 评估迪士尼漫反射BRDF: 在法线入射时值为1，掠射角时降至0.5
Spectrum DisneyDiffuse::f(const Vector3f &wo, const Vector3f &wi) const {
    Float Fo = SchlickWeight(AbsCosTheta(wo)),
          Fi = SchlickWeight(AbsCosTheta(wi));

    // 漫反射菲涅耳: 法线入射时1，掠射时0.5，见Burley 2015公式(4)
    return R * InvPi * (1 - Fo / 2) * (1 - Fi / 2);
}

std::string DisneyDiffuse::ToString() const {
    return StringPrintf("[ DisneyDiffuse R: %s ]", R.ToString().c_str());
}

///////////////////////////////////////////////////////////////////////////
// DisneyFakeSS

// "Fake"次表面散射lobe, 基于Hanrahan-Krueger BRDF对BSSRDF的近似
// 用于薄表面情况下(thin=true)模拟次表面散射效果
class DisneyFakeSS : public BxDF {
  public:
    DisneyFakeSS(const Spectrum &R, Float roughness)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
          R(R),
          roughness(roughness) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
    Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
    std::string ToString() const;

  private:
    Spectrum R;
    Float roughness;
};

// 评估"假"次表面散射: 通过粗糙度参数控制回反射的平坦程度
Spectrum DisneyFakeSS::f(const Vector3f &wo, const Vector3f &wi) const {
    Vector3f wh = wi + wo;
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
    wh = Normalize(wh);
    Float cosThetaD = Dot(wi, wh);

    // Fss90用于根据粗糙度"展平"回反射
    Float Fss90 = cosThetaD * cosThetaD * roughness;
    Float Fo = SchlickWeight(AbsCosTheta(wo)),
          Fi = SchlickWeight(AbsCosTheta(wi));
    Float Fss = Lerp(Fo, 1.0, Fss90) * Lerp(Fi, 1.0, Fss90);
    // 1.25 scale is used to (roughly) preserve albedo
    Float ss =
        1.25f * (Fss * (1 / (AbsCosTheta(wo) + AbsCosTheta(wi)) - .5f) + .5f);

    return R * InvPi * ss;
}

std::string DisneyFakeSS::ToString() const {
    return StringPrintf("[ DisneyFakeSS R: %s roughness: %f ]",
                        R.ToString().c_str(), roughness);
}

///////////////////////////////////////////////////////////////////////////
// DisneyRetro - 迪士尼回反射BRDF
// 模拟表面微观结构导致的回反射(光沿入射方向返回)效应

class DisneyRetro : public BxDF {
  public:
    DisneyRetro(const Spectrum &R, Float roughness)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
          R(R),
          roughness(roughness) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
    Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
    std::string ToString() const;

  private:
    Spectrum R;
    Float roughness;
};

Spectrum DisneyRetro::f(const Vector3f &wo, const Vector3f &wi) const {
    Vector3f wh = wi + wo;
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
    wh = Normalize(wh);
    Float cosThetaD = Dot(wi, wh);

    Float Fo = SchlickWeight(AbsCosTheta(wo)),
          Fi = SchlickWeight(AbsCosTheta(wi));
    Float Rr = 2 * roughness * cosThetaD * cosThetaD;

    // Burley 2015, eq (4).
    return R * InvPi * Rr * (Fo + Fi + Fo * Fi * (Rr - 1));
}

std::string DisneyRetro::ToString() const {
    return StringPrintf("[ DisneyRetro R: %s roughness: %f ]",
                        R.ToString().c_str(), roughness);
}

///////////////////////////////////////////////////////////////////////////
// DisneySheen - 迪士尼光泽BRDF
// 模拟织物等表面上的彩色光泽(绒毛)效果

class DisneySheen : public BxDF {
  public:
    DisneySheen(const Spectrum &R)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
    Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
    std::string ToString() const;

  private:
    Spectrum R;
};

Spectrum DisneySheen::f(const Vector3f &wo, const Vector3f &wi) const {
    Vector3f wh = wi + wo;
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
    wh = Normalize(wh);
    Float cosThetaD = Dot(wi, wh);

    return R * SchlickWeight(cosThetaD);
}

std::string DisneySheen::ToString() const {
    return StringPrintf("[ DisneySheen R: %s]", R.ToString().c_str());
}

///////////////////////////////////////////////////////////////////////////
// DisneyClearcoat - 迪士尼清漆BRDF
// 使用GTR1分布(比Trowbridge-Reitz更宽的尾部)和固定IOR=1.5的额外高光层

class DisneyClearcoat : public BxDF {
  public:
    DisneyClearcoat(Float weight, Float gloss)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
          weight(weight),
          gloss(gloss) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    Float weight, gloss;
};

// GTR1分布: 具有比Trowbridge-Reitz更宽尾部的微表面法线分布
inline Float GTR1(Float cosTheta, Float alpha) {
    Float alpha2 = alpha * alpha;
    return (alpha2 - 1) /
           (Pi * std::log(alpha2) * (1 + (alpha2 - 1) * cosTheta * cosTheta));
}

// Smith几何遮挡/阴影函数(使用GGX模型)
inline Float smithG_GGX(Float cosTheta, Float alpha) {
    Float alpha2 = alpha * alpha;
    Float cosTheta2 = cosTheta * cosTheta;
    return 1 / (cosTheta + sqrt(alpha2 + cosTheta2 - alpha2 * cosTheta2));
}

// 评估清漆BRDF: 使用GTR1分布、固定IOR=1.5的菲涅耳和alpha=0.25的几何项
Spectrum DisneyClearcoat::f(const Vector3f &wo, const Vector3f &wi) const {
    Vector3f wh = wi + wo;
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
    wh = Normalize(wh);

    // 清漆硬编码IOR=1.5 -> F0=0.04。使用GTR1分布(尾部比Trowbridge-Reitz/GTR2更肥)
    Float Dr = GTR1(AbsCosTheta(wh), gloss);
    Float Fr = FrSchlick(.04, Dot(wo, wh));
    // The geometric term always based on alpha = 0.25.
    Float Gr =
        smithG_GGX(AbsCosTheta(wo), .25) * smithG_GGX(AbsCosTheta(wi), .25);

    return weight * Gr * Fr * Dr / 4;
}

Spectrum DisneyClearcoat::Sample_f(const Vector3f &wo, Vector3f *wi,
                                   const Point2f &u, Float *pdf,
                                   BxDFType *sampledType) const {
    // TODO: double check all this: there still seem to be some very
    // occasional fireflies with clearcoat; presumably there is a bug
    // somewhere.
    if (wo.z == 0) return 0.;

    Float alpha2 = gloss * gloss;
    Float cosTheta = std::sqrt(
        std::max(Float(0), (1 - std::pow(alpha2, 1 - u[0])) / (1 - alpha2)));
    Float sinTheta = std::sqrt(std::max((Float)0, 1 - cosTheta * cosTheta));
    Float phi = 2 * Pi * u[1];
    Vector3f wh = SphericalDirection(sinTheta, cosTheta, phi);
    if (!SameHemisphere(wo, wh)) wh = -wh;

    *wi = Reflect(wo, wh);
    if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);

    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

Float DisneyClearcoat::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (!SameHemisphere(wo, wi)) return 0;

    Vector3f wh = wi + wo;
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return 0;
    wh = Normalize(wh);

    // The sampling routine samples wh exactly from the GTR1 distribution.
    // Thus, the final value of the PDF is just the value of the
    // distribution for wh converted to a mesure with respect to the
    // surface normal.
    Float Dr = GTR1(AbsCosTheta(wh), gloss);
    return Dr * AbsCosTheta(wh) / (4 * Dot(wo, wh));
}

std::string DisneyClearcoat::ToString() const {
    return StringPrintf("[ DisneyClearcoat weight: %f gloss: %f ]", weight,
                        gloss);
}

///////////////////////////////////////////////////////////////////////////
// DisneyFresnel - 迪士尼菲涅耳函数
// 高光分量的专用菲涅耳函数，在介电质菲涅耳和Schlick菲涅耳近似之间混合
// 通过metallic参数控制混合比例
class DisneyFresnel : public Fresnel {
  public:
    DisneyFresnel(const Spectrum &R0, Float metallic, Float eta)
        : R0(R0), metallic(metallic), eta(eta) {}
    Spectrum Evaluate(Float cosI) const {
        return Lerp(metallic, Spectrum(FrDielectric(cosI, 1, eta)),
                    FrSchlick(R0, cosI));
    }
    std::string ToString() const {
        return StringPrintf("[ DisneyFresnel R0: %s metallic: %f eta: %f ]",
                            R0.ToString().c_str(), metallic, eta);
    }

  private:
    const Spectrum R0;
    const Float metallic, eta;
};

///////////////////////////////////////////////////////////////////////////
// DisneyMicrofacetDistribution - 迪士尼微表面分布
// 继承Trowbridge-Reitz分布，使用可分离的Smith几何遮挡模型
public:
    DisneyMicrofacetDistribution(Float alphax, Float alphay)
        : TrowbridgeReitzDistribution(alphax, alphay) {}

    Float G(const Vector3f &wo, const Vector3f &wi) const {
        // Disney uses the separable masking-shadowing model.
        return G1(wo) * G1(wi);
    }
};

///////////////////////////////////////////////////////////////////////////
// DisneyBSSRDF

// 经验BSSRDF实现，基于Burley的"将迪士尼BRDF扩展到集成次表面散射的BSDF"
// 以及Christensen和Burley的"高效次表面散射的近似反射率剖面"
class DisneyBSSRDF : public SeparableBSSRDF {
  public:
    DisneyBSSRDF(const Spectrum &R, const Spectrum &d,
                 const SurfaceInteraction &po, Float eta,
                 const Material *material, TransportMode mode)
        // 0.2因子来自Brent Burley和Matt Chiang的个人通信
        : SeparableBSSRDF(po, eta, material, mode), R(R), d(0.2 * d) {}

    Spectrum S(const SurfaceInteraction &pi, const Vector3f &wi);
    Spectrum Sr(Float d) const;
    Float Sample_Sr(int ch, Float u) const;
    Float Pdf_Sr(int ch, Float r) const;

  private:
    Spectrum R, d;
};

// 重写BSSRDF::S()以访问完整的击点信息，从而根据表面法线方向进行调制
Spectrum DisneyBSSRDF::S(const SurfaceInteraction &pi, const Vector3f &wi) {
    ProfilePhase pp(Prof::BSSRDFEvaluation);
    // 根据两个表面法线的相对方向进行淡入淡出，以更好地处理表面凹陷
    // (细节来自Brent Burley的个人通信，未在课程笔记中发布)
    Vector3f a = Normalize(pi.p - po.p);
    Float fade = 1;
    Vector3f n = Vector3f(po.shading.n);
    Float cosTheta = Dot(a, n);
    if (cosTheta > 0) {
        // Point on or above surface plane
        Float sinTheta = std::sqrt(std::max(Float(0), 1 - cosTheta * cosTheta));
        Vector3f a2 = n * sinTheta - (a - n * cosTheta) * cosTheta / sinTheta;
        fade = std::max(Float(0), Dot(pi.shading.n, a2));
    }

    Float Fo = SchlickWeight(AbsCosTheta(po.wo)),
          Fi = SchlickWeight(AbsCosTheta(wi));
    return fade * (1 - Fo / 2) * (1 - Fi / 2) * Sp(pi) / Pi;
}

// 扩散剖面函数，来自Burley 2015公式(5): 双指数衰减
Spectrum DisneyBSSRDF::Sr(Float r) const {
    ProfilePhase pp(Prof::BSSRDFEvaluation);
    if (r < 1e-6f) r = 1e-6f;  // 避免在r==0处的奇异性
    return R * (Exp(-Spectrum(r) / d) + Exp(-Spectrum(r) / (3 * d))) /
           (8 * Pi * d * r);
}

Float DisneyBSSRDF::Sample_Sr(int ch, Float u) const {
    // Sr中实现的扩散剖面已归一化(极坐标积分=1)
    // 使用MIS在两个指数项之间进行采样
    //
    // int_0^2pi int_0^Infinity Sr(r) r dr dphi == 1.
    //
    // The CDF can be found in closed-form. It is:
    //
    // 1 - e^(-x/d) / 4 - (3 / 4) e^(-x / (3d)).
    //
    // Unfortunately, inverting the CDF requires solving a cubic, which
    // would be nice to sidestep. Therefore, following Christensen and
    // Burley's suggestion (section 6), we will sample from each of the two
    // exponential terms individually (which can be done directly) and then
    // compute an overall PDF using MIS.  There are a few details to work
    // through...
    //
    // For the first exponential term, we can find:
    // normalized PDF: e^(-r/d) / (2 Pi d r)
    // CDF: 1 - e^(-r/d)
    // sampling recipe: r = d log(1 / (1 - u))
    //
    // For the second:
    // PDF: e^(-r/(3d)) / (6 Pi d r)
    // CDF: 1 - e^(-r/(3d))
    // sampling: r = 3 d log(1 / (1 - u))
    //
    // The last question is what fraction of samples to use for each
    // technique.  The second exponential has 3x the contribution to the
    // final value as the first does, so therefore we'll take three samples
    // from that for every one sample we take from the first.
    if (u < .25f) {
        // 从第一个指数项采样(占1/4)
        u = std::min<Float>(u * 4, OneMinusEpsilon);  // renormalize to [0,1)
        return d[ch] * std::log(1 / (1 - u));
    } else {
        // 从第二个指数项采样(占3/4，权重更大)
        u = std::min<Float>((u - .25f) / .75f, OneMinusEpsilon);  // normalize to [0,1)
        return 3 * d[ch] * std::log(1 / (1 - u));
    }
}

Float DisneyBSSRDF::Pdf_Sr(int ch, Float r) const {
    if (r < 1e-6f) r = 1e-6f;  // 避免在r==0处的奇异性

    // 按照Sample_Sr()中的采样频率加权两个PDF
    return (.25f * std::exp(-r / d[ch]) / (2 * Pi * d[ch] * r) +
            .75f * std::exp(-r / (3 * d[ch])) / (6 * Pi * d[ch] * r));
}

///////////////////////////////////////////////////////////////////////////
// DisneyMaterial

// DisneyMaterial Method Definitions
// 计算散射函数: 整合所有迪士尼BSDF/BSSRDF组件
// 包括漫反射、次表面散射、回反射、光泽、清漆、微表面高光和透射
void DisneyMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                                MemoryArena &arena,
                                                TransportMode mode,
                                                bool allowMultipleLobes) const {
    // 如果存在凹凸贴图则执行凹凸映射
    if (bumpMap) Bump(bumpMap, si);

    // 评估迪士尼材质纹理参数并分配BRDF
    si->bsdf = ARENA_ALLOC(arena, BSDF)(*si);

    // 漫反射分量: 计算颜色和材质属性参数
    Spectrum c = color->Evaluate(*si).Clamp();
    Float metallicWeight = metallic->Evaluate(*si);
    Float e = eta->Evaluate(*si);
    Float strans = specTrans->Evaluate(*si);
    Float diffuseWeight = (1 - metallicWeight) * (1 - strans);
    Float dt = diffTrans->Evaluate(*si) /
               2;  // 0: all diffuse is reflected -> 1, transmitted
    Float rough = roughness->Evaluate(*si);
    Float lum = c.y();
    // 归一化亮度以分离色调和饱和度
    Spectrum Ctint = lum > 0 ? (c / lum) : Spectrum(1.);

    Float sheenWeight = sheen->Evaluate(*si);
    Spectrum Csheen;
    if (sheenWeight > 0) {
        Float stint = sheenTint->Evaluate(*si);
        Csheen = Lerp(stint, Spectrum(1.), Ctint);
    }

    // 构建漫反射组件(包含薄表面和次表面散射分支)
    if (diffuseWeight > 0) {
        if (thin) {
            Float flat = flatness->Evaluate(*si);
            // Blend between DisneyDiffuse and fake subsurface based on
            // flatness.  Additionally, weight using diffTrans.
            si->bsdf->Add(ARENA_ALLOC(arena, DisneyDiffuse)(
                diffuseWeight * (1 - flat) * (1 - dt) * c));
            si->bsdf->Add(ARENA_ALLOC(arena, DisneyFakeSS)(
                diffuseWeight * flat * (1 - dt) * c, rough));
        } else {
            Spectrum sd = scatterDistance->Evaluate(*si);
            if (sd.IsBlack())
                // No subsurface scattering; use regular (Fresnel modified)
                // diffuse.
                si->bsdf->Add(
                    ARENA_ALLOC(arena, DisneyDiffuse)(diffuseWeight * c));
            else {
                // Use a BSSRDF instead.
                si->bsdf->Add(ARENA_ALLOC(arena, SpecularTransmission)(
                    1.f, 1.f, e, mode));
                si->bssrdf = ARENA_ALLOC(arena, DisneyBSSRDF)(
                    c * diffuseWeight, sd, *si, e, this, mode);
            }
        }

        // 回反射分量
        si->bsdf->Add(
            ARENA_ALLOC(arena, DisneyRetro)(diffuseWeight * c, rough));

        // 光泽分量(如果启用)
        if (sheenWeight > 0)
            si->bsdf->Add(ARENA_ALLOC(arena, DisneySheen)(
                diffuseWeight * sheenWeight * Csheen));
    }

    // 创建金属和/或镜面透射的微表面分布
    // 支持各向异性: aspect控制椭圆率
    Float aspect = std::sqrt(1 - anisotropic->Evaluate(*si) * .9);
    Float ax = std::max(Float(.001), sqr(rough) / aspect);
    Float ay = std::max(Float(.001), sqr(rough) * aspect);
    MicrofacetDistribution *distrib =
        ARENA_ALLOC(arena, DisneyMicrofacetDistribution)(ax, ay);

    // 高光反射: Trowbridge-Reitz分布 + 改进的迪士尼菲涅耳函数
    Float specTint = specularTint->Evaluate(*si);
    Spectrum Cspec0 =
        Lerp(metallicWeight,
             SchlickR0FromEta(e) * Lerp(specTint, Spectrum(1.), Ctint), c);
    Fresnel *fresnel =
        ARENA_ALLOC(arena, DisneyFresnel)(Cspec0, metallicWeight, e);
    si->bsdf->Add(
        ARENA_ALLOC(arena, MicrofacetReflection)(Spectrum(1.), distrib, fresnel));

    // 清漆层: 独立的GTR1高光层
    Float cc = clearcoat->Evaluate(*si);
    if (cc > 0) {
        si->bsdf->Add(ARENA_ALLOC(arena, DisneyClearcoat)(
            cc, Lerp(clearcoatGloss->Evaluate(*si), .1, .001)));
    }

    // 透射分量(BTDF): 使用微表面透射模型处理镜面透射
    if (strans > 0) {
        // Walter等人模型，透射项使用sqrt(color)缩放，
        // 使得经过两次折射后回到指定的颜色
        Spectrum T = strans * Sqrt(c);
        if (thin) {
            // Scale roughness based on IOR (Burley 2015, Figure 15).
            Float rscaled = (0.65f * e - 0.35f) * rough;
            Float ax = std::max(Float(.001), sqr(rscaled) / aspect);
            Float ay = std::max(Float(.001), sqr(rscaled) * aspect);
            MicrofacetDistribution *scaledDistrib =
                ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(ax, ay);
            si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
                T, scaledDistrib, 1., e, mode));
        } else
            si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
                T, distrib, 1., e, mode));
    }
    if (thin) {
        // 薄表面: Lambertian透射，由diffTrans加权
        si->bsdf->Add(ARENA_ALLOC(arena, LambertianTransmission)(dt * c));
    }
}

// 创建迪士尼材质对象的工厂函数
DisneyMaterial *CreateDisneyMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> color =
        mp.GetSpectrumTexture("color", Spectrum(0.5f));
    std::shared_ptr<Texture<Float>> metallic =
        mp.GetFloatTexture("metallic", 0.f);
    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTexture("eta", 1.5f);
    std::shared_ptr<Texture<Float>> roughness =
        mp.GetFloatTexture("roughness", .5f);
    std::shared_ptr<Texture<Float>> specularTint =
        mp.GetFloatTexture("speculartint", 0.f);
    std::shared_ptr<Texture<Float>> anisotropic =
        mp.GetFloatTexture("anisotropic", 0.f);
    std::shared_ptr<Texture<Float>> sheen = mp.GetFloatTexture("sheen", 0.f);
    std::shared_ptr<Texture<Float>> sheenTint =
        mp.GetFloatTexture("sheentint", .5f);
    std::shared_ptr<Texture<Float>> clearcoat =
        mp.GetFloatTexture("clearcoat", 0.f);
    std::shared_ptr<Texture<Float>> clearcoatGloss =
        mp.GetFloatTexture("clearcoatgloss", 1.f);
    std::shared_ptr<Texture<Float>> specTrans =
        mp.GetFloatTexture("spectrans", 0.f);
    std::shared_ptr<Texture<Spectrum>> scatterDistance =
        mp.GetSpectrumTexture("scatterdistance", Spectrum(0.));
    bool thin = mp.FindBool("thin", false);
    std::shared_ptr<Texture<Float>> flatness =
        mp.GetFloatTexture("flatness", 0.f);
    std::shared_ptr<Texture<Float>> diffTrans =
        mp.GetFloatTexture("difftrans", 1.f);
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    return new DisneyMaterial(color, metallic, eta, roughness, specularTint,
                              anisotropic, sheen, sheenTint, clearcoat,
                              clearcoatGloss, specTrans, scatterDistance, thin,
                              flatness, diffTrans, bumpMap);
}

}  // namespace pbrt
