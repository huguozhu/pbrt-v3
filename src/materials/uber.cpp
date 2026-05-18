
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


// materials/uber.cpp*
// 文件描述: Uber(万能)材质的实现。综合了漫反射、微表面高光、
// 镜面反射、镜面透射和透明度控制。
#include "materials/uber.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "interaction.h"
#include "paramset.h"

namespace pbrt {

// UberMaterial Method Definitions
// 计算散射函数: 创建包含透明度、漫反射、高光、镜面反射和透射的综合BSDF
void UberMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                              MemoryArena &arena,
                                              TransportMode mode,
                                              bool allowMultipleLobes) const {
    // 如果存在凹凸贴图则执行凹凸映射
    if (bumpMap) Bump(bumpMap, si);
    Float e = eta->Evaluate(*si);

    Spectrum op = opacity->Evaluate(*si).Clamp();
    // 透明度处理: 创建SpecularTransmission用于透明效果
    Spectrum t = (-op + Spectrum(1.f)).Clamp();
    if (!t.IsBlack()) {
        // 有透明度时: 创建SpecularTransmission并设置折射率为1
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, 1.f);
        BxDF *tr = ARENA_ALLOC(arena, SpecularTransmission)(t, 1.f, 1.f, mode);
        si->bsdf->Add(tr);
    } else
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, e);

    // 漫反射分量: 使用不透明度加权的Kd
    Spectrum kd = op * Kd->Evaluate(*si).Clamp();
    if (!kd.IsBlack()) {
        // 添加Lambertian漫反射BSDF
        BxDF *diff = ARENA_ALLOC(arena, LambertianReflection)(kd);
        si->bsdf->Add(diff);
    }

    // 高光反射分量(Ks): 使用微表面模型(MicrofacetReflection)
    Spectrum ks = op * Ks->Evaluate(*si).Clamp();
    if (!ks.IsBlack()) {
        Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, e);
        Float roughu, roughv;
        if (roughnessu)
            roughu = roughnessu->Evaluate(*si);
        else
            roughu = roughness->Evaluate(*si);
        if (roughnessv)
            roughv = roughnessv->Evaluate(*si);
        else
            roughv = roughu;
        if (remapRoughness) {
            roughu = TrowbridgeReitzDistribution::RoughnessToAlpha(roughu);
            roughv = TrowbridgeReitzDistribution::RoughnessToAlpha(roughv);
        }
        MicrofacetDistribution *distrib =
            ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(roughu, roughv);
        BxDF *spec =
            ARENA_ALLOC(arena, MicrofacetReflection)(ks, distrib, fresnel);
        si->bsdf->Add(spec);
    }

    // 镜面反射分量(Kr): 使用SpecularReflection
    Spectrum kr = op * Kr->Evaluate(*si).Clamp();
    if (!kr.IsBlack()) {
        Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, e);
        si->bsdf->Add(ARENA_ALLOC(arena, SpecularReflection)(kr, fresnel));
    }

    // 镜面透射分量(Kt): 使用SpecularTransmission
    Spectrum kt = op * Kt->Evaluate(*si).Clamp();
    if (!kt.IsBlack())
        si->bsdf->Add(
            ARENA_ALLOC(arena, SpecularTransmission)(kt, 1.f, e, mode));
}

// 创建Uber材质对象的工厂函数
UberMaterial *CreateUberMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd =
        mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    std::shared_ptr<Texture<Spectrum>> Ks =
        mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    std::shared_ptr<Texture<Spectrum>> Kr =
        mp.GetSpectrumTexture("Kr", Spectrum(0.f));
    std::shared_ptr<Texture<Spectrum>> Kt =
        mp.GetSpectrumTexture("Kt", Spectrum(0.f));
    std::shared_ptr<Texture<Float>> roughness =
        mp.GetFloatTexture("roughness", .1f);
    std::shared_ptr<Texture<Float>> uroughness =
        mp.GetFloatTextureOrNull("uroughness");
    std::shared_ptr<Texture<Float>> vroughness =
        mp.GetFloatTextureOrNull("vroughness");
    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTextureOrNull("eta");
    if (!eta) eta = mp.GetFloatTexture("index", 1.5f);
    std::shared_ptr<Texture<Spectrum>> opacity =
        mp.GetSpectrumTexture("opacity", 1.f);
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);
    return new UberMaterial(Kd, Ks, Kr, Kt, roughness, uroughness, vroughness,
                            opacity, eta, bumpMap, remapRoughness);
}

}  // namespace pbrt
