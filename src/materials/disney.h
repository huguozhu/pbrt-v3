
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_DISNEY_H
#define PBRT_MATERIALS_DISNEY_H

// materials/disney.h*
// 文件描述: 迪士尼风格材质(Burley BRDF/BSDF)的头文件。
// 实现了基于物理的迪士尼着色模型，支持漫反射、金属、清漆、
// 次表面散射等多种效果。
#include "material.h"
#include "pbrt.h"

namespace pbrt {

// DisneyMaterial Declarations
// 迪士尼材质类，基于Burley等人的迪士尼BRDF/BSDF模型。
// 支持颜色、金属度、折射率、粗糙度、高光染色、各向异性、
// 光泽、清漆、镜面透射、次表面散射等多种参数。
class DisneyMaterial : public Material {
  public:
    // DisneyMaterial Public Methods
    DisneyMaterial(const std::shared_ptr<Texture<Spectrum>> &color,
                   const std::shared_ptr<Texture<Float>> &metallic,
                   const std::shared_ptr<Texture<Float>> &eta,
                   const std::shared_ptr<Texture<Float>> &roughness,
                   const std::shared_ptr<Texture<Float>> &specularTint,
                   const std::shared_ptr<Texture<Float>> &anisotropic,
                   const std::shared_ptr<Texture<Float>> &sheen,
                   const std::shared_ptr<Texture<Float>> &sheenTint,
                   const std::shared_ptr<Texture<Float>> &clearcoat,
                   const std::shared_ptr<Texture<Float>> &clearcoatGloss,
                   const std::shared_ptr<Texture<Float>> &specTrans,
                   const std::shared_ptr<Texture<Spectrum>> &scatterDistance,
                   bool thin,
                   const std::shared_ptr<Texture<Float>> &flatness,
                   const std::shared_ptr<Texture<Float>> &diffTrans,
                   const std::shared_ptr<Texture<Float>> &bumpMap)
        : color(color),
          metallic(metallic),
          eta(eta),
          roughness(roughness),
          specularTint(specularTint),
          anisotropic(anisotropic),
          sheen(sheen),
          sheenTint(sheenTint),
          clearcoat(clearcoat),
          clearcoatGloss(clearcoatGloss),
          specTrans(specTrans),
          scatterDistance(scatterDistance),
          thin(thin),
          flatness(flatness),
          diffTrans(diffTrans),
          bumpMap(bumpMap) {}
    // 计算散射函数: 根据材质参数创建BSDF/BSSRDF分布
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // DisneyMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> color;             // 基础颜色
    std::shared_ptr<Texture<Float>> metallic, eta;        // 金属度, 折射率
    std::shared_ptr<Texture<Float>> roughness, specularTint, anisotropic, sheen;  // 粗糙度, 高光染色, 各向异性, 光泽
    std::shared_ptr<Texture<Float>> sheenTint, clearcoat, clearcoatGloss;  // 光泽染色, 清漆, 清漆光泽度
    std::shared_ptr<Texture<Float>> specTrans;            // 镜面透射
    std::shared_ptr<Texture<Spectrum>> scatterDistance;   // 散射距离
    bool thin;                                            // 是否为薄表面
    std::shared_ptr<Texture<Float>> flatness, diffTrans, bumpMap;  // 平坦度, 漫透射, 凹凸贴图
};

// 创建Disney材质对象的工厂函数
DisneyMaterial *CreateDisneyMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_DISNEY_H
