
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

#ifndef PBRT_MATERIALS_SUBSTRATE_H
#define PBRT_MATERIALS_SUBSTRATE_H

// materials/substrate.h*
// 文件描述: 基底层材质(FresnelBlend)的头文件。实现了漫反射基底加微表面高光的模型，
// 使用FresnelBlend BSDF混合漫反射和高光反射。
#include "pbrt.h"
#include "material.h"

namespace pbrt {

// SubstrateMaterial Declarations
// 基底层材质类，使用FresnelBlend BSDF融合漫反射(Kd)和高光(Ks)分量，
// 并使用各向异性粗糙度(nu, nv)控制微表面分布。
class SubstrateMaterial : public Material {
  public:
    // SubstrateMaterial Public Methods
    SubstrateMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
                      const std::shared_ptr<Texture<Spectrum>> &Ks,
                      const std::shared_ptr<Texture<Float>> &nu,
                      const std::shared_ptr<Texture<Float>> &nv,
                      const std::shared_ptr<Texture<Float>> &bumpMap,
                      bool remapRoughness)
        : Kd(Kd),
          Ks(Ks),
          nu(nu),
          nv(nv),
          bumpMap(bumpMap),
          remapRoughness(remapRoughness) {}
    // 计算散射函数: 创建基底层材质的FresnelBlend BSDF
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // SubstrateMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> Kd, Ks;  // 漫反射和高光颜色
    std::shared_ptr<Texture<Float>> nu, nv;      // U和V方向的粗糙度
    std::shared_ptr<Texture<Float>> bumpMap;     // 凹凸贴图
    bool remapRoughness;                         // 是否重映射粗糙度
};

// 创建基底层材质对象的工厂函数
SubstrateMaterial *CreateSubstrateMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSTRATE_H
