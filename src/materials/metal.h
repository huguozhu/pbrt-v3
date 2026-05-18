
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

#ifndef PBRT_MATERIALS_METAL_H
#define PBRT_MATERIALS_METAL_H

// materials/metal.h*
// 文件描述: 金属材质的头文件。使用微表面模型和菲涅耳导体反射来描述金属表面。
#include "pbrt.h"
#include "material.h"
#include "spectrum.h"

namespace pbrt {

// MetalMaterial Declarations
// 金属材质类，使用Trowbridge-Reitz微表面分布和导体菲涅耳反射模型。
// 支持各向异性粗糙度和折射率(eta)、消光系数(k)的纹理控制。
class MetalMaterial : public Material {
  public:
    // MetalMaterial Public Methods
    MetalMaterial(const std::shared_ptr<Texture<Spectrum>> &eta,
                  const std::shared_ptr<Texture<Spectrum>> &k,
                  const std::shared_ptr<Texture<Float>> &rough,
                  const std::shared_ptr<Texture<Float>> &urough,
                  const std::shared_ptr<Texture<Float>> &vrough,
                  const std::shared_ptr<Texture<Float>> &bump,
                  bool remapRoughness);
    // 计算散射函数: 创建金属材质的微表面反射BSDF
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // MetalMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> eta, k;
    std::shared_ptr<Texture<Float>> roughness, uRoughness, vRoughness;
    std::shared_ptr<Texture<Float>> bumpMap;
    bool remapRoughness;
};

// 创建金属材质对象的工厂函数
MetalMaterial *CreateMetalMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_METAL_H
