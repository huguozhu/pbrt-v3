
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

#ifndef PBRT_MATERIALS_MIXMAT_H
#define PBRT_MATERIALS_MIXMAT_H

// materials/mixmat.h*
// 文件描述: 混合材质的头文件。实现两种材质的加权混合，
// 通过纹理权重控制混合比例。
#include "pbrt.h"
#include "material.h"

namespace pbrt {

// MixMaterial Declarations
// 混合材质类，将两种材质(m1, m2)按纹理权重(scale)进行混合，
// 用于创建复杂的组合材质效果。
class MixMaterial : public Material {
  public:
    // MixMaterial Public Methods
    MixMaterial(const std::shared_ptr<Material> &m1,
                const std::shared_ptr<Material> &m2,
                const std::shared_ptr<Texture<Spectrum>> &scale)
        : m1(m1), m2(m2), scale(scale) {}
    // 计算散射函数: 混合两种材质的BSDF
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // MixMaterial Private Data
    std::shared_ptr<Material> m1, m2;         // 被混合的两种材质
    std::shared_ptr<Texture<Spectrum>> scale; // 混合权重纹理
};

// 创建混合材质对象的工厂函数
MixMaterial *CreateMixMaterial(const TextureParams &mp,
                               const std::shared_ptr<Material> &m1,
                               const std::shared_ptr<Material> &m2);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_MIXMAT_H
