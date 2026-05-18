
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

#ifndef PBRT_MATERIALS_MIRROR_H
#define PBRT_MATERIALS_MIRROR_H

// materials/mirror.h*
// 文件描述: 镜面材质的头文件。实现完美镜面反射(Snell反射)。
#include "pbrt.h"
#include "material.h"

namespace pbrt {

// MirrorMaterial Declarations
// 镜面材质类，使用完美镜面反射模型(SpecularReflection)，
// 配合FresnelNoOp(无菲涅耳效应)实现理想的镜面反射。
class MirrorMaterial : public Material {
  public:
    // MirrorMaterial Public Methods
    MirrorMaterial(const std::shared_ptr<Texture<Spectrum>> &r,
                   const std::shared_ptr<Texture<Float>> &bump) {
        Kr = r;
        bumpMap = bump;
    }
    // 计算散射函数: 创建镜面反射BSDF
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // MirrorMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> Kr;    // 镜面反射率
    std::shared_ptr<Texture<Float>> bumpMap;  // 凹凸贴图
};

// 创建镜面材质对象的工厂函数
MirrorMaterial *CreateMirrorMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_MIRROR_H
