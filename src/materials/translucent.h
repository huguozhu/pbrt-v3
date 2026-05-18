
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

#ifndef PBRT_MATERIALS_TRANSLUCENT_H
#define PBRT_MATERIALS_TRANSLUCENT_H

// materials/translucent.h*
// 文件描述: 半透明材质的头文件。实现了同时具有漫反射/透射和微表面高光反射/透射的材质。
#include "pbrt.h"
#include "material.h"

namespace pbrt {

// TranslucentMaterial Declarations
// 半透明材质类，支持漫反射(Kd)、高光(Ks)、反射比例(reflect)、透射比例(transmit)，
// 以及粗糙度(roughness)控制，用于模拟半透明表面。
class TranslucentMaterial : public Material {
  public:
    // TranslucentMaterial Public Methods
    TranslucentMaterial(const std::shared_ptr<Texture<Spectrum>> &kd,
                        const std::shared_ptr<Texture<Spectrum>> &ks,
                        const std::shared_ptr<Texture<Float>> &rough,
                        const std::shared_ptr<Texture<Spectrum>> &refl,
                        const std::shared_ptr<Texture<Spectrum>> &trans,
                        const std::shared_ptr<Texture<Float>> &bump,
                        bool remap) {
        Kd = kd;
        Ks = ks;
        roughness = rough;
        reflect = refl;
        transmit = trans;
        bumpMap = bump;
        remapRoughness = remap;
    }
    // 计算散射函数: 创建半透明材质的BSDF
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // TranslucentMaterial Private Data
    std::shared_ptr<Texture<Spectrum>> Kd, Ks;   // 漫反射和高光颜色
    std::shared_ptr<Texture<Float>> roughness;    // 表面粗糙度
    std::shared_ptr<Texture<Spectrum>> reflect, transmit;  // 反射和透射比例
    std::shared_ptr<Texture<Float>> bumpMap;      // 凹凸贴图
    bool remapRoughness;                          // 是否重映射粗糙度
};

// 创建半透明材质对象的工厂函数
TranslucentMaterial *CreateTranslucentMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_TRANSLUCENT_H
