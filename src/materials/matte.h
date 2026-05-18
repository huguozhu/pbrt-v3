
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

#ifndef PBRT_MATERIALS_MATTE_H
#define PBRT_MATERIALS_MATTE_H

// materials/matte.h*
#include "pbrt.h"
#include "material.h"

namespace pbrt {

// MatteMaterial Declarations
// 漫反射/哑光材质：漫反射表面的渲染实现
// 在内部，根据sigma参数决定BSDF类型：若sigma==0则创建Lambertian反射；
// 否则创建OrenNayar反射模型(sigma控制OrenNayar的粗糙度参数)
class MatteMaterial : public Material {
  public:
    // MatteMaterial Public Methods
    MatteMaterial(const std::shared_ptr<Texture<Spectrum>> &Kd,
                  const std::shared_ptr<Texture<Float>> &sigma,
                  const std::shared_ptr<Texture<Float>> &bumpMap)
        : Kd(Kd), sigma(sigma), bumpMap(bumpMap) {}
    // 计算散射函数: 创建漫反射BSDF(Lambertian或OrenNayar)
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // MatteMaterial Private Data
	// Kd-->漫反射颜色
    std::shared_ptr<Texture<Spectrum>> Kd;
	// sigma控制表面粗糙度: sigma==0 --> Lambertian反射，否则OrenNayar反射
	// bumpMap-->用于产生表面法线方向扰动
    std::shared_ptr<Texture<Float>> sigma, bumpMap;
};

// 创建哑光材质对象的工厂函数
MatteMaterial *CreateMatteMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_MATTE_H
