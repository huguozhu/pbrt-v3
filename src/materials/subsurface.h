
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

#ifndef PBRT_MATERIALS_SUBSURFACE_H
#define PBRT_MATERIALS_SUBSURFACE_H

// materials/subsurface.h*
// 文件描述: 次表面散射材质的头文件。实现了基于BSSRDF的次表面散射模型，
// 支持漫反射和微表面高光反射。
#include "pbrt.h"
#include "material.h"
#include "reflection.h"
#include "bssrdf.h"

namespace pbrt {

// SubsurfaceMaterial Declarations
// 次表面散射材质类，支持通过sigma_a(吸收系数)和sigma_s(散射系数)控制。
// 包含TabulatedBSSRDF用于次表面散射模拟，以及可选的镜面/微表面反射/透射。
class SubsurfaceMaterial : public Material {
  public:
    // SubsurfaceMaterial Public Methods
    SubsurfaceMaterial(Float scale,
                       const std::shared_ptr<Texture<Spectrum>> &Kr,
                       const std::shared_ptr<Texture<Spectrum>> &Kt,
                       const std::shared_ptr<Texture<Spectrum>> &sigma_a,
                       const std::shared_ptr<Texture<Spectrum>> &sigma_s,
                       Float g, Float eta,
                       const std::shared_ptr<Texture<Float>> &uRoughness,
                       const std::shared_ptr<Texture<Float>> &vRoughness,
                       const std::shared_ptr<Texture<Float>> &bumpMap,
                       bool remapRoughness)
        : scale(scale),
          Kr(Kr),
          Kt(Kt),
          sigma_a(sigma_a),
          sigma_s(sigma_s),
          uRoughness(uRoughness),
          vRoughness(vRoughness),
          bumpMap(bumpMap),
          eta(eta),
          remapRoughness(remapRoughness),
          table(100, 64) {
        ComputeBeamDiffusionBSSRDF(g, eta, &table);
    }
    // 计算散射函数: 创建BSDF和次表面散射BSSRDF
    void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena,
                                    TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // SubsurfaceMaterial Private Data
    const Float scale;                                     // 缩放因子
    std::shared_ptr<Texture<Spectrum>> Kr, Kt, sigma_a, sigma_s;  // 反射/透射/吸收/散射系数
    std::shared_ptr<Texture<Float>> uRoughness, vRoughness;       // 表面粗糙度
    std::shared_ptr<Texture<Float>> bumpMap;                      // 凹凸贴图
    const Float eta;                                        // 折射率
    const bool remapRoughness;                              // 是否重映射粗糙度
    BSSRDFTable table;                                      // BSSRDF预计算表
};

// 创建次表面散射材质对象的工厂函数
SubsurfaceMaterial *CreateSubsurfaceMaterial(const TextureParams &mp);

}  // namespace pbrt

#endif  // PBRT_MATERIALS_SUBSURFACE_H
