
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

#ifndef PBRT_CORE_MATERIAL_H
#define PBRT_CORE_MATERIAL_H

// core/material.h*
// 材质: 定义材质接口，负责计算表面散射函数(BSDF/BSSRDF)，
// 所有具体材质类型继承此类并实现散射函数计算
#include "pbrt.h"
#include "memory.h"

namespace pbrt {

// TransportMode Declarations
// TransportMode: 表示光传输模式，Radiance（从光源到相机）或 Importance（从相机到光源），用于路径追踪中的重要性驱动
enum class TransportMode { Radiance, Importance };

// Material Declarations
class Material {
  public:
    // Material Interface
    virtual void ComputeScatteringFunctions(SurfaceInteraction *si,
                                            MemoryArena &arena,
                                            TransportMode mode,
                                            bool allowMultipleLobes) const = 0;
    virtual ~Material();
    // Bump: 凹凸映射，通过纹理扰动表面法线模拟凹凸效果
    // 公式: p'(u,v) = p(u,v) + d(u,v)n(u,v)
    static void Bump(const std::shared_ptr<Texture<Float>> &d,
                     SurfaceInteraction *si);
};

}  // namespace pbrt

#endif  // PBRT_CORE_MATERIAL_H
