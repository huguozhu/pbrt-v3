
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

#ifndef PBRT_TEXTURES_UV_H
#define PBRT_TEXTURES_UV_H

// textures/uv.h*
// 模块功能：UV纹理（UVTexture），将纹理坐标(u,v)直接映射为颜色值，用于可视化UV映射
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

namespace pbrt {

// UVTexture Declarations
// UV纹理类声明：将纹理坐标直接显示为RGB颜色（R=u, G=v, B=0）
class UVTexture : public Texture<Spectrum> {
  public:
    // UVTexture Public Methods
    UVTexture(std::unique_ptr<TextureMapping2D> mapping)
        : mapping(std::move(mapping)) {}
    Spectrum Evaluate(const SurfaceInteraction &si) const {
        Vector2f dstdx, dstdy;
        Point2f st = mapping->Map(si, &dstdx, &dstdy);
        Float rgb[3] = {st[0] - std::floor(st[0]), st[1] - std::floor(st[1]),
                        0};
        return Spectrum::FromRGB(rgb);
    }

  private:
    std::unique_ptr<TextureMapping2D> mapping;
};

Texture<Float> *CreateUVFloatTexture(const Transform &tex2world,
                                     const TextureParams &tp);
UVTexture *CreateUVSpectrumTexture(const Transform &tex2world,
                                   const TextureParams &tp);

}  // namespace pbrt

#endif  // PBRT_TEXTURES_UV_H
