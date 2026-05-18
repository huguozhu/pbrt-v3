
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

#ifndef PBRT_TEXTURES_MIX_H
#define PBRT_TEXTURES_MIX_H

// textures/mix.h*
// 模块功能：混合纹理（MixTexture），根据混合因子在两个纹理之间进行线性插值
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"

namespace pbrt {

// MixTexture Declarations
// 混合纹理类声明：通过amount参数在两个纹理tex1和tex2之间线性混合
template <typename T>
class MixTexture : public Texture<T> {
  public:
    // MixTexture Public Methods
    MixTexture(const std::shared_ptr<Texture<T>> &tex1,
               const std::shared_ptr<Texture<T>> &tex2,
               const std::shared_ptr<Texture<Float>> &amount)
        : tex1(tex1), tex2(tex2), amount(amount) {}
    T Evaluate(const SurfaceInteraction &si) const {
        T t1 = tex1->Evaluate(si), t2 = tex2->Evaluate(si);
        Float amt = amount->Evaluate(si);
        return (1 - amt) * t1 + amt * t2;
    }

  private:
    std::shared_ptr<Texture<T>> tex1, tex2;
    std::shared_ptr<Texture<Float>> amount;
};

MixTexture<Float> *CreateMixFloatTexture(const Transform &tex2world,
                                         const TextureParams &tp);
MixTexture<Spectrum> *CreateMixSpectrumTexture(const Transform &tex2world,
                                               const TextureParams &tp);

}  // namespace pbrt

#endif  // PBRT_TEXTURES_MIX_H
