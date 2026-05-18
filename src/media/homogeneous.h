
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

#ifndef PBRT_MEDIA_HOMOGENEOUS_H
#define PBRT_MEDIA_HOMOGENEOUS_H

// media/homogeneous.h*
// 文件描述: 均匀参与介质的头文件。实现空间均匀分布的散射和吸收介质，
// 支持各向异性相位函数(Henyey-Greenstein)。
#include "medium.h"

namespace pbrt {

// HomogeneousMedium Declarations
// 均匀介质类，介质的吸收系数(sigma_a)和散射系数(sigma_s)在空间中为常数。
// 使用解析方法计算透射率和采样散射事件。
class HomogeneousMedium : public Medium {
  public:
    // HomogeneousMedium Public Methods
    HomogeneousMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g)
        : sigma_a(sigma_a),
          sigma_s(sigma_s),
          sigma_t(sigma_s + sigma_a),
          g(g) {}
    // 计算光线穿过均匀介质的透射率(解析解)
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;
    // 采样介质中的散射事件(解析采样)
    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;

  private:
    // HomogeneousMedium Private Data
    const Spectrum sigma_a, sigma_s, sigma_t;  // 吸收系数、散射系数、总衰减系数
    const Float g;                              // 相位函数各向异性参数
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_HOMOGENEOUS_H
