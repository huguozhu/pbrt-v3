
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

#ifndef PBRT_FILTERS_SINC_H
#define PBRT_FILTERS_SINC_H

// filters/sinc.h*
// LanczosSincFilter: Lanczos窗口化的Sinc滤波器，使用sinc函数作为理想低通滤波器，
// 通过Lanczos窗口（tau参数）截断滤波器的空间范围以平衡振铃和模糊
#include "filter.h"

namespace pbrt {

// Sinc Filter Declarations
class LanczosSincFilter : public Filter {
  public:
    // LanczosSincFilter Public Methods
    // 构造函数，初始化滤波器半径和Lanczos窗口参数tau
    LanczosSincFilter(const Vector2f &radius, Float tau)
        : Filter(radius), tau(tau) {}
    // Evaluate: 计算LanczosSinc滤波器在给定位置p处的权重值
    Float Evaluate(const Point2f &p) const;
    // Sinc: 标准的sinc函数 sin(pi*x)/(pi*x)
    Float Sinc(Float x) const {
        x = std::abs(x);
        if (x < 1e-5) return 1;  // 避免除零
        return std::sin(Pi * x) / (Pi * x);
    }
    // WindowedSinc: 使用Lanczos窗口截断的sinc函数
    Float WindowedSinc(Float x, Float radius) const {
        x = std::abs(x);
        if (x > radius) return 0;
        Float lanczos = Sinc(x / tau);
        return Sinc(x) * lanczos;
    }

  private:
    const Float tau;  // Lanczos窗口参数，控制窗口大小
};

LanczosSincFilter *CreateSincFilter(const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_FILTERS_SINC_H
