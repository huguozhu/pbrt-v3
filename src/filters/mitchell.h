
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

#ifndef PBRT_FILTERS_MITCHELL_H
#define PBRT_FILTERS_MITCHELL_H

// filters/mitchell.h*
// MitchellFilter: Mitchell-Netravali滤波器，使用分段三次多项式逼近理想重建滤波器，
// 通过B和C参数控制滤波器的行为，B=1/3 C=1/3时为经典Mitchell滤波
#include "filter.h"

namespace pbrt {

// Mitchell Filter Declarations
class MitchellFilter : public Filter {
  public:
    // MitchellFilter Public Methods
    // 构造函数，初始化滤波器半径和Mitchell参数B、C
    MitchellFilter(const Vector2f &radius, Float B, Float C)
        : Filter(radius), B(B), C(C) {}
    // Evaluate: 计算Mitchell滤波器在给定位置p处的权重值
    Float Evaluate(const Point2f &p) const;
    // Mitchell1D: 一维Mitchell滤波函数，使用分段三次多项式
    Float Mitchell1D(Float x) const {
        x = std::abs(2 * x);
        if (x > 1)
            // 在[1,2)区间上的三次多项式
            return ((-B - 6 * C) * x * x * x + (6 * B + 30 * C) * x * x +
                    (-12 * B - 48 * C) * x + (8 * B + 24 * C)) *
                   (1.f / 6.f);
        else
            // 在[0,1]区间上的三次多项式
            return ((12 - 9 * B - 6 * C) * x * x * x +
                    (-18 + 12 * B + 6 * C) * x * x + (6 - 2 * B)) *
                   (1.f / 6.f);
    }

  private:
    const Float B, C;  // Mitchell-Netravali滤波器的两个形状参数
};

MitchellFilter *CreateMitchellFilter(const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_FILTERS_MITCHELL_H
