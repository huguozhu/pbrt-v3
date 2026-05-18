
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

#ifndef PBRT_FILTERS_GAUSSIAN_H
#define PBRT_FILTERS_GAUSSIAN_H

// filters/gaussian.h*
// GaussianFilter: 高斯滤波器，使用高斯函数作为图像重建滤波器，
// 通过alpha参数控制高斯函数的衰减速度，产生平滑的重建结果
#include "filter.h"

namespace pbrt {

// Gaussian Filter Declarations
class GaussianFilter : public Filter {
  public:
    // GaussianFilter Public Methods
    // 构造函数，初始化滤波器半径、alpha值及边界处的指数值
    GaussianFilter(const Vector2f &radius, Float alpha)
        : Filter(radius),
          alpha(alpha),
          expX(std::exp(-alpha * radius.x * radius.x)),
          expY(std::exp(-alpha * radius.y * radius.y)) {}
    // Evaluate: 计算高斯滤波器在给定位置p处的权重值
    Float Evaluate(const Point2f &p) const;

  private:
    // GaussianFilter Private Data
    const Float alpha;       // 高斯函数的衰减速率参数
    const Float expX, expY;  // x和y方向滤波器边界处的高斯函数值，用于归一化

    // GaussianFilter Utility Functions
    // Gaussian: 一维高斯函数，返回exp(-alpha*d^2)减去边界偏移量的结果
    Float Gaussian(Float d, Float expv) const {
        return std::max((Float)0, Float(std::exp(-alpha * d * d) - expv));
    }
};

GaussianFilter *CreateGaussianFilter(const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_FILTERS_GAUSSIAN_H
