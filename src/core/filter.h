
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

#ifndef PBRT_CORE_FILTER_H
#define PBRT_CORE_FILTER_H

// core/filter.h*
// Filter: 图像重建滤波器的抽象基类，定义像素滤波接口，
// 所有具体的滤波器（盒状、高斯、Mitchell等）都继承此类
#include "pbrt.h"
#include "geometry.h"

namespace pbrt {

// Filter Declarations
class Filter {
  public:
    // Filter Interface
    virtual ~Filter();
    // 构造函数：初始化滤波器半径和半径倒数
    Filter(const Vector2f &radius)
        : radius(radius), invRadius(Vector2f(1 / radius.x, 1 / radius.y)) {}
    // Evaluate: 计算滤波器在给定位置p处的权重值（纯虚函数）
    virtual Float Evaluate(const Point2f &p) const = 0;

    // Filter Public Data
    const Vector2f radius, invRadius;  // 滤波器半径及倒数
};

}  // namespace pbrt

#endif  // PBRT_CORE_FILTER_H
