
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

#ifndef PBRT_SHAPES_NURBS_H
#define PBRT_SHAPES_NURBS_H

// shapes/nurbs.h*
/**
 * @file nurbs.h
 * @brief NURBS曲面(Non-Uniform Rational B-Spline)模块
 *
 * NURBS是非均匀有理B样条曲面的实现，通过控制点网格、节点向量和阶数
 * 来定义光滑曲面。最终NURBS曲面被细分为三角形网格进行渲染。
 */
#include "pbrt.h"
#include "shape.h"
#include "geometry.h"

namespace pbrt {

/**
 * @brief 创建NURBS曲面形状
 * @param o2w 对象到世界变换
 * @param w2o 世界到对象变换
 * @param reverseOrientation 是否反转法线方向
 * @param params 参数集
 * @return 三角形网格表示的NURBS曲面
 */
std::vector<std::shared_ptr<Shape>> CreateNURBS(const Transform *o2w,
                                                const Transform *w2o,
                                                bool reverseOrientation,
                                                const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_NURBS_H
