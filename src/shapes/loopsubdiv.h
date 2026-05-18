
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

#ifndef PBRT_SHAPES_LOOPSUBDIV_H
#define PBRT_SHAPES_LOOPSUBDIV_H

// shapes/loopsubdiv.h*
/**
 * @file loopsubdiv.h
 * @brief Loop subdivision曲面(LoopSubdiv)模块
 *
 * 实现了Loop细分曲面算法，将粗糙的控制网格通过多次细分生成光滑曲面。
 * 最终输出为三角形网格。
 */
#include "shape.h"

namespace pbrt {

// LoopSubdiv Declarations / Loop细分曲面声明
/**
 * @brief 创建Loop细分曲面形状
 * @param o2w 对象到世界变换
 * @param w2o 世界到对象变换
 * @param reverseOrientation 是否反转法线方向
 * @param params 参数集(包含细分层级、控制网格顶点和索引)
 * @return 细分后的三角形网格
 */
std::vector<std::shared_ptr<Shape>> CreateLoopSubdiv(const Transform *o2w,
                                                     const Transform *w2o,
                                                     bool reverseOrientation,
                                                     const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_LOOPSUBDIV_H
