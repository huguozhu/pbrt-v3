
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

#ifndef PBRT_SHAPES_HEIGHTFIELD_H
#define PBRT_SHAPES_HEIGHTFIELD_H

// shapes/heightfield.h*
/**
 * @file heightfield.h
 * @brief 高度场(Heightfield)几何体模块
 *
 * 通过高度图(heigh map)数据生成三角形网格地形。
 * 输入是一个nx*ny的网格高度值数组，自动生成三角形网格。
 */
#include "shape.h"

namespace pbrt {

// Heightfield Declarations / 高度场声明
/**
 * @brief 创建高度场形状
 * @param o2w 对象到世界变换
 * @param w2o 世界到对象变换
 * @param ro 是否反转法线方向
 * @param params 参数集(包含nu, nv, Pz等)
 * @return 生成的三角形网格
 */
std::vector<std::shared_ptr<Shape>> CreateHeightfield(const Transform *o2w,
                                                      const Transform *w2o,
                                                      bool ro,
                                                      const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_HEIGHTFIELD_H
