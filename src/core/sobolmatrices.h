
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

#ifndef PBRT_CORE_SOBOLMATRICES_H
#define PBRT_CORE_SOBOLMATRICES_H

// core/sobolmatrices.h*
// Sobol序列矩阵数据模块：包含用于生成Sobol准随机序列的
// 预计算生成矩阵。Sobol序列是一种低差异序列(LDS)，
// 在渲染中用于重要性采样和蒙特卡罗积分，能够减少方差。
#include "pbrt.h"

namespace pbrt {

// Sobol Matrix Declarations
static PBRT_CONSTEXPR int NumSobolDimensions = 1024;
static PBRT_CONSTEXPR int SobolMatrixSize = 52;
extern const uint32_t SobolMatrices32[NumSobolDimensions * SobolMatrixSize];
extern const uint64_t SobolMatrices64[NumSobolDimensions * SobolMatrixSize];
extern const uint64_t VdCSobolMatrices[][SobolMatrixSize];
extern const uint64_t VdCSobolMatricesInv[][SobolMatrixSize];

}  // namespace pbrt

#endif  // PBRT_CORE_SOBOLMATRICES_H
