
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

#ifndef PBRT_CORE_INTERPOLATION_H
#define PBRT_CORE_INTERPOLATION_H

// core/interpolation.h*
// 插值工具: 提供Catmull-Rom样条插值和Fourier级数插值函数，
// 用于BRDF/BSSRDF等数据的平滑插值和采样
#include "pbrt.h"

namespace pbrt {

// Spline Interpolation Declarations
// CatmullRom: 在给定节点和值上进行Catmull-Rom样条插值
Float CatmullRom(int size, const Float *nodes, const Float *values, Float x);
// CatmullRomWeights: 计算插值权重和偏移，用于快速插值
bool CatmullRomWeights(int size, const Float *nodes, Float x, int *offset,
                       Float *weights);
// SampleCatmullRom: 从一维Catmull-Rom样条表示的概率分布中采样
Float SampleCatmullRom(int size, const Float *nodes, const Float *f,
                       const Float *cdf, Float sample, Float *fval = nullptr,
                       Float *pdf = nullptr);
// SampleCatmullRom2D: 从二维Catmull-Rom样条表示的概率分布中采样
Float SampleCatmullRom2D(int size1, int size2, const Float *nodes1,
                         const Float *nodes2, const Float *values,
                         const Float *cdf, Float alpha, Float sample,
                         Float *fval = nullptr, Float *pdf = nullptr);
// IntegrateCatmullRom: 对Catmull-Rom样条进行数值积分，生成CDF
Float IntegrateCatmullRom(int n, const Float *nodes, const Float *values,
                          Float *cdf);
// InvertCatmullRom: 对Catmull-Rom样条的CDF进行逆变换采样
Float InvertCatmullRom(int n, const Float *x, const Float *values, Float u);

// Fourier Interpolation Declarations
// Fourier: 计算Fourier级数在给定角度的值，用于各向异性BRDF
Float Fourier(const Float *a, int m, double cosPhi);
// SampleFourier: 根据Fourier级数表示的分布进行采样
Float SampleFourier(const Float *ak, const Float *recip, int m, Float u,
                    Float *pdf, Float *phiPtr);

}  // namespace pbrt

#endif  // PBRT_CORE_INTERPOLATION_H
