
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

#ifndef PBRT_SAMPLERS_ZEROTWOSEQUENCE_H
#define PBRT_SAMPLERS_ZEROTWOSEQUENCE_H

// samplers/zerotwosequence.h*
/**
 * @file zerotwosequence.h
 * @brief (0,2)-序列采样器(ZeroTwoSequenceSampler)模块
 *
 * 使用Van Der Corput序列（1D）和Sobol' 2D序列生成
 * 具有(0,2)-低差异性质的样本点。适用于需要良好
 * 多维均匀分布的蒙特卡洛渲染。
 */
#include "sampler.h"

namespace pbrt {

// ZeroTwoSequenceSampler Declarations / (0,2)-序列采样器声明
/**
 * @brief (0,2)-序列采样器
 *
 * 使用(0,2)-序列生成低差异样本。1D维度使用Van Der Corput序列
 * （基数为2的Radical Inverse），2D维度使用Sobol' 2D序列。
 * 每像素样本数必须为2的幂。
 */
class ZeroTwoSequenceSampler : public PixelSampler {
  public:
    // ZeroTwoSequenceSampler Public Methods / 公有方法
    /**
     * @brief 构造函数
     * @param samplesPerPixel 每像素样本数（将被向上取整为2的幂）
     * @param nSampledDimensions 需要独立采样的维度数
     */
    ZeroTwoSequenceSampler(int64_t samplesPerPixel, int nSampledDimensions = 4);
    /**
     * @brief 开始新像素的采样，使用(0,2)-序列生成样本
     */
    void StartPixel(const Point2i &);
    /**
     * @brief 克隆采样器
     */
    std::unique_ptr<Sampler> Clone(int seed);
    /**
     * @brief 将样本数向上取整到2的幂（(0,2)-序列要求）
     */
    int RoundCount(int count) const { return RoundUpPow2(count); }
};

/**
 * @brief 创建ZeroTwoSequenceSampler的工厂函数
 * @param params 参数集
 * @return 创建的ZeroTwoSequenceSampler实例指针
 */
ZeroTwoSequenceSampler *CreateZeroTwoSequenceSampler(const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SAMPLERS_ZEROTWOSEQUENCE_H
