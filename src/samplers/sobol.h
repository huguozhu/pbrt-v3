
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

#ifndef PBRT_SAMPLERS_SOBOL_H
#define PBRT_SAMPLERS_SOBOL_H

// samplers/sobol.h*
/**
 * @file sobol.h
 * @brief Sobol'序列采样器(SobolSampler)模块
 *
 * 使用Sobol'低差异序列生成多维样本点。
 * Sobol'序列是一种(0,2)-序列，具有良好的低差异性质，
 * 适合用于蒙特卡洛渲染中的重要性采样。
 */
#include "sampler.h"

namespace pbrt {

// SobolSampler Declarations / Sobol采样器声明
/**
 * @brief Sobol'低差异序列采样器
 *
 * 使用Sobol'序列生成低差异样本点。Sobol'序列是基数为2的
 * (0,2)-序列，样本在联合维度上具有良好的分布均匀性。
 * 通过SobolIntervalToIndex函数将像素坐标映射到Sobol'序列索引。
 */
class SobolSampler : public GlobalSampler {
  public:
    // SobolSampler Public Methods / 公有方法
    /**
     * @brief 克隆采样器
     */
    std::unique_ptr<Sampler> Clone(int seed);
    /**
     * @brief 构造函数
     * @param samplesPerPixel 每像素样本数（将被向上取整为2的幂）
     * @param sampleBounds 采样区域边界
     */
    SobolSampler(int64_t samplesPerPixel, const Bounds2i &sampleBounds)
        : GlobalSampler(RoundUpPow2(samplesPerPixel)),
          sampleBounds(sampleBounds) {
        if (!IsPowerOf2(samplesPerPixel))
            Warning("Non power-of-two sample count rounded up to %" PRId64
                    " for SobolSampler.",
                    this->samplesPerPixel);
        resolution = RoundUpPow2(
            std::max(sampleBounds.Diagonal().x, sampleBounds.Diagonal().y));
        log2Resolution = Log2Int(resolution);
        if (resolution > 0) CHECK_EQ(1 << log2Resolution, resolution);
    }
    /**
     * @brief 计算样本的全局Sobol'序列索引
     */
    int64_t GetIndexForSample(int64_t sampleNum) const;
    /**
     * @brief 获取指定维度的样本值
     */
    Float SampleDimension(int64_t index, int dimension) const;

  private:
    // SobolSampler Private Data / 私有数据
    const Bounds2i sampleBounds; /**< 采样区域边界 */
    int resolution;              /**< 采样分辨率（2的幂） */
    int log2Resolution;          /**< 分辨率的以2为底的对数 */
};

/**
 * @brief 创建SobolSampler的工厂函数
 * @param params 参数集
 * @param sampleBounds 采样区域边界
 */
SobolSampler *CreateSobolSampler(const ParamSet &params,
                                 const Bounds2i &sampleBounds);

}  // namespace pbrt

#endif  // PBRT_SAMPLERS_SOBOL_H
