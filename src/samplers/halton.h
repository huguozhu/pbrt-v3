
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

#ifndef PBRT_SAMPLERS_HALTON_H
#define PBRT_SAMPLERS_HALTON_H

// samplers/halton.h*
#include "sampler.h"
#include "lowdiscrepancy.h"

namespace pbrt {

/**
 * @brief Halton低差异序列采样器
 *
 * 使用Radical Inversion(基数反转)生成低差异序列。
 * 对第0维使用基数2，第1维使用基数3，第2维及以上使用质数表中的质数。
 * 支持加扰(Scrambling)和像素中心采样。
 */
// HaltonSampler Declarations
// �����ٲ�����ʹ�õ����㷨(Radical Inversion)�õ���ƫ������
// �͹�Ĭ˹��(Hammersley)��������
class HaltonSampler : public GlobalSampler {
  public:
    // HaltonSampler Public Methods / Halton采样器公有方法
    /**
     * @brief 构造函数
     * @param nsamp 每像素样本数
     * @param sampleBounds 采样区域边界
     * @param sampleAtCenter 是否强制在像素中心采样
     */
    HaltonSampler(int nsamp, const Bounds2i &sampleBounds,
                  bool sampleAtCenter = false);
    /**
     * @brief 获取样本的全局索引
     */
    int64_t GetIndexForSample(int64_t sampleNum) const;
    /**
     * @brief 获取指定维度的样本值
     */
    Float SampleDimension(int64_t index, int dimension) const;
    /**
     * @brief 克隆采样器
     */
    std::unique_ptr<Sampler> Clone(int seed);

  private:
    // HaltonSampler Private Data / Halton采样器私有数据
    static std::vector<uint16_t> radicalInversePermutations; /**< 基数反转置换表 */
    Point2i baseScales, baseExponents;  /**< 基数2和3的缩放指数 */
    int sampleStride;                   /**< 像素访问的样本步长 */
    int multInverse[2];                 /**< 乘法逆元 */
    mutable Point2i pixelForOffset = Point2i(std::numeric_limits<int>::max(),
                                             std::numeric_limits<int>::max());
    mutable int64_t offsetForCurrentPixel;
    // Added after book publication: force all image samples to be at the
    // center of the pixel area. / 书后添加：强制所有图像样本在像素中心
    bool sampleAtPixelCenter;

    // HaltonSampler Private Methods / 私有方法
    /**
     * @brief 获取指定维度的置换表
     */
    const uint16_t *PermutationForDimension(int dim) const {
        if (dim >= PrimeTableSize)
            LOG(FATAL) << StringPrintf("HaltonSampler can only sample %d "
                                       "dimensions.", PrimeTableSize);
        return &radicalInversePermutations[PrimeSums[dim]];
    }
};

/**
 * @brief 创建Halton采样器的工厂函数
 */
HaltonSampler *CreateHaltonSampler(const ParamSet &params,
                                   const Bounds2i &sampleBounds);

}  // namespace pbrt

#endif  // PBRT_SAMPLERS_HALTON_H
