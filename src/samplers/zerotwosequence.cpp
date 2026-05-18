
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


// samplers/zerotwosequence.cpp*
/**
 * @file zerotwosequence.cpp
 * @brief (0,2)-序列采样器(ZeroTwoSequenceSampler)的实现
 *
 * 使用(0,2)-低差异序列生成像素样本。
 * 1D维度使用Van Der Corput序列（基数为2的Radical Inverse），
 * 2D维度使用Sobol' 2D序列，确保良好的多维均匀分布。
 */
#include "samplers/zerotwosequence.h"
#include "lowdiscrepancy.h"
#include "paramset.h"
#include "stats.h"

namespace pbrt {

// ZeroTwoSequenceSampler Method Definitions / (0,2)-序列采样器方法定义
/**
 * @brief 构造函数
 *
 * 将每像素样本数向上取整为2的幂，因为(0,2)-序列要求
 * 样本数为2的幂次。
 * @param samplesPerPixel 每像素样本数（将被向上取整为2的幂）
 * @param nSampledDimensions 需要独立采样的维度数
 */
ZeroTwoSequenceSampler::ZeroTwoSequenceSampler(int64_t samplesPerPixel,
                                               int nSampledDimensions)
    : PixelSampler(RoundUpPow2(samplesPerPixel), nSampledDimensions) {
    if (!IsPowerOf2(samplesPerPixel))
        Warning(
            "Pixel samples being rounded up to power of 2 "
            "(from %" PRId64 " to %" PRId64 ").",
            samplesPerPixel, RoundUpPow2(samplesPerPixel));
}

/**
 * @brief 开始新像素的采样
 *
 * 使用(0,2)-序列生成当前像素的所有样本：
 * - 1D样本使用Van Der Corput序列（基数为2的Radical Inverse）
 * - 2D样本使用Sobol' 2D序列
 * 数组样本同样使用对应的低差异序列生成。
 * @param p 当前像素坐标
 */
void ZeroTwoSequenceSampler::StartPixel(const Point2i &p) {
    ProfilePhase _(Prof::StartPixel);
    // Generate 1D and 2D pixel sample components using $(0,2)$-sequence
    // 使用(0,2)-序列生成1D和2D像素样本分量
    for (size_t i = 0; i < samples1D.size(); ++i)
        VanDerCorput(1, samplesPerPixel, &samples1D[i][0], rng);
    for (size_t i = 0; i < samples2D.size(); ++i)
        Sobol2D(1, samplesPerPixel, &samples2D[i][0], rng);

    // Generate 1D and 2D array samples using $(0,2)$-sequence
    // 使用(0,2)-序列生成1D和2D数组样本
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i)
        VanDerCorput(samples1DArraySizes[i], samplesPerPixel,
                     &sampleArray1D[i][0], rng);
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i)
        Sobol2D(samples2DArraySizes[i], samplesPerPixel, &sampleArray2D[i][0],
                rng);
    PixelSampler::StartPixel(p);
}

/**
 * @brief 克隆采样器
 *
 * 创建当前采样器的副本，为新副本设置独立的随机数序列
 * （用于随机打乱(Scramble)操作）。
 * @param seed 新采样器的随机数种子
 * @return 新采样器的唯一指针
 */
std::unique_ptr<Sampler> ZeroTwoSequenceSampler::Clone(int seed) {
    ZeroTwoSequenceSampler *lds = new ZeroTwoSequenceSampler(*this);
    lds->rng.SetSequence(seed);
    return std::unique_ptr<Sampler>(lds);
}

/**
 * @brief 创建ZeroTwoSequenceSampler的工厂函数
 * @param params 参数集，从中读取"pixelsamples"和"dimensions"参数
 * @return 创建的ZeroTwoSequenceSampler实例指针
 */
ZeroTwoSequenceSampler *CreateZeroTwoSequenceSampler(const ParamSet &params) {
    int nsamp = params.FindOneInt("pixelsamples", 16);
    int sd = params.FindOneInt("dimensions", 4);
    if (PbrtOptions.quickRender) nsamp = 1;
    return new ZeroTwoSequenceSampler(nsamp, sd);
}

}  // namespace pbrt
