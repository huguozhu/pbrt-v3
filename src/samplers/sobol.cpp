
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


// samplers/sobol.cpp*
/**
 * @file sobol.cpp
 * @brief Sobol'序列采样器(SobolSampler)的实现
 *
 * 使用Sobol'低差异序列生成多维样本点。
 * 通过SobolIntervalToIndex函数将像素坐标和样本序号映射到
 * Sobol'序列的全局索引。
 */
#include "samplers/sobol.h"
#include "lowdiscrepancy.h"
#include "paramset.h"

namespace pbrt {

// SobolSampler Method Definitions / Sobol采样器方法定义
/**
 * @brief 计算样本的全局Sobol'序列索引
 *
 * 将当前像素坐标和样本序号映射为Sobol'序列中的全局索引。
 * 使用SobolIntervalToIndex函数实现像素空间到序列空间的转换。
 * @param sampleNum 当前像素内的样本序号
 * @return Sobol'序列中的全局索引
 */
int64_t SobolSampler::GetIndexForSample(int64_t sampleNum) const {
    return SobolIntervalToIndex(log2Resolution, sampleNum,
                                Point2i(currentPixel - sampleBounds.pMin));
}

/**
 * @brief 获取指定维度的Sobol'样本值
 *
 * 对于前两个维度（像素采样维度），将Sobol'样本值从
 * [0,1)范围映射到像素坐标范围并计算像素内偏移。
 * 其余维度直接返回归一化的Sobol'样本值。
 * @param index Sobol'序列中的全局索引
 * @param dim 采样维度
 * @return 该维度的样本值
 */
Float SobolSampler::SampleDimension(int64_t index, int dim) const {
    if (dim >= NumSobolDimensions)
        LOG(FATAL) << StringPrintf("SobolSampler can only sample up to %d "
                                   "dimensions! Exiting.",
                                   NumSobolDimensions);
    Float s = SobolSample(index, dim);
    // Remap Sobol dimensions used for pixel samples / 重映射用于像素采样的Sobol维度
    if (dim == 0 || dim == 1) {
        s = s * resolution + sampleBounds.pMin[dim];
        s = Clamp(s - currentPixel[dim], (Float)0, OneMinusEpsilon);
    }
    return s;
}

/**
 * @brief 克隆采样器
 *
 * SobolSampler的克隆不包含随机种子重新初始化，
 * 因为Sobol'序列是确定性的低差异序列。
 * @param seed 种子参数（不使用，Sobol'序列无随机性）
 * @return 新采样器的唯一指针
 */
std::unique_ptr<Sampler> SobolSampler::Clone(int seed) {
    return std::unique_ptr<Sampler>(new SobolSampler(*this));
}

/**
 * @brief 创建SobolSampler的工厂函数
 * @param params 参数集，从中读取"pixelsamples"参数
 * @param sampleBounds 采样区域边界
 * @return 创建的SobolSampler实例指针
 */
SobolSampler *CreateSobolSampler(const ParamSet &params,
                                 const Bounds2i &sampleBounds) {
    int nsamp = params.FindOneInt("pixelsamples", 16);
    if (PbrtOptions.quickRender) nsamp = 1;
    return new SobolSampler(nsamp, sampleBounds);
}

}  // namespace pbrt
