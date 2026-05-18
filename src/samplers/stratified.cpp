
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


// samplers/stratified.cpp*
/**
 * @file stratified.cpp
 * @brief 分层采样器(StratifiedSampler)的实现
 *
 * 使用分层采样(Stratified Sampling)策略生成像素样本。
 * 将像素区域划分为xPixelSamples x yPixelSamples的网格，
 * 在每个单元内进行抖动采样。
 */
#include "samplers/stratified.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// StratifiedSampler Method Definitions / 分层采样器方法定义
/**
 * @brief 开始新像素的采样
 *
 * 为当前像素生成所有维度的分层样本：
 * - 1D样本使用分层采样+随机打乱(Shuffle)
 * - 2D样本使用2D分层采样+随机打乱
 * - 数组样本使用1D分层采样+打乱，或Latin Hypercube采样
 * @param p 当前像素坐标
 */
void StratifiedSampler::StartPixel(const Point2i &p) {
    ProfilePhase _(Prof::StartPixel);
    // Generate single stratified samples for the pixel / 为像素生成单层分层样本
    for (size_t i = 0; i < samples1D.size(); ++i) {
        StratifiedSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng,
                           jitterSamples);
        Shuffle(&samples1D[i][0], xPixelSamples * yPixelSamples, 1, rng);
    }
    for (size_t i = 0; i < samples2D.size(); ++i) {
        StratifiedSample2D(&samples2D[i][0], xPixelSamples, yPixelSamples, rng,
                           jitterSamples);
        Shuffle(&samples2D[i][0], xPixelSamples * yPixelSamples, 1, rng);
    }

    // Generate arrays of stratified samples for the pixel / 为像素生成数组形式的分层样本
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i)
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            int count = samples1DArraySizes[i];
            StratifiedSample1D(&sampleArray1D[i][j * count], count, rng,
                               jitterSamples);
            Shuffle(&sampleArray1D[i][j * count], count, 1, rng);
        }
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i)
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            int count = samples2DArraySizes[i];
            // 对于2D数组样本使用Latin Hypercube采样以保持良好分布
            // Use Latin Hypercube sampling for 2D array samples
            LatinHypercube(&sampleArray2D[i][j * count].x, count, 2, rng);
        }
    PixelSampler::StartPixel(p);
}

/**
 * @brief 克隆采样器
 *
 * 创建当前采样器的副本，为新副本设置独立的随机数序列。
 * @param seed 新采样器的随机数种子
 * @return 新采样器的唯一指针
 */
std::unique_ptr<Sampler> StratifiedSampler::Clone(int seed) {
    StratifiedSampler *ss = new StratifiedSampler(*this);
    ss->rng.SetSequence(seed);
    return std::unique_ptr<Sampler>(ss);
}

/**
 * @brief 创建StratifiedSampler的工厂函数
 * @param params 参数集，从中读取采样参数
 * @return 创建的StratifiedSampler实例指针
 */
StratifiedSampler *CreateStratifiedSampler(const ParamSet &params) {
    bool jitter = params.FindOneBool("jitter", true);
    int xsamp = params.FindOneInt("xsamples", 4);
    int ysamp = params.FindOneInt("ysamples", 4);
    int sd = params.FindOneInt("dimensions", 4);
    if (PbrtOptions.quickRender) xsamp = ysamp = 1;
    return new StratifiedSampler(xsamp, ysamp, jitter, sd);
}

}  // namespace pbrt
