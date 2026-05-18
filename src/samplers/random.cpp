
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


// samplers/random.cpp*
/**
 * @file random.cpp
 * @brief 纯随机采样器(RandomSampler)的实现
 *
 * 提供基于伪随机数生成器的完全随机采样策略。
 * 每个像素的样本通过RNG独立均匀生成。
 */
#include "samplers/random.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

/**
 * @brief 构造函数，初始化基类Sampler和随机数生成器
 * @param ns 每像素样本数
 * @param seed 随机数种子
 */
RandomSampler::RandomSampler(int ns, int seed) : Sampler(ns), rng(seed) {}

/**
 * @brief 获取下一个1D随机样本值
 *
 * 从RNG中直接获取均匀分布的随机浮点数。
 * @return [0,1)范围内的随机浮点数
 */
Float RandomSampler::Get1D() {
    ProfilePhase _(Prof::GetSample);
    CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
    return rng.UniformFloat();
}

/**
 * @brief 获取下一个2D随机样本值
 *
 * 从RNG中独立生成两个均匀随机数组成2D样本点。
 * @return [0,1)x[0,1)范围内的随机2D点
 */
Point2f RandomSampler::Get2D() {
    ProfilePhase _(Prof::GetSample);
    CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
    return {rng.UniformFloat(), rng.UniformFloat()};
}

/**
 * @brief 克隆采样器
 *
 * 创建当前采样器的副本，并为新副本设置独立的随机数序列。
 * @param seed 新采样器的随机数种子
 * @return 新采样器的唯一指针
 */
std::unique_ptr<Sampler> RandomSampler::Clone(int seed) {
    RandomSampler *rs = new RandomSampler(*this);
    rs->rng.SetSequence(seed);
    return std::unique_ptr<Sampler>(rs);
}

/**
 * @brief 开始新像素的采样
 *
 * 预生成该像素所有数组维度所需的随机样本。
 * 1D和2D数组样本均从RNG中直接获取。
 * @param p 当前像素坐标
 */
void RandomSampler::StartPixel(const Point2i &p) {
    ProfilePhase _(Prof::StartPixel);
    // 生成1D数组样本 / Generate 1D array samples
    for (size_t i = 0; i < sampleArray1D.size(); ++i)
        for (size_t j = 0; j < sampleArray1D[i].size(); ++j)
            sampleArray1D[i][j] = rng.UniformFloat();

    // 生成2D数组样本 / Generate 2D array samples
    for (size_t i = 0; i < sampleArray2D.size(); ++i)
        for (size_t j = 0; j < sampleArray2D[i].size(); ++j)
            sampleArray2D[i][j] = {rng.UniformFloat(), rng.UniformFloat()};
    Sampler::StartPixel(p);
}

/**
 * @brief 创建RandomSampler的工厂函数
 * @param params 参数集
 * @return 创建的RandomSampler实例指针
 */
Sampler *CreateRandomSampler(const ParamSet &params) {
    int ns = params.FindOneInt("pixelsamples", 4);
    return new RandomSampler(ns);
}

}  // namespace pbrt
