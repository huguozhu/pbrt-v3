
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

#ifndef PBRT_SAMPLERS_RANDOM_H
#define PBRT_SAMPLERS_RANDOM_H

// samplers/random.h*
/**
 * @file random.h
 * @brief 纯随机采样器(RandomSampler)模块
 *
 * 使用伪随机数生成器(RNG)生成完全随机的样本。
 * 这是最简单的采样策略，每个样本独立均匀随机分布。
 */
#include "sampler.h"
#include "rng.h"

namespace pbrt {

/**
 * @brief 纯随机采样器
 *
 * 使用伪随机数生成器(RNG)生成完全随机的样本点。
 * 虽然实现简单，但由于样本聚集(clumping)效应，收敛速度较慢。
 * 适用于快速预览或对采样质量要求不高的场景。
 */
class RandomSampler : public Sampler {
  public:
    /**
     * @brief 构造函数
     * @param ns 每像素样本数
     * @param seed 随机数种子
     */
    RandomSampler(int ns, int seed = 0);
    /**
     * @brief 开始新像素的采样，预生成该像素所需的数组样本
     */
    void StartPixel(const Point2i &);
    /**
     * @brief 获取下一个1D随机样本值
     * @return [0,1)范围内的均匀随机浮点数
     */
    Float Get1D();
    /**
     * @brief 获取下一个2D随机样本值
     * @return [0,1)x[0,1)范围内的均匀随机2D点
     */
    Point2f Get2D();
    /**
     * @brief 克隆采样器，为新线程创建独立的副本
     * @param seed 新副本的随机数种子
     */
    std::unique_ptr<Sampler> Clone(int seed);

  private:
    RNG rng; /**< 伪随机数生成器 */
};

/**
 * @brief 创建RandomSampler的工厂函数
 * @param params 参数集，从中读取"pixelsamples"参数
 * @return 创建的RandomSampler实例指针
 */
Sampler *CreateRandomSampler(const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SAMPLERS_RANDOM_H
