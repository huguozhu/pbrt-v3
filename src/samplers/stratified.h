
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

#ifndef PBRT_SAMPLERS_STRATIFIED_H
#define PBRT_SAMPLERS_STRATIFIED_H

// samplers/stratified.h*
/**
 * @file stratified.h
 * @brief 分层采样器(StratifiedSampler)模块
 *
 * 将像素区域划分为规则网格(单元)，在每个单元内独立采样。
 * 支持抖动采样(Jittering)和N-Rooks采样策略，
 * 通过空间均匀分布减少样本聚集，提高收敛速度。
 */
#include "sampler.h"
#include "rng.h"

namespace pbrt {

// StratifiedSampler Declarations / 分层采样器声明
/**
 * @brief 分层采样器
 *
 * 将像素区域细分为xPixelSamples x yPixelSamples的规则网格，
 * 在每个网格单元内生成一个样本点。当启用抖动(jitter)时，
 * 样本在单元内随机偏移，否则位于单元中心。
 * 分层采样有效减少样本聚集，改善蒙特卡洛估计的收敛速度。
 */
// �ֲ��������ͼ��ƽ�滮��Ϊ�������򣬲��ڸ����������ɵ�һ����
class StratifiedSampler : public PixelSampler {
  public:
    // StratifiedSampler Public Methods / 公有方法
    /**
     * @brief 构造函数
     * @param xPixelSamples X方向的分层数
     * @param yPixelSamples Y方向的分层数
     * @param jitterSamples 是否启用抖动偏移
     * @param nSampledDimensions 需要独立采样的维度数
     */
    StratifiedSampler(int xPixelSamples, int yPixelSamples, bool jitterSamples,
                      int nSampledDimensions)
        : PixelSampler(xPixelSamples * yPixelSamples, nSampledDimensions),
          xPixelSamples(xPixelSamples),
          yPixelSamples(yPixelSamples),
          jitterSamples(jitterSamples) {}
    /**
     * @brief 开始新像素的采样，生成分层样本
     */
    void StartPixel(const Point2i &);
    /**
     * @brief 克隆采样器
     */
    std::unique_ptr<Sampler> Clone(int seed);

  private:
    // StratifiedSampler Private Data / 私有数据
    const int xPixelSamples, yPixelSamples; /**< X和Y方向的分层数 */
	// n-rooks sampling:��nxn��С�����У�����n�������㣬ÿ����������������ж�û������������
	// jitter sampling:��nxn��С�����У�ÿ��С��������һ��������
    const bool jitterSamples; /**< 是否启用抖动偏移 */
};

/**
 * @brief 创建StratifiedSampler的工厂函数
 * @param params 参数集，从中读取"jitter"/"xsamples"/"ysamples"/"dimensions"参数
 * @return 创建的StratifiedSampler实例指针
 */
StratifiedSampler *CreateStratifiedSampler(const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SAMPLERS_STRATIFIED_H
