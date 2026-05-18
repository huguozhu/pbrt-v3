
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


// samplers/halton.cpp*
/**
 * @file halton.cpp
 * @brief Halton采样器的实现
 *
 * 实现了Halton低差异序列采样器，使用Radical Inverse(基数反转)生成
 * 在空间中均匀分布的样本。通过扩展欧几里得算法计算乘法逆元，
 * 使用可选的置换(Scrambling)提高采样质量。
 */
#include "samplers/halton.h"
#include "paramset.h"
#include "rng.h"

namespace pbrt {

// HaltonSampler Local Constants / Halton采样器局部常量
static PBRT_CONSTEXPR int kMaxResolution = 128;

// HaltonSampler Utility Functions / Halton采样器工具函数
static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y);  /**< 扩展欧几里得算法声明 */
/**
 * @brief 计算模n下的乘法逆元
 * 使用扩展欧几里得算法
 */
static uint64_t multiplicativeInverse(int64_t a, int64_t n) {
    int64_t x, y;
    extendedGCD(a, n, &x, &y);
    return Mod(x, n);
}

/**
 * @brief 扩展欧几里得算法
 * 计算gcd(a,b)以及系数x,y使得 a*x + b*y = gcd(a,b)
 */
static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y) {
    if (b == 0) {
        *x = 1;
        *y = 0;
        return;
    }
    int64_t d = a / b, xp, yp;
    extendedGCD(b, a % b, &xp, &yp);
    *x = yp;
    *y = xp - (d * yp);
}

// HaltonSampler Method Definitions / Halton采样器方法实现
/**
 * @brief Halton采样器构造函数
 *
 * 预计算Radical Inverse置换表，确定基数和指数以覆盖采样区域，
 * 计算像素间样本步长和乘法逆元。
 */
HaltonSampler::HaltonSampler(int samplesPerPixel, const Bounds2i &sampleBounds,
                             bool sampleAtPixelCenter)
    : GlobalSampler(samplesPerPixel), sampleAtPixelCenter(sampleAtPixelCenter) {
    // Generate random digit permutations for Halton sampler / 生成Halton采样器的随机数位置换
    if (radicalInversePermutations.empty()) {
        RNG rng;
        radicalInversePermutations = ComputeRadicalInversePermutations(rng);
    }

    // Find radical inverse base scales and exponents that cover sampling area / 计算覆盖采样区域的基数和指数
    Vector2i res = sampleBounds.pMax - sampleBounds.pMin;
    for (int i = 0; i < 2; ++i) {
        int base = (i == 0) ? 2 : 3;
        int scale = 1, exp = 0;
        while (scale < std::min(res[i], kMaxResolution)) {
            scale *= base;
            ++exp;
        }
        baseScales[i] = scale;
        baseExponents[i] = exp;
    }

    // Compute stride in samples for visiting each pixel area / 计算访问每个像素区域的样本步长
    sampleStride = baseScales[0] * baseScales[1];

    // Compute multiplicative inverses for _baseScales_ / 计算基数之间的乘法逆元
    multInverse[0] = multiplicativeInverse(baseScales[1], baseScales[0]);
    multInverse[1] = multiplicativeInverse(baseScales[0], baseScales[1]);
}

std::vector<uint16_t> HaltonSampler::radicalInversePermutations;
/**
 * @brief 计算样本的全局索引
 * 根据当前像素位置和样本序号，计算Halton序列中的全局索引
 */
int64_t HaltonSampler::GetIndexForSample(int64_t sampleNum) const {
    if (currentPixel != pixelForOffset) {
        // Compute Halton sample offset for _currentPixel_ / 计算当前像素的Halton样本偏移量
        offsetForCurrentPixel = 0;
        if (sampleStride > 1) {
            Point2i pm(Mod(currentPixel[0], kMaxResolution),
                       Mod(currentPixel[1], kMaxResolution));
            for (int i = 0; i < 2; ++i) {
                uint64_t dimOffset =
                    (i == 0)
                        ? InverseRadicalInverse<2>(pm[i], baseExponents[i])
                        : InverseRadicalInverse<3>(pm[i], baseExponents[i]);
                offsetForCurrentPixel +=
                    dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
            }
            offsetForCurrentPixel %= sampleStride;
        }
        pixelForOffset = currentPixel;
    }
    return offsetForCurrentPixel + sampleNum * sampleStride;
}

/**
 * @brief 获取指定维度的样本值
 * 对前两维使用Radical Inverse，其它维度使用加扰Radical Inverse
 */
Float HaltonSampler::SampleDimension(int64_t index, int dim) const {
    if (sampleAtPixelCenter && (dim == 0 || dim == 1)) return 0.5f;
    if (dim == 0)
        return RadicalInverse(dim, index >> baseExponents[0]);
    else if (dim == 1)
        return RadicalInverse(dim, index / baseScales[1]);
    else
        return ScrambledRadicalInverse(dim, index,
                                       PermutationForDimension(dim));
}

/**
 * @brief 克隆Halton采样器
 */
std::unique_ptr<Sampler> HaltonSampler::Clone(int seed) {
    return std::unique_ptr<Sampler>(new HaltonSampler(*this));
}

/**
 * @brief 创建Halton采样器的工厂函数
 */
HaltonSampler *CreateHaltonSampler(const ParamSet &params,
                                   const Bounds2i &sampleBounds) {
    int nsamp = params.FindOneInt("pixelsamples", 16);
    if (PbrtOptions.quickRender) nsamp = 1;
    bool sampleAtCenter = params.FindOneBool("samplepixelcenter", false);
    return new HaltonSampler(nsamp, sampleBounds, sampleAtCenter);
}

}  // namespace pbrt
