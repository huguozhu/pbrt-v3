
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

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

// core/sampler.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include <inttypes.h>

namespace pbrt {

// Sampler Declarations
class Sampler {
  public:
    // Sampler Interface
    virtual ~Sampler();
    Sampler(int64_t samplesPerPixel);
    virtual void StartPixel(const Point2i &p);				// 开始采样时调用
    virtual Float Get1D() = 0;								// 返回下一个一维采样点的采样位置
    virtual Point2f Get2D() = 0;							// 返回下一个二维采样点的采样位置
    CameraSample GetCameraSample(const Point2i &pRaster);	// 
    void Request1DArray(int n);
    void Request2DArray(int n);
    virtual int RoundCount(int n) const { return n; }
    const Float *Get1DArray(int n);
    const Point2f *Get2DArray(int n);
    virtual bool StartNextSample();							// 当前像素是否还有下个采样点
    virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
    virtual bool SetSampleNumber(int64_t sampleNum);		// 设置当前采样索引，因为有些算法不需要所有采样点，所以调用该函数跳过某些采样点
    std::string StateString() const {
      return StringPrintf("(%d,%d), sample %" PRId64, currentPixel.x,
                          currentPixel.y, currentPixelSampleIndex);
    }
    int64_t CurrentSampleNumber() const { return currentPixelSampleIndex; }

    // Sampler Public Data
    const int64_t samplesPerPixel;		// 每个像素的采样个数(该值在pbrt资源文件中的sampler字段指定）

  protected:
    // Sampler Protected Data
    Point2i currentPixel;				// 当前被采样像素点的位置
    int64_t currentPixelSampleIndex;	// 当前被采样像素点的采样索引号
    std::vector<int> samples1DArraySizes, samples2DArraySizes;
    std::vector<std::vector<Float>> sampleArray1D;
    std::vector<std::vector<Point2f>> sampleArray2D;

  private:
    // Sampler Private Data
    size_t array1DOffset, array2DOffset;
};
// 一个像素的采样
class PixelSampler : public Sampler {
  public:
    // PixelSampler Public Methods
    PixelSampler(int64_t samplesPerPixel, int nSampledDimensions);
    bool StartNextSample();
    bool SetSampleNumber(int64_t);
    Float Get1D();
    Point2f Get2D();

  protected:
    // PixelSampler Protected Data
    std::vector<std::vector<Float>> samples1D;		// 一个像素中的所有采样值（float），外部的vector表示。。。
    std::vector<std::vector<Point2f>> samples2D;	// sampler1D对应的采样值在像素中的二维坐标
    int current1DDimension = 0, current2DDimension = 0;
    RNG rng;
};

class GlobalSampler : public Sampler {
  public:
    // GlobalSampler Public Methods
    bool StartNextSample();
    void StartPixel(const Point2i &);
    bool SetSampleNumber(int64_t sampleNum);
    Float Get1D();
    Point2f Get2D();
    GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) {}
    virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
    virtual Float SampleDimension(int64_t index, int dimension) const = 0;

  private:
    // GlobalSampler Private Data
    int dimension;
    int64_t intervalSampleIndex;
    static const int arrayStartDim = 5;
    int arrayEndDim;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SAMPLER_H
