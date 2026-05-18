
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


// core/sampler.cpp*
// Sampler: 采样器实现，提供像素采样和维度采样的基础设施，
// 包含Sampler基类、PixelSampler(像素采样器，预计算所有维度的样本)、
// GlobalSampler(全局采样器，支持按需计算任意维度的样本值)
#include "sampler.h"
#include "sampling.h"
#include "camera.h"
#include "stats.h"

namespace pbrt {

// Sampler Method Definitions
Sampler::~Sampler() {}  // Sampler虚析构函数

Sampler::Sampler(int64_t samplesPerPixel) : samplesPerPixel(samplesPerPixel) {}
// 获取相机样本：为当前像素生成包含胶片位置、时间和透镜位置的相机采样点
CameraSample Sampler::GetCameraSample(const Point2i &pRaster) {
    CameraSample cs;
    cs.pFilm = (Point2f)pRaster + Get2D();
    cs.time = Get1D();
    cs.pLens = Get2D();
    return cs;
}

// 开始处理新像素：重置样本索引和数组偏移量
void Sampler::StartPixel(const Point2i &p) {
    currentPixel = p;
    currentPixelSampleIndex = 0;
    // Reset array offsets for next pixel sample
    array1DOffset = array2DOffset = 0;
}

// 开始下一个样本：递增样本索引，返回是否还有更多样本
bool Sampler::StartNextSample() {
    // Reset array offsets for next pixel sample
    array1DOffset = array2DOffset = 0;
    return ++currentPixelSampleIndex < samplesPerPixel;
}

// 设置当前样本编号：用于随机访问特定样本
bool Sampler::SetSampleNumber(int64_t sampleNum) {
    // Reset array offsets for next pixel sample
    array1DOffset = array2DOffset = 0;
    currentPixelSampleIndex = sampleNum;
    return currentPixelSampleIndex < samplesPerPixel;
}

// 请求1D数组样本：注册一个1D样本数组维度
void Sampler::Request1DArray(int n) {
    CHECK_EQ(RoundCount(n), n);
    samples1DArraySizes.push_back(n);
    sampleArray1D.push_back(std::vector<Float>(n * samplesPerPixel));
}

// 请求2D数组样本：注册一个2D样本数组维度
void Sampler::Request2DArray(int n) {
    CHECK_EQ(RoundCount(n), n);
    samples2DArraySizes.push_back(n);
    sampleArray2D.push_back(std::vector<Point2f>(n * samplesPerPixel));
}

// 获取1D数组样本：返回已注册的1D样本数组中的下一个数组
const Float *Sampler::Get1DArray(int n) {
    if (array1DOffset == sampleArray1D.size()) return nullptr;
    CHECK_EQ(samples1DArraySizes[array1DOffset], n);
    CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
    return &sampleArray1D[array1DOffset++][currentPixelSampleIndex * n];
}

// 获取2D数组样本：返回已注册的2D样本数组中的下一个数组
const Point2f *Sampler::Get2DArray(int n) {
    if (array2DOffset == sampleArray2D.size()) return nullptr;
    CHECK_EQ(samples2DArraySizes[array2DOffset], n);
    CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
    return &sampleArray2D[array2DOffset++][currentPixelSampleIndex * n];
}

// PixelSampler构造：预分配指定数量的1D和2D样本维度
PixelSampler::PixelSampler(int64_t samplesPerPixel, int nSampledDimensions)
    : Sampler(samplesPerPixel) {
    for (int i = 0; i < nSampledDimensions; ++i) {
        samples1D.push_back(std::vector<Float>(samplesPerPixel));
        samples2D.push_back(std::vector<Point2f>(samplesPerPixel));
    }
}

// 开始PixelSampler的下一个样本：重置维度计数
bool PixelSampler::StartNextSample() {
    current1DDimension = current2DDimension = 0;
    return Sampler::StartNextSample();
}

// 设置PixelSampler的当前样本编号：重置维度计数
bool PixelSampler::SetSampleNumber(int64_t sampleNum) {
    current1DDimension = current2DDimension = 0;
    return Sampler::SetSampleNumber(sampleNum);
}

// PixelSampler获取1D样本：优先使用预计算的样本，超出则用随机数
Float PixelSampler::Get1D() {
    ProfilePhase _(Prof::GetSample);
    CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
    if (current1DDimension < samples1D.size())
        return samples1D[current1DDimension++][currentPixelSampleIndex];
    else
        return rng.UniformFloat();
}

// PixelSampler获取2D样本：优先使用预计算的样本，超出则用随机数
Point2f PixelSampler::Get2D() {
    ProfilePhase _(Prof::GetSample);
    CHECK_LT(currentPixelSampleIndex, samplesPerPixel);
    if (current2DDimension < samples2D.size())
        return samples2D[current2DDimension++][currentPixelSampleIndex];
    else
        return Point2f(rng.UniformFloat(), rng.UniformFloat());
}

// 开始处理GlobalSampler的新像素：计算所有数组样本的维度值
void GlobalSampler::StartPixel(const Point2i &p) {
    ProfilePhase _(Prof::StartPixel);
    Sampler::StartPixel(p);
    dimension = 0;
    intervalSampleIndex = GetIndexForSample(0);
    // Compute _arrayEndDim_ for dimensions used for array samples
    arrayEndDim =
        arrayStartDim + sampleArray1D.size() + 2 * sampleArray2D.size();

    // Compute 1D array samples for _GlobalSampler_
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
        int nSamples = samples1DArraySizes[i] * samplesPerPixel;
        for (int j = 0; j < nSamples; ++j) {
            int64_t index = GetIndexForSample(j);
            sampleArray1D[i][j] = SampleDimension(index, arrayStartDim + i);
        }
    }

    // Compute 2D array samples for _GlobalSampler_
    int dim = arrayStartDim + samples1DArraySizes.size();
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
        int nSamples = samples2DArraySizes[i] * samplesPerPixel;
        for (int j = 0; j < nSamples; ++j) {
            int64_t idx = GetIndexForSample(j);
            sampleArray2D[i][j].x = SampleDimension(idx, dim);
            sampleArray2D[i][j].y = SampleDimension(idx, dim + 1);
        }
        dim += 2;
    }
    CHECK_EQ(arrayEndDim, dim);
}

// 开始GlobalSampler的下一个样本：重置维度并更新样本索引
bool GlobalSampler::StartNextSample() {
    dimension = 0;
    intervalSampleIndex = GetIndexForSample(currentPixelSampleIndex + 1);
    return Sampler::StartNextSample();
}

// 设置GlobalSampler的当前样本编号
bool GlobalSampler::SetSampleNumber(int64_t sampleNum) {
    dimension = 0;
    intervalSampleIndex = GetIndexForSample(sampleNum);
    return Sampler::SetSampleNumber(sampleNum);
}

// GlobalSampler获取1D样本：通过SampleDimension计算指定维度的样本值
Float GlobalSampler::Get1D() {
    ProfilePhase _(Prof::GetSample);
    if (dimension >= arrayStartDim && dimension < arrayEndDim)
        dimension = arrayEndDim;
    return SampleDimension(intervalSampleIndex, dimension++);
}

// GlobalSampler获取2D样本：连续获取两个维度的样本值
Point2f GlobalSampler::Get2D() {
    ProfilePhase _(Prof::GetSample);
    if (dimension + 1 >= arrayStartDim && dimension < arrayEndDim)
        dimension = arrayEndDim;
    Point2f p(SampleDimension(intervalSampleIndex, dimension),
              SampleDimension(intervalSampleIndex, dimension + 1));
    dimension += 2;
    return p;
}

}  // namespace pbrt
