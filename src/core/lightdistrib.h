
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

#ifndef PBRT_CORE_LIGHTDISTRIB_H
#define PBRT_CORE_LIGHTDISTRIB_H

// core/lightdistrib.h*
// 光源分布: 提供光源重要性采样的概率分布，支持均匀分布、功率分布和空间变化分布
#include "pbrt.h"
#include "geometry.h"
#include "sampling.h"
#include <atomic>
#include <functional>
#include <mutex>
#include <unordered_map>
#include <vector>

namespace pbrt {

// LightDistribution 定义了一个通用接口，用于为空间中给定点提供光源采样的概率分布。
// 光源采样分布接口：提供空间中某点处光源采样概率分布的通用接口
class LightDistribution {
  public:
    virtual ~LightDistribution();

    // Lookup: 返回空间中一点p处用于光源采样的概率分布
    // Given a point |p| in space, this method returns a (hopefully
    // effective) sampling distribution for light sources at that point.
    virtual const Distribution1D *Lookup(const Point3f &p) const = 0;
};

// CreateLightSampleDistribution: 根据名称创建相应的光源采样分布（均匀/功率/空间变化）
std::unique_ptr<LightDistribution> CreateLightSampleDistribution(
    const std::string &name, const Scene &scene);

// UniformLightDistribution: 最简单的光源分布实现，对所有光源返回均匀分布（忽略位置点p）。
// 适用于简单场景，但光源较多时效率较低。
// The simplest possible implementation of LightDistribution: this returns
// a uniform distribution over all light sources, ignoring the provided
// point. This approach works well for very simple scenes, but is quite
// ineffective for scenes with more than a handful of light sources. (This
// was the sampling method originally used for the PathIntegrator and the
// VolPathIntegrator in the printed book, though without the
// UniformLightDistribution class.)
class UniformLightDistribution : public LightDistribution {
  public:
    UniformLightDistribution(const Scene &scene);
    const Distribution1D *Lookup(const Point3f &p) const;

  private:
    std::unique_ptr<Distribution1D> distrib;
};

// PowerLightDistribution: 按光源功率分配采样概率（忽略位置p），功率越大的光源被采样的概率越高。
// 适用于主要光照来自少数强光源的场景。
// PowerLightDistribution returns a distribution with sampling probability
// proportional to the total emitted power for each light. (It also ignores
// the provided point |p|.)  This approach works well for scenes where
// there the most powerful lights are also the most important contributors
// to lighting in the scene, but doesn't do well if there are many lights
// and if different lights are relatively important in some areas of the
// scene and unimportant in others. (This was the default sampling method
// used for the BDPT integrator and MLT integrator in the printed book,
// though also without the PowerLightDistribution class.)
class PowerLightDistribution : public LightDistribution {
  public:
    PowerLightDistribution(const Scene &scene);
    const Distribution1D *Lookup(const Point3f &p) const;

  private:
    std::unique_ptr<Distribution1D> distrib;
};

// SpatialLightDistribution: 空间变化的光源分布，根据光源对空间中不同区域的贡献调整采样概率。
// 将场景划分为体素网格，为每个体素单独计算采样分布。
// A spatially-varying light distribution that adjusts the probability of
// sampling a light source based on an estimate of its contribution to a
// region of space.  A fixed voxel grid is imposed over the scene bounds
// and a sampling distribution is computed as needed for each voxel.
class SpatialLightDistribution : public LightDistribution {
  public:
    SpatialLightDistribution(const Scene &scene, int maxVoxels = 64);
    ~SpatialLightDistribution();
    const Distribution1D *Lookup(const Point3f &p) const;

  private:
    // Compute the sampling distribution for the voxel with integer
    // coordiantes given by "pi".
    Distribution1D *ComputeDistribution(Point3i pi) const;

    const Scene &scene;
    int nVoxels[3];

    // The hash table is a fixed number of HashEntry structs (where we
    // allocate more than enough entries in the SpatialLightDistribution
    // constructor). During rendering, the table is allocated without
    // locks, using atomic operations. (See the Lookup() method
    // implementation for details.)
    struct HashEntry {
        std::atomic<uint64_t> packedPos;
        std::atomic<Distribution1D *> distribution;
    };
    mutable std::unique_ptr<HashEntry[]> hashTable;
    size_t hashTableSize;
};

}  // namespace pbrt

#endif  // PBRT_CORE_LIGHTDISTRIB_H
