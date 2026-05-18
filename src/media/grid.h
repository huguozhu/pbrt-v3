
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

#ifndef PBRT_MEDIA_GRID_H
#define PBRT_MEDIA_GRID_H

// media/grid.h*
// 文件描述: 网格密度介质(体积雾/云)的头文件。基于三维密度网格的参与介质，
// 使用三线性插值计算空间中的密度值，支持光线步进和delta-tracking采样。
#include "medium.h"
#include "transform.h"
#include "stats.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Volume density grid", densityBytes);

// GridDensityMedium Declarations
// 网格密度介质类，使用三维体素网格存储密度值。
// 通过WorldToMedium变换将光线转换到介质空间进行采样，
// 采用delta-tracking(比率追踪)方法计算透射率和散射事件。
class GridDensityMedium : public Medium {
  public:
    // GridDensityMedium Public Methods
    GridDensityMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g,
                      int nx, int ny, int nz, const Transform &mediumToWorld,
                      const Float *d)
        : sigma_a(sigma_a),
          sigma_s(sigma_s),
          g(g),
          nx(nx),
          ny(ny),
          nz(nz),
          WorldToMedium(Inverse(mediumToWorld)),
          density(new Float[nx * ny * nz]) {
        densityBytes += nx * ny * nz * sizeof(Float);
        memcpy((Float *)density.get(), d, sizeof(Float) * nx * ny * nz);
        // Precompute values for Monte Carlo sampling of _GridDensityMedium_
        sigma_t = (sigma_a + sigma_s)[0];
        if (Spectrum(sigma_t) != sigma_a + sigma_s)
            Error(
                "GridDensityMedium requires a spectrally uniform attenuation "
                "coefficient!");
        Float maxDensity = 0;
        for (int i = 0; i < nx * ny * nz; ++i)
            maxDensity = std::max(maxDensity, density[i]);
        invMaxDensity = 1 / maxDensity;
    }

    // 计算给定位置p处的密度值(使用三线性插值)
    Float Density(const Point3f &p) const;
    Float D(const Point3i &p) const {
        Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
        if (!InsideExclusive(p, sampleBounds)) return 0;
        return density[(p.z * ny + p.y) * nx + p.x];
    }
    // 采样介质中的散射事件(使用delta-tracking)
    Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena,
                    MediumInteraction *mi) const;
    // 计算光线穿过介质的透射率(使用比率追踪)
    Spectrum Tr(const Ray &ray, Sampler &sampler) const;

  private:
    // GridDensityMedium Private Data
    const Spectrum sigma_a, sigma_s;  // 吸收和散射系数
    const Float g;                     // 相位函数各向异性参数
    const int nx, ny, nz;              // 体素网格维度
    const Transform WorldToMedium;     // 世界到介质空间的变换
    std::unique_ptr<Float[]> density;  // 密度值数组
    Float sigma_t;                     // 总衰减系数(sigma_a+sigma_s)
    Float invMaxDensity;               // 最大密度的倒数
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_GRID_H
