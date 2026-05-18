
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


// core/sampling.cpp*
// Sampling: 采样函数实现，提供各种蒙特卡洛采样方法，
// 包括分层采样、拉丁超立方采样、圆盘采样、半球/球面均匀采样、
// 圆锥采样以及二维分布采样等
#include "sampling.h"
#include "geometry.h"
#include "shape.h"

namespace pbrt {

// Sampling Function Definitions
// 分层1D采样：将[0,1)区间等分为nSamples层，每层内可选抖动(jitter)或取中心
void StratifiedSample1D(Float *samp, int nSamples, RNG &rng, bool jitter) {
    Float invNSamples = (Float)1 / nSamples;
    for (int i = 0; i < nSamples; ++i) {
        Float delta = jitter ? rng.UniformFloat() : 0.5f;
        samp[i] = std::min((i + delta) * invNSamples, OneMinusEpsilon);
    }
}

// 分层2D采样：将[0,1)^2区域等分为nx×ny个网格，每个网格内可选抖动采样
void StratifiedSample2D(Point2f *samp, int nx, int ny, RNG &rng, bool jitter) {
    Float dx = (Float)1 / nx, dy = (Float)1 / ny;
    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x) {
            Float jx = jitter ? rng.UniformFloat() : 0.5f;
            Float jy = jitter ? rng.UniformFloat() : 0.5f;
            samp->x = std::min((x + jx) * dx, OneMinusEpsilon);
            samp->y = std::min((y + jy) * dy, OneMinusEpsilon);
            ++samp;
        }
}

// 拉丁超立方采样(LHS)：在每个维度上均匀分层后随机排列，保证每个维度的投影均匀覆盖
void LatinHypercube(Float *samples, int nSamples, int nDim, RNG &rng) {
    // Generate LHS samples along diagonal
    Float invNSamples = (Float)1 / nSamples;
    for (int i = 0; i < nSamples; ++i)
        for (int j = 0; j < nDim; ++j) {
            Float sj = (i + (rng.UniformFloat())) * invNSamples;
            samples[nDim * i + j] = std::min(sj, OneMinusEpsilon);
        }

    // Permute LHS samples in each dimension
    for (int i = 0; i < nDim; ++i) {
        for (int j = 0; j < nSamples; ++j) {
            int other = j + rng.UniformUInt32(nSamples - j);
            std::swap(samples[nDim * j + i], samples[nDim * other + i]);
        }
    }
}

// 拒绝采样法在单位圆盘内采样：在[-1,1]^2正方形内生成点，拒绝落在圆外的点
Point2f RejectionSampleDisk(RNG &rng) {
    Point2f p;
    do {
        p.x = 1 - 2 * rng.UniformFloat();
        p.y = 1 - 2 * rng.UniformFloat();
    } while (p.x * p.x + p.y * p.y > 1);
    return p;
}

// 在半球面上均匀采样：使用余弦方法，z = u[0] (cosθ)，φ = 2πu[1]
Vector3f UniformSampleHemisphere(const Point2f &u) {
    Float z = u[0];
    Float r = std::sqrt(std::max((Float)0, (Float)1. - z * z));
    Float phi = 2 * Pi * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

Float UniformHemispherePdf() { return Inv2Pi; }  // 半球均匀采样的PDF：1/(2π)

// 在球面上均匀采样：z = 1-2u[0] (cosθ)，φ = 2πu[1]
Vector3f UniformSampleSphere(const Point2f &u) {
    Float z = 1 - 2 * u[0];
    Float r = std::sqrt(std::max((Float)0, (Float)1 - z * z));
    Float phi = 2 * Pi * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

Float UniformSpherePdf() { return Inv4Pi; }  // 球面均匀采样的PDF：1/(4π)

// 在圆盘上均匀采样：使用极坐标映射，r = √u[0]，θ = 2πu[1]
Point2f UniformSampleDisk(const Point2f &u) {
    Float r = std::sqrt(u[0]);
    Float theta = 2 * Pi * u[1];
    return Point2f(r * std::cos(theta), r * std::sin(theta));
}

// 同心映射圆盘采样：将正方形[-1,1]^2映射到圆盘，比极坐标映射畸变更小
Point2f ConcentricSampleDisk(const Point2f &u) {
    // Map uniform random numbers to $[-1,1]^2$
    Point2f uOffset = 2.f * u - Vector2f(1, 1);

    // Handle degeneracy at the origin
    if (uOffset.x == 0 && uOffset.y == 0) return Point2f(0, 0);

    // Apply concentric mapping to point
    Float theta, r;
    if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
        r = uOffset.x;
        theta = PiOver4 * (uOffset.y / uOffset.x);
    } else {
        r = uOffset.y;
        theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
    }
    return r * Point2f(std::cos(theta), std::sin(theta));
}

// 圆锥均匀采样的PDF：1/(2π(1-cosθ_max))
Float UniformConePdf(Float cosThetaMax) {
    return 1 / (2 * Pi * (1 - cosThetaMax));
}

// 在z轴方向的圆锥内均匀采样，cosThetaMax为圆锥半角余弦
Vector3f UniformSampleCone(const Point2f &u, Float cosThetaMax) {
    Float cosTheta = ((Float)1 - u[0]) + u[0] * cosThetaMax;
    Float sinTheta = std::sqrt((Float)1 - cosTheta * cosTheta);
    Float phi = u[1] * 2 * Pi;
    return Vector3f(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta,
                    cosTheta);
}

// 在任意方向(z轴)的圆锥内均匀采样，使用给定的正交基(x,y,z)表示方向
Vector3f UniformSampleCone(const Point2f &u, Float cosThetaMax,
                           const Vector3f &x, const Vector3f &y,
                           const Vector3f &z) {
    Float cosTheta = Lerp(u[0], cosThetaMax, 1.f);
    Float sinTheta = std::sqrt((Float)1. - cosTheta * cosTheta);
    Float phi = u[1] * 2 * Pi;
    return std::cos(phi) * sinTheta * x + std::sin(phi) * sinTheta * y +
           cosTheta * z;
}

// 在三角形上均匀采样：使用重心坐标，将单位正方形映射到三角形
Point2f UniformSampleTriangle(const Point2f &u) {
    Float su0 = std::sqrt(u[0]);
    return Point2f(1 - su0, u[1] * su0);
}

// 2D分布构造：对每个v方向计算条件分布p[u|v]，再计算边缘分布p[v]，用于二维重要性采样
Distribution2D::Distribution2D(const Float *func, int nu, int nv) {
    pConditionalV.reserve(nv);
    for (int v = 0; v < nv; ++v) {
        // Compute conditional sampling distribution for $\tilde{v}$
        pConditionalV.emplace_back(new Distribution1D(&func[v * nu], nu));
    }
    // Compute marginal sampling distribution $p[\tilde{v}]$
    std::vector<Float> marginalFunc;
    marginalFunc.reserve(nv);
    for (int v = 0; v < nv; ++v)
        marginalFunc.push_back(pConditionalV[v]->funcInt);
    pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));
}

}  // namespace pbrt
