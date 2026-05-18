
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

// integrators/ao.cpp*
// AOIntegrator实现：通过向半球方向发射大量光线并检测遮挡来计算每个着色点的环境光照值
#include "integrators/ao.h"
#include "sampling.h"
#include "interaction.h"
#include "paramset.h"
#include "camera.h"
#include "film.h"
#include "scene.h"

namespace pbrt {

// AOIntegrator Method Definitions
// AOIntegrator构造函数：初始化采样方式和采样数量
AOIntegrator::AOIntegrator(bool cosSample, int ns,
                           std::shared_ptr<const Camera> camera,
                           std::shared_ptr<Sampler> sampler,
                           const Bounds2i &pixelBounds)
    : SamplerIntegrator(camera, sampler, pixelBounds),
      cosSample(cosSample) {
    nSamples = sampler->RoundCount(ns);
    if (ns != nSamples)
        Warning("Taking %d samples, not %d as specified", nSamples, ns);
    sampler->Request2DArray(nSamples);
}

Spectrum AOIntegrator::Li(const RayDifferential &r, const Scene &scene,
                          Sampler &sampler, MemoryArena &arena,
                          int depth) const {
    // AOIntegrator::Li：计算单个像素的环境光照值，通过采样半球方向并检测可见性来实现
    ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum L(0.f);
    RayDifferential ray(r);

    // Intersect _ray_ with scene and store intersection in _isect_
    // 与场景求交，获取表面交点
    SurfaceInteraction isect;
 retry:
    if (scene.Intersect(ray, &isect)) {
        isect.ComputeScatteringFunctions(ray, arena, true);
        if (!isect.bsdf) {
            VLOG(2) << "Skipping intersection due to null bsdf";
            ray = isect.SpawnRay(ray.d);
            goto retry;
        }

        // Compute coordinate frame based on true geometry, not shading
        // geometry.
        // 基于几何法线构建局部坐标系（而非着色法线），用于生成半球采样方向
        Normal3f n = Faceforward(isect.n, -ray.d);
        Vector3f s = Normalize(isect.dpdu);
        Vector3f t = Cross(isect.n, s);

        const Point2f *u = sampler.Get2DArray(nSamples);
        for (int i = 0; i < nSamples; ++i) {
            Vector3f wi;
            Float pdf;
            if (cosSample) {
                wi = CosineSampleHemisphere(u[i]);
                pdf = CosineHemispherePdf(std::abs(wi.z));
            } else {
                wi = UniformSampleHemisphere(u[i]);
                pdf = UniformHemispherePdf();
            }

            // Transform wi from local frame to world space.
            wi = Vector3f(s.x * wi.x + t.x * wi.y + n.x * wi.z,
                          s.y * wi.x + t.y * wi.y + n.y * wi.z,
                          s.z * wi.x + t.z * wi.y + n.z * wi.z);

            // 检查该方向是否被遮挡，未被遮挡则累积贡献值
            if (!scene.IntersectP(isect.SpawnRay(wi)))
                L += Dot(wi, n) / (pdf * nSamples);
        }
    }
    return L;
}

// CreateAOIntegrator：根据参数集创建AOIntegrator实例
AOIntegrator *CreateAOIntegrator(const ParamSet &params,
                                 std::shared_ptr<Sampler> sampler,
                                 std::shared_ptr<const Camera> camera) {
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    bool cosSample = params.FindOneBool("cossample", true);
    int nSamples = params.FindOneInt("nsamples", 64);
    if (PbrtOptions.quickRender) nSamples = 1;
    return new AOIntegrator(cosSample, nSamples, camera, sampler, pixelBounds);
}

}  // namespace pbrt
