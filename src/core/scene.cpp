
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


// core/scene.cpp*
// Scene: 场景类实现，提供光线与场景图元的求交(Intersect)、遮挡查询(IntersectP)
// 以及透射率计算(IntersectTr)接口
#include "scene.h"
#include "stats.h"

namespace pbrt {

STAT_COUNTER("Intersections/Regular ray intersection tests",
             nIntersectionTests);
STAT_COUNTER("Intersections/Shadow ray intersection tests", nShadowTests);

// Scene Method Definitions
// Intersect: 对光线与场景求交，返回最近的交点信息，统计光线求交测试次数
bool Scene::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ++nIntersectionTests;
    DCHECK_NE(ray.d, Vector3f(0,0,0));
    return aggregate->Intersect(ray, isect);
}

// IntersectP: 遮挡查询，只判断光线是否与场景相交（不返回交点信息），统计阴影测试次数
bool Scene::IntersectP(const Ray &ray) const {
    ++nShadowTests;
    DCHECK_NE(ray.d, Vector3f(0,0,0));
    return aggregate->IntersectP(ray);
}

// IntersectTr: 计算光线在场景中的透射率，考虑参与介质的衰减，
// 如果遇到不透明表面则返回该表面，否则继续传播直到离开场景
bool Scene::IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *isect,
                        Spectrum *Tr) const {
    *Tr = Spectrum(1.f);
    while (true) {
        bool hitSurface = Intersect(ray, isect);
        // Accumulate beam transmittance for ray segment
        if (ray.medium) *Tr *= ray.medium->Tr(ray, sampler);

        // Initialize next ray segment or terminate transmittance computation
        if (!hitSurface) return false;
        if (isect->primitive->GetMaterial() != nullptr) return true;
        ray = isect->SpawnRay(ray.d);
    }
}

}  // namespace pbrt
