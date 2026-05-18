
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


// lights/distant.cpp*
// DistantLight实现：平行光，所有光线沿固定方向从无限远处照射，用于模拟太阳光等方向性光源
#include "lights/distant.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// DistantLight Method Definitions
// DistantLight构造函数：初始化光辐射度和光照方向
DistantLight::DistantLight(const Transform &LightToWorld, const Spectrum &L,
                           const Vector3f &wLight)
    : Light((int)LightFlags::DeltaDirection, LightToWorld, MediumInterface()),
      L(L),
      wLight(Normalize(LightToWorld(wLight))) {}

// Sample_Li：在场景边界球外采样一个光源点，返回沿固定方向入射的辐射度
Spectrum DistantLight::Sample_Li(const Interaction &ref, const Point2f &u,
                                 Vector3f *wi, Float *pdf,
                                 VisibilityTester *vis) const {
    ProfilePhase _(Prof::LightSample);
    *wi = wLight;
    *pdf = 1;
    Point3f pOutside = ref.p + wLight * (2 * worldRadius);
    *vis =
        VisibilityTester(ref, Interaction(pOutside, ref.time, mediumInterface));
    return L;
}

// Power：计算平行光的总功率
Spectrum DistantLight::Power() const {
    return L * Pi * worldRadius * worldRadius;
}

// Pdf_Li：平行光的PDF为0（Delta方向分布，不能通过连续采样得到）
Float DistantLight::Pdf_Li(const Interaction &, const Vector3f &) const {
    return 0.f;
}

// Sample_Le：从平行光发射光线，在垂直光方向的圆盘上采样起点
Spectrum DistantLight::Sample_Le(const Point2f &u1, const Point2f &u2,
                                 Float time, Ray *ray, Normal3f *nLight,
                                 Float *pdfPos, Float *pdfDir) const {
    ProfilePhase _(Prof::LightSample);
    // Choose point on disk oriented toward infinite light direction
    Vector3f v1, v2;
    CoordinateSystem(wLight, &v1, &v2);
    Point2f cd = ConcentricSampleDisk(u1);
    Point3f pDisk = worldCenter + worldRadius * (cd.x * v1 + cd.y * v2);

    // Set ray origin and direction for infinite light ray
    *ray = Ray(pDisk + worldRadius * wLight, -wLight, Infinity, time);
    *nLight = (Normal3f)ray->d;
    *pdfPos = 1 / (Pi * worldRadius * worldRadius);
    *pdfDir = 1;
    return L;
}

void DistantLight::Pdf_Le(const Ray &, const Normal3f &, Float *pdfPos,
                          Float *pdfDir) const {
    ProfilePhase _(Prof::LightPdf);
    *pdfPos = 1 / (Pi * worldRadius * worldRadius);
    *pdfDir = 0;
}

std::shared_ptr<DistantLight> CreateDistantLight(const Transform &light2world,
                                                 const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    Point3f from = paramSet.FindOnePoint3f("from", Point3f(0, 0, 0));
    Point3f to = paramSet.FindOnePoint3f("to", Point3f(0, 0, 1));
    Vector3f dir = from - to;
    return std::make_shared<DistantLight>(light2world, L * sc, dir);
}

}  // namespace pbrt
