
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

#ifndef PBRT_CORE_LIGHT_H
#define PBRT_CORE_LIGHT_H

// core/light.h*
#include "pbrt.h"
#include "memory.h"
#include "interaction.h"

namespace pbrt {

// LightFlags Declarations
enum class LightFlags : int {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8
};

inline bool IsDeltaLight(int flags) {
    return flags & (int)LightFlags::DeltaPosition ||
           flags & (int)LightFlags::DeltaDirection;
}

// Light Declarations
class Light {
  public:
    // Light Interface
    virtual ~Light();
    Light(int flags, const Transform &LightToWorld,
          const MediumInterface &mediumInterface, int nSamples = 1);
    virtual Spectrum Sample_Li(const Interaction &ref, const Point2f &u,
                               Vector3f *wi, Float *pdf,
                               VisibilityTester *vis) const = 0;
    // 光源已发射的总能量(total emitted power)
	virtual Spectrum Power() const = 0;
	// 光源在使用前的预处理，大多数实现都是空值
    virtual void Preprocess(const Scene &scene) {}
    virtual Spectrum Le(const RayDifferential &r) const;
    virtual Float Pdf_Li(const Interaction &ref, const Vector3f &wi) const = 0;
    virtual Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
                               Ray *ray, Normal3f *nLight, Float *pdfPos,
                               Float *pdfDir) const = 0;
    virtual void Pdf_Le(const Ray &ray, const Normal3f &nLight, Float *pdfPos,
                        Float *pdfDir) const = 0;

    // Light Public Data
    const int flags;			// 光线类型
    const int nSamples;			// 用于area light，用于技术软阴影
    const MediumInterface mediumInterface;	// inside和outside的介质
  protected:
    // Light Protected Data
	// 光源坐标系到世界坐标系即及其反向的变换矩阵
	// 光源坐标系:光源位置位于原点(0,0,0),方向为(0,0,1)
    const Transform LightToWorld, WorldToLight;		
};

// 用于判断两点之间是否有遮挡
class VisibilityTester {
  public:
    VisibilityTester() {}
    // VisibilityTester Public Methods
    VisibilityTester(const Interaction &p0, const Interaction &p1)
        : p0(p0), p1(p1) {}
    const Interaction &P0() const { return p0; }
    const Interaction &P1() const { return p1; }
    bool Unoccluded(const Scene &scene) const;					// 在p0和p1之间是否有遮挡物
    Spectrum Tr(const Scene &scene, Sampler &sampler) const;	// 从p0透过介质传播到p1透射比

  private:
    Interaction p0, p1;
};

class AreaLight : public Light {
  public:
    // AreaLight Interface
    AreaLight(const Transform &LightToWorld, const MediumInterface &medium, int nSamples);
	// 光源在面上的交点intr处，在方向为w处的放射光(emitted radiance)
    virtual Spectrum L(const Interaction &intr, const Vector3f &w) const = 0;
};

}  // namespace pbrt

#endif  // PBRT_CORE_LIGHT_H
