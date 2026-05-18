
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

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

// core/shape.h*
// Shape: 几何形状的抽象基类，定义所有几何体（球体、三角形、圆柱等）的公共接口，
// 包括求交、包围盒、面积计算和表面采样
#include "pbrt.h"
#include "geometry.h"
#include "interaction.h"
#include "memory.h"
#include "transform.h"

namespace pbrt {

// Shape Declarations
class Shape {
  public:
    // Shape Interface
    Shape(const Transform *ObjectToWorld, const Transform *WorldToObject,
          bool reverseOrientation);
    virtual ~Shape();
    virtual Bounds3f ObjectBound() const = 0;						// 物体坐标系中的包围盒
    virtual Bounds3f WorldBound() const;							// 世界坐标系中的包围盒
    virtual bool Intersect(const Ray &ray, Float *tHit,				// 判断与ray是否相交，填充相交信息到isect
                           SurfaceInteraction *isect,
                           bool testAlphaTexture = true) const = 0;
    virtual bool IntersectP(const Ray &ray,							// 只判断是否相交（不记录详细信息）
                            bool testAlphaTexture = true) const {
        return Intersect(ray, nullptr, nullptr, testAlphaTexture);
    }
    virtual Float Area() const = 0;		// 计算表面积，用于将Shape用作面光源时的采样
    // 在形状表面采样一个点，返回相对于面积的概率密度
    // Sample a point on the surface of the shape and return the PDF with
    // respect to area on the surface.
    virtual Interaction Sample(const Point2f &u, Float *pdf) const = 0;
    virtual Float Pdf(const Interaction &) const { return 1 / Area(); }

    // 给定参考点，在形状上采样一个点，返回相对于参考点立体角的概率密度
    // Sample a point on the shape given a reference point |ref| and
    // return the PDF with respect to solid angle from |ref|.
    virtual Interaction Sample(const Interaction &ref, const Point2f &u,
                               Float *pdf) const;
    virtual Float Pdf(const Interaction &ref, const Vector3f &wi) const;

    // 计算形状相对于参考点p所张的立体角（世界坐标系）
    // 某些形状可解析计算，默认实现使用蒙特卡洛积分
    // Returns the solid angle subtended by the shape w.r.t. the reference
    // point p, given in world space. Some shapes compute this value in
    // closed-form, while the default implementation uses Monte Carlo
    // integration; the nSamples parameter determines how many samples are
    // used in this case.
    virtual Float SolidAngle(const Point3f &p, int nSamples = 512) const;

    // Shape Public Data
    const Transform *ObjectToWorld, *WorldToObject;	// 物体到世界、世界到物体的变换矩阵
    const bool reverseOrientation;					// 法线方向是否翻转，用于确定shape的外部/内部朝向，发光体只在法线方向发射光线
    const bool transformSwapsHandedness;			// 默认使用左手坐标系，如果当前shape使用右手坐标系，则该值为true
};

}  // namespace pbrt

#endif  // PBRT_CORE_SHAPE_H
