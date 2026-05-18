
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

#ifndef PBRT_SHAPES_SPHERE_H
#define PBRT_SHAPES_SPHERE_H

// shapes/sphere.h*
/**
 * @file sphere.h
 * @brief 球体(Sphere)几何体模块
 *
 * 定义了球体几何体，通过半径、Z轴裁剪范围和角度范围来定义。
 * 支持部分球体(球冠、球带)的创建，支持重要性采样。
 */
#include "shape.h"

namespace pbrt {

// Sphere Declarations / 球体声明
class Sphere : public Shape {
  public:
    // Sphere Public Methods / 球体公有方法
    Sphere(const Transform *ObjectToWorld, const Transform *WorldToObject,
           bool reverseOrientation, Float radius, Float zMin, Float zMax,
           Float phiMax)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
          radius(radius),
          zMin(Clamp(std::min(zMin, zMax), -radius, radius)),
          zMax(Clamp(std::max(zMin, zMax), -radius, radius)),
          thetaMin(std::acos(Clamp(std::min(zMin, zMax) / radius, -1, 1))),
          thetaMax(std::acos(Clamp(std::max(zMin, zMax) / radius, -1, 1))),
          phiMax(Radians(Clamp(phiMax, 0, 360))) {}
    /**
     * @brief 计算对象空间包围盒
     */
    Bounds3f ObjectBound() const;
    /**
     * @brief 光线-球体求交(完整求交)
     */
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    /**
     * @brief 光线-球体求交测试(仅判断是否相交)
     */
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    /**
     * @brief 计算球体表面积
     */
    Float Area() const;
    /**
     * @brief 在球体表面均匀采样一个点
     */
    Interaction Sample(const Point2f &u, Float *pdf) const;
    /**
     * @brief 从参考点对球体进行重要性采样
     */
    Interaction Sample(const Interaction &ref, const Point2f &u,
                       Float *pdf) const;
    /**
     * @brief 计算从参考点采样球体的概率密度
     */
    Float Pdf(const Interaction &ref, const Vector3f &wi) const;
    /**
     * @brief 计算球体相对于某点的立体角
     */
    Float SolidAngle(const Point3f &p, int nSamples) const;

  private:
    // Sphere Private Data / 球体私有数据
    const Float radius;           /**< 球体半径 */
    const Float zMin, zMax;      /**< Z轴裁剪范围 */
    const Float thetaMin, thetaMax, phiMax; /**< 角度范围参数 */
};

/**
 * @brief 创建球体形状的工厂函数
 */
std::shared_ptr<Shape> CreateSphereShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_SPHERE_H
