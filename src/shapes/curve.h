
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

#ifndef PBRT_SHAPES_CURVE_H
#define PBRT_SHAPES_CURVE_H

// shapes/curve.h*
/**
 * @file curve.h
 * @brief 曲线(Curve)几何体模块
 *
 * 定义了贝塞尔曲线几何体，支持三种曲线类型：Flat(平面)、Cylinder(圆柱)和Ribbon(带形)。
 * 曲线由四个控制点定义的贝塞尔曲线和宽度参数来确定。
 */
#include "shape.h"

namespace pbrt {
struct CurveCommon;

// CurveType Declarations / 曲线类型声明
enum class CurveType { Flat, Cylinder, Ribbon };

// CurveCommon Declarations / 曲线共享数据结构声明
struct CurveCommon {
    CurveCommon(const Point3f c[4], Float w0, Float w1, CurveType type,
                const Normal3f *norm);
    const CurveType type;   /**< 曲线类型 */
    Point3f cpObj[4];       /**< 对象空间中的4个贝塞尔控制点 */
    Float width[2];         /**< 曲线端点宽度 */
    Normal3f n[2];          /**< 法线(仅Ribbon类型使用) */
    Float normalAngle, invSinNormalAngle;  /**< 法线角度参数 */
};

// Curve Declarations / 曲线声明
class Curve : public Shape {
  public:
    // Curve Public Methods / 曲线公有方法
    Curve(const Transform *ObjectToWorld, const Transform *WorldToObject,
          bool reverseOrientation, const std::shared_ptr<CurveCommon> &common,
          Float uMin, Float uMax)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
          common(common),
          uMin(uMin),
          uMax(uMax) {}
    /**
     * @brief 计算对象空间包围盒
     */
    Bounds3f ObjectBound() const;
    /**
     * @brief 光线-曲线求交测试
     */
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    /**
     * @brief 计算曲线面积
     */
    Float Area() const;
    /**
     * @brief 在曲线表面采样一个点
     */
    Interaction Sample(const Point2f &u, Float *pdf) const;

  private:
    // Curve Private Methods / 曲线私有方法
    /**
     * @brief 递归光线-曲线求交(逐层细分贝塞尔曲线段)
     */
    bool recursiveIntersect(const Ray &r, Float *tHit,
                            SurfaceInteraction *isect, const Point3f cp[4],
                            const Transform &rayToObject, Float u0, Float u1,
                            int depth) const;

    // Curve Private Data / 曲线私有数据
    const std::shared_ptr<CurveCommon> common; /**< 曲线共享数据 */
    const Float uMin, uMax;                    /**< 曲线段的参数范围 */
};

/**
 * @brief 创建曲线形状的工厂函数
 */
std::vector<std::shared_ptr<Shape>> CreateCurveShape(const Transform *o2w,
                                                     const Transform *w2o,
                                                     bool reverseOrientation,
                                                     const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_CURVE_H
