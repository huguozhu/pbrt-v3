
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

#ifndef PBRT_SHAPES_CONE_H
#define PBRT_SHAPES_CONE_H

// shapes/cone.h*
/**
 * @file cone.h
 * @brief 圆锥体(Cone)几何体模块
 *
 * 该模块定义了圆锥体几何体，用于表示沿Z轴延伸的圆锥形状。
 * 圆锥体通过底面半径(radius)、高度(height)和最大角度(phiMax)来定义，
 * 支持裁剪( clipping)参数来创建部分圆锥。
 */
#include "shape.h"

namespace pbrt {

// Cone Declarations / 圆锥体声明
class Cone : public Shape {
  public:
    // Cone Public Methods / 圆锥体公有方法
    /**
     * @brief 构造函数：创建圆锥体
     * @param o2w 对象到世界变换
     * @param w2o 世界到对象变换
     * @param reverseOrientation 是否反转法线方向
     * @param height 圆锥高度
     * @param radius 圆锥底面半径
     * @param phiMax 最大角度范围(度)，用于创建部分圆锥
     */
    Cone(const Transform *o2w, const Transform *w2o, bool reverseOrientation,
         Float height, Float radius, Float phiMax);
    /**
     * @brief 计算对象空间包围盒
     * @return 对象空间中的包围盒
     */
    Bounds3f ObjectBound() const;
    /**
     * @brief 光线-圆锥体求交测试(完整求交)
     * @param ray 光线
     * @param tHit 返回交点距离参数
     * @param isect 返回交点表面交互信息
     * @param testAlphaTexture 是否测试Alpha纹理
     * @return 是否相交
     */
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture) const;
    /**
     * @brief 光线-圆锥体求交测试(仅判断是否相交)
     * @param ray 光线
     * @param testAlphaTexture 是否测试Alpha纹理
     * @return 是否相交
     */
    bool IntersectP(const Ray &ray, bool testAlphaTexture) const;
    /**
     * @brief 计算圆锥体表面积
     * @return 表面积
     */
    Float Area() const;
    /**
     * @brief 在圆锥体表面采样一个点
     * @param u 均匀随机数
     * @param pdf 返回概率密度
     * @return 采样得到的交互点
     */
    Interaction Sample(const Point2f &u, Float *pdf) const;

  protected:
    // Cone Private Data / 圆锥体私有数据
    const Float radius;   /**< 圆锥底面半径 */
    const Float height;   /**< 圆锥高度 */
    const Float phiMax;   /**< 最大角度范围(弧度) */
};

/**
 * @brief 创建圆锥体形状的工厂函数
 * @param o2w 对象到世界变换
 * @param w2o 世界到对象变换
 * @param reverseOrientation 是否反转法线方向
 * @param params 参数集
 * @return 创建的圆锥体共享指针
 */
std::shared_ptr<Cone> CreateConeShape(const Transform *o2w,
                                      const Transform *w2o,
                                      bool reverseOrientation,
                                      const ParamSet &params);

}  // namespace pbrt

#endif  // PBRT_SHAPES_CONE_H
