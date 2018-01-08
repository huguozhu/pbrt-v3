
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

#ifndef PBRT_CORE_PRIMITIVE_H
#define PBRT_CORE_PRIMITIVE_H

// core/primitive.h*
#include "pbrt.h"
#include "shape.h"
#include "material.h"
#include "medium.h"

namespace pbrt {

// Primitive Declarations
// Shape类：用于表示物体的形状、mesh等信息
// Primitive类：包括了物体的mesh、材质、光线等信息（自身是发光体）
class Primitive {
  public:
    // Primitive Interface
    virtual ~Primitive();
    virtual Bounds3f WorldBound() const = 0;								// 在世界坐标系的包围盒
    virtual bool Intersect(const Ray &r, SurfaceInteraction *) const = 0;	// 与射线ray是否相交，且保存相交信息
    virtual bool IntersectP(const Ray &r) const = 0;						// 只简单判断是否与ray相交
    virtual const AreaLight *GetAreaLight() const = 0;	// 如果图元本身表示一个光源，则返回非null，AreaLight: 描述图元的光源发射特征
    virtual const Material *GetMaterial() const = 0;
    virtual void ComputeScatteringFunctions(SurfaceInteraction *isect,		// 计算散射的函数
                                            MemoryArena &arena,
                                            TransportMode mode,
                                            bool allowMultipleLobes) const = 0;
};

// GeometricPrimitive Declarations
// 几何图元：表示场景中的单一形状，比如球体或者平面等
class GeometricPrimitive : public Primitive {
  public:
    // GeometricPrimitive Public Methods
    virtual Bounds3f WorldBound() const;
    virtual bool Intersect(const Ray &r, SurfaceInteraction *isect) const;
    virtual bool IntersectP(const Ray &r) const;
    GeometricPrimitive(const std::shared_ptr<Shape> &shape,
                       const std::shared_ptr<Material> &material,
                       const std::shared_ptr<AreaLight> &areaLight,
                       const MediumInterface &mediumInterface);
    const AreaLight *GetAreaLight() const;
    const Material *GetMaterial() const;
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // GeometricPrimitive Private Data
    std::shared_ptr<Shape> shape;
    std::shared_ptr<Material> material;
    std::shared_ptr<AreaLight> areaLight;
    MediumInterface mediumInterface;
};

// TransformedPrimitive Declarations
class TransformedPrimitive : public Primitive {
  public:
    // TransformedPrimitive Public Methods
    TransformedPrimitive(std::shared_ptr<Primitive> &primitive,
                         const AnimatedTransform &PrimitiveToWorld);
    bool Intersect(const Ray &r, SurfaceInteraction *in) const;
    bool IntersectP(const Ray &r) const;
    const AreaLight *GetAreaLight() const { return nullptr; }
    const Material *GetMaterial() const { return nullptr; }
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const {
        LOG(FATAL) <<
            "TransformedPrimitive::ComputeScatteringFunctions() shouldn't be "
            "called";
    }
    Bounds3f WorldBound() const {
        return PrimitiveToWorld.MotionBounds(primitive->WorldBound());
    }

  private:
    // TransformedPrimitive Private Data
    std::shared_ptr<Primitive> primitive;
    const AnimatedTransform PrimitiveToWorld;
};

// Aggregate Declarations
// 一个把多个物体归为一组的接口,有两种方案：空间划分和对象划分
// 空间划分：将3D空间分解为多个区域，算法：GridAccel、KdTreeAccel
// 对象划分：将场景分解为较小的对象集，例如房间包括墙面、底板和桌椅等，算法：BVHAccel
class Aggregate : public Primitive {
  public:
    // Aggregate Public Methods
    const AreaLight *GetAreaLight() const;
    const Material *GetMaterial() const;
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const;
};

}  // namespace pbrt

#endif  // PBRT_CORE_PRIMITIVE_H
