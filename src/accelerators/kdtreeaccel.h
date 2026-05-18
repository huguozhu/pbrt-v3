
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

#ifndef PBRT_ACCELERATORS_KDTREEACCEL_H
#define PBRT_ACCELERATORS_KDTREEACCEL_H

// accelerators/kdtreeaccel.h*
// KdTreeAccel: Kd树加速结构，用于加速光线与场景图元的相交测试。
// 通过递归地将空间沿坐标轴分割，构建平衡的Kd树以减少不必要的相交测试
#include "pbrt.h"
#include "primitive.h"

namespace pbrt {

// KdTreeAccel Declarations
struct KdAccelNode;   // Kd树节点（内部节点或叶子节点）
struct BoundEdge;     // 图元包围盒的边，用于SAH代价计算
class KdTreeAccel : public Aggregate {
  public:
    // KdTreeAccel Public Methods
    // 构造函数：构建Kd树加速结构，参数包括相交测试代价、遍历代价、空区域奖励等
    KdTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
                int isectCost = 80, int traversalCost = 1,
                Float emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
    // WorldBound: 返回整个加速结构的包围盒
    Bounds3f WorldBound() const { return bounds; }
    ~KdTreeAccel();
    // Intersect: 求光线与场景的最近交点
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
    // IntersectP: 判断光线是否与场景相交（阴影检测用）
    bool IntersectP(const Ray &ray) const;

  private:
    // KdTreeAccel Private Methods
    // buildTree: 递归构建Kd树
    void buildTree(int nodeNum, const Bounds3f &bounds,
                   const std::vector<Bounds3f> &primBounds, int *primNums,
                   int nprims, int depth,
                   const std::unique_ptr<BoundEdge[]> edges[3], int *prims0,
                   int *prims1, int badRefines = 0);

    // KdTreeAccel Private Data
    const int isectCost, traversalCost, maxPrims;  // SAH代价参数和每个叶子节点最大图元数
    const Float emptyBonus;                         // 空区域奖励系数
    std::vector<std::shared_ptr<Primitive>> primitives;  // 所有图元
    std::vector<int> primitiveIndices;                   // 叶子节点中图元的索引
    KdAccelNode *nodes;                                  // Kd树节点数组
    int nAllocedNodes, nextFreeNode;                     // 已分配节点数和下一个空闲节点
    Bounds3f bounds;                                     // 场景总包围盒
};

// KdToDo: 遍历Kd树时的待处理节点栈元素
struct KdToDo {
    const KdAccelNode *node;  // 待处理节点
    Float tMin, tMax;          // 光线在该节点空间范围内的参数区间
};

// CreateKdTreeAccelerator: 根据参数集创建Kd树加速器工厂函数
std::shared_ptr<KdTreeAccel> CreateKdTreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps);

}  // namespace pbrt

#endif  // PBRT_ACCELERATORS_KDTREEACCEL_H
