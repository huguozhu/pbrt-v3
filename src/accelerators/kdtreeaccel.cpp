
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


// accelerators/kdtreeaccel.cpp*
// Kd树加速结构实现：用于加速光线追踪中的求交测试，
// 采用SAH(Surface Area Heuristic)启发式策略选择最优分割平面
#include "accelerators/kdtreeaccel.h"
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>

namespace pbrt {

// KdTreeAccel Local Declarations
// KdAccelNode: Kd树节点，使用联合体存储内部节点和叶子节点的不同数据
struct KdAccelNode {
    // KdAccelNode Methods
    // InitLeaf: 初始化叶子节点，存储图元编号或指向图元索引数组的偏移
    void InitLeaf(int *primNums, int np, std::vector<int> *primitiveIndices);
    // InitInterior: 初始化内部节点，设置分割轴、子节点偏移和分割位置
    void InitInterior(int axis, int ac, Float s) {
        split = s;
        flags = axis;
        aboveChild |= (ac << 2);
    }
    // SplitPos: 获取分割平面的位置
    Float SplitPos() const { return split; }
    // nPrimitives: 获取叶子节点中的图元数量
    int nPrimitives() const { return nPrims >> 2; }
    // SplitAxis: 获取内部节点的分割轴（0=x, 1=y, 2=z）
    int SplitAxis() const { return flags & 3; }
    // IsLeaf: 判断是否为叶子节点（flags低2位=3表示叶子）
    bool IsLeaf() const { return (flags & 3) == 3; }
    // AboveChild: 获取内部节点的"上方"子节点索引
    int AboveChild() const { return aboveChild >> 2; }
    union {
        Float split;                 // Interior: 分割位置
        int onePrimitive;            // Leaf: 只有一个图元时直接存储图元编号
        int primitiveIndicesOffset;  // Leaf: 多个图元时存储在全局索引数组中的偏移
    };

  private:
    union {
        int flags;       // Both: 低2位为标志位，高30位为数据
        int nPrims;      // Leaf: 图元数量（编码在高30位）
        int aboveChild;  // Interior: 另一个子节点的索引
    };
};

// EdgeType: 边的类型（图元包围盒的起始边或结束边）
enum class EdgeType { Start, End };
// BoundEdge: 图元包围盒在分割轴上的投影边，用于SAH代价计算
struct BoundEdge {
    // BoundEdge Public Methods
    BoundEdge() {}
    BoundEdge(Float t, int primNum, bool starting) : t(t), primNum(primNum) {
        type = starting ? EdgeType::Start : EdgeType::End;
    }
    Float t;         // 分割位置
    int primNum;     // 对应的图元编号
    EdgeType type;   // 边类型：起始边或结束边
};

// KdTreeAccel Method Definitions
// KdTreeAccel构造函数：构建Kd树，使用SAH启发式选择分割平面
KdTreeAccel::KdTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
                         int isectCost, int traversalCost, Float emptyBonus,
                         int maxPrims, int maxDepth)
    : isectCost(isectCost),
      traversalCost(traversalCost),
      maxPrims(maxPrims),
      emptyBonus(emptyBonus),
      primitives(std::move(p)) {
    // Build kd-tree for accelerator
    ProfilePhase _(Prof::AccelConstruction);
    nextFreeNode = nAllocedNodes = 0;
    if (maxDepth <= 0)
        maxDepth = std::round(8 + 1.3f * Log2Int(int64_t(primitives.size())));

    // Compute bounds for kd-tree construction
    std::vector<Bounds3f> primBounds;
    primBounds.reserve(primitives.size());
    for (const std::shared_ptr<Primitive> &prim : primitives) {
        Bounds3f b = prim->WorldBound();
        bounds = Union(bounds, b);
        primBounds.push_back(b);
    }

    // Allocate working memory for kd-tree construction
    std::unique_ptr<BoundEdge[]> edges[3];
    for (int i = 0; i < 3; ++i)
        edges[i].reset(new BoundEdge[2 * primitives.size()]);
    std::unique_ptr<int[]> prims0(new int[primitives.size()]);
    std::unique_ptr<int[]> prims1(new int[(maxDepth + 1) * primitives.size()]);

    // Initialize _primNums_ for kd-tree construction
    std::unique_ptr<int[]> primNums(new int[primitives.size()]);
    for (size_t i = 0; i < primitives.size(); ++i) primNums[i] = i;

    // Start recursive construction of kd-tree
    buildTree(0, bounds, primBounds, primNums.get(), primitives.size(),
              maxDepth, edges, prims0.get(), prims1.get());
}

// InitLeaf: 初始化叶子节点，根据图元数量选择不同的存储方式
void KdAccelNode::InitLeaf(int *primNums, int np,
                           std::vector<int> *primitiveIndices) {
    flags = 3;  // 低2位=3表示叶子节点
    nPrims |= (np << 2);  // 将图元数量编码到高30位
    // Store primitive ids for leaf node
    if (np == 0)
        onePrimitive = 0;
    else if (np == 1)
        onePrimitive = primNums[0];  // 只有一个图元时直接存储编号
    else {
        // 多个图元时，在primitiveIndices数组中分配空间存储所有编号
        primitiveIndicesOffset = primitiveIndices->size();
        for (int i = 0; i < np; ++i) primitiveIndices->push_back(primNums[i]);
    }
}

KdTreeAccel::~KdTreeAccel() { FreeAligned(nodes); }

// buildTree: 递归构建Kd树，使用SAH启发式选择最优分割平面
void KdTreeAccel::buildTree(int nodeNum, const Bounds3f &nodeBounds,
                            const std::vector<Bounds3f> &allPrimBounds,
                            int *primNums, int nPrimitives, int depth,
                            const std::unique_ptr<BoundEdge[]> edges[3],
                            int *prims0, int *prims1, int badRefines) {
    CHECK_EQ(nodeNum, nextFreeNode);
    // Get next free node from _nodes_ array
    // 如果节点数组空间不足，则分配更大的数组并复制原有数据
    if (nextFreeNode == nAllocedNodes) {
        int nNewAllocNodes = std::max(2 * nAllocedNodes, 512);
        KdAccelNode *n = AllocAligned<KdAccelNode>(nNewAllocNodes);
        if (nAllocedNodes > 0) {
            memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
            FreeAligned(nodes);
        }
        nodes = n;
        nAllocedNodes = nNewAllocNodes;
    }
    ++nextFreeNode;

    // Initialize leaf node if termination criteria met
    // 递归终止条件：图元数达到下限或深度达到上限
    if (nPrimitives <= maxPrims || depth == 0) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
    // 使用SAH评估所有可能的切割位置，选择代价最小的
    int bestAxis = -1, bestOffset = -1;
    Float bestCost = Infinity;
    Float oldCost = isectCost * Float(nPrimitives);  // 不分割时的代价
    Float totalSA = nodeBounds.SurfaceArea();
    Float invTotalSA = 1 / totalSA;
    Vector3f d = nodeBounds.pMax - nodeBounds.pMin;

    // Choose which axis to split along
    // 首先尝试包围盒最长的轴
    int axis = nodeBounds.MaximumExtent();
    int retries = 0;
retrySplit:

    // Initialize edges for _axis_
    // 计算所有图元包围盒在当前轴上的投影边
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const Bounds3f &bounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(bounds.pMin[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(bounds.pMax[axis], pn, false);
    }

    // Sort _edges_ for _axis_
    // 按边的位置排序，起始边排在结束边之前
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int)e0.type < (int)e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
    // 扫描所有分割位置，计算每个分割的SAH代价
    int nBelow = 0, nAbove = nPrimitives;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End) --nAbove;
        Float edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.pMin[axis] && edgeT < nodeBounds.pMax[axis]) {
            // Compute cost for split at _i_th edge

            // Compute child surface areas for split at _edgeT_
            // 计算分割后两个子空间的表面积
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
            Float belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                 (edgeT - nodeBounds.pMin[axis]) *
                                     (d[otherAxis0] + d[otherAxis1]));
            Float aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                 (nodeBounds.pMax[axis] - edgeT) *
                                     (d[otherAxis0] + d[otherAxis1]));
            Float pBelow = belowSA * invTotalSA;
            Float pAbove = aboveSA * invTotalSA;
            Float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;
            Float cost =
                traversalCost +
                isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

            // Update best split if this is lowest cost so far
            if (cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestOffset = i;
            }
        }
        if (edges[axis][i].type == EdgeType::Start) ++nBelow;
    }
    CHECK(nBelow == nPrimitives && nAbove == 0);

    // Create leaf if no good splits were found
    // 如果当前轴没有找到好的分割，尝试其他轴
    if (bestAxis == -1 && retries < 2) {
        ++retries;
        axis = (axis + 1) % 3;
        goto retrySplit;
    }
    if (bestCost > oldCost) ++badRefines;
    // 如果分割代价过高或已经过多欠佳分割，创建叶子节点
    if ((bestCost > 4 * oldCost && nPrimitives < 16) || bestAxis == -1 ||
        badRefines == 3) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Classify primitives with respect to split
    // 根据最优分割平面分类图元：左侧和右侧
    int n0 = 0, n1 = 0;
    for (int i = 0; i < bestOffset; ++i)
        if (edges[bestAxis][i].type == EdgeType::Start)
            prims0[n0++] = edges[bestAxis][i].primNum;
    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)
        if (edges[bestAxis][i].type == EdgeType::End)
            prims1[n1++] = edges[bestAxis][i].primNum;

    // Recursively initialize children nodes
    // 递归构建左右子树
    Float tSplit = edges[bestAxis][bestOffset].t;
    Bounds3f bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tSplit;
    buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
              prims0, prims1 + nPrimitives, badRefines);
    int aboveChild = nextFreeNode;
    nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
    buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
              prims0, prims1 + nPrimitives, badRefines);
}

// Intersect: 求光线与场景中图元的最近交点，通过遍历Kd树加速相交测试
bool KdTreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    // Compute initial parametric range of ray inside kd-tree extent
    // 计算光线在场景包围盒内的参数范围
    Float tMin, tMax;
    if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }

    // Prepare to traverse kd-tree for ray
    // 准备遍历：计算方向倒数，初始化待处理栈
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    PBRT_CONSTEXPR int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    // 按顺序遍历Kd树节点
    bool hit = false;
    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {
        // Bail out if we found a hit closer than the current node
        if (ray.tMax < tMin) break;
        if (!node->IsLeaf()) {
            // Process kd-tree interior node
            // 处理内部节点：计算光线到分割平面的距离，确定遍历顺序

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            Float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            // 根据光线方向确定先遍历哪个子节点
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst =
                (ray.o[axis] < node->SplitPos()) ||
                (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;           // 近子节点（在内存中相邻）
                secondChild = &nodes[node->AboveChild()];  // 远子节点
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            // 根据tPlane决定遍历路径：将远子节点压入栈，先遍历近子节点
            if (tPlane > tMax || tPlane <= 0)
                node = firstChild;       // 分割平面在光线范围外，只需遍历近子节点
            else if (tPlane < tMin)
                node = secondChild;      // 光线从分割平面后开始，直接遍历远子节点
            else {
                // Enqueue _secondChild_ in todo list
                // 将远子节点入栈，先遍历近子节点
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;
            }
        } else {
            // Check for intersections inside leaf node
            // 处理叶子节点：测试光线与所有图元的相交
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                    primitives[node->onePrimitive];
                // Check one primitive inside leaf node
                if (p->Intersect(ray, isect)) hit = true;
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->Intersect(ray, isect)) hit = true;
                }
            }

            // Grab next node to process from todo list
            // 从栈中取出下一个待处理节点
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
        }
    }
    return hit;
}

// IntersectP: 判断光线是否与场景相交（阴影光线测试用），找到任意交点即返回
bool KdTreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    // Compute initial parametric range of ray inside kd-tree extent
    Float tMin, tMax;
    if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    PBRT_CONSTEXPR int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;
    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {
        if (node->IsLeaf()) {
            // Check for shadow ray intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                    primitives[node->onePrimitive];
                if (p->IntersectP(ray)) {
                    return true;  // 找到任意交点即可返回
                }
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int primitiveIndex =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &prim =
                        primitives[primitiveIndex];
                    if (prim->IntersectP(ray)) {
                        return true;
                    }
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
        } else {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            Float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst =
                (ray.o[axis] < node->SplitPos()) ||
                (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tPlane > tMax || tPlane <= 0)
                node = firstChild;
            else if (tPlane < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in todo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;
            }
        }
    }
    return false;
}

// CreateKdTreeAccelerator: 工厂函数，从参数集中读取配置并创建KdTreeAccel加速器
std::shared_ptr<KdTreeAccel> CreateKdTreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    int isectCost = ps.FindOneInt("intersectcost", 80);
    int travCost = ps.FindOneInt("traversalcost", 1);
    Float emptyBonus = ps.FindOneFloat("emptybonus", 0.5f);
    int maxPrims = ps.FindOneInt("maxprims", 1);
    int maxDepth = ps.FindOneInt("maxdepth", -1);
    return std::make_shared<KdTreeAccel>(std::move(prims), isectCost, travCost, emptyBonus,
                                         maxPrims, maxDepth);
}

}  // namespace pbrt
