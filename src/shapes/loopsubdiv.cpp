
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


// shapes/loopsubdiv.cpp*
/**
 * @file loopsubdiv.cpp
 * @brief Loop细分曲面算法的实现
 *
 * 实现了Catmull-Clark风格的Loop细分算法，将任意三角网格通过多次细分
 * 生成光滑曲面。包含顶点权重计算、一邻域顶点收集、极限曲面推点等功能。
 */
#include "shapes/loopsubdiv.h"
#include "shapes/triangle.h"
#include "paramset.h"
#include <set>
#include <map>

namespace pbrt {

struct SDFace;
struct SDVertex;

// LoopSubdiv Macros / Loop细分宏定义
#define NEXT(i) (((i) + 1) % 3)
#define PREV(i) (((i) + 2) % 3)

// LoopSubdiv Local Structures / Loop细分局部数据结构
/**
 * @brief 细分曲面顶点结构
 * 存储顶点位置、所属面、子顶点以及规则/边界标志
 */
struct SDVertex {
    // SDVertex Constructor / 顶点构造函数
    SDVertex(const Point3f &p = Point3f(0, 0, 0)) : p(p) {}

    // SDVertex Methods / 顶点方法
    int valence();              /**< 计算顶点的价(valence，即邻接边数) */
    void oneRing(Point3f *p);   /**< 收集一邻域顶点 */
    Point3f p;                  /**< 顶点位置 */
    SDFace *startFace = nullptr;   /**< 起始面指针 */
    SDVertex *child = nullptr;     /**< 细分后的子顶点 */
    bool regular = false, boundary = false; /**< 是否正则顶点/是否为边界顶点 */
};

/**
 * @brief 细分曲面面片结构
 * 存储三角形面的三个顶点、相邻面以及细分后的四个子面
 */
struct SDFace {
    // SDFace Constructor / 面片构造函数
    SDFace() {
        for (int i = 0; i < 3; ++i) {
            v[i] = nullptr;
            f[i] = nullptr;
        }
        for (int i = 0; i < 4; ++i) children[i] = nullptr;
    }

    // SDFace Methods / 面片方法
    int vnum(SDVertex *vert) const {  /**< 返回顶点在面中的索引(0,1,2) */
        for (int i = 0; i < 3; ++i)
            if (v[i] == vert) return i;
        LOG(FATAL) << "Basic logic error in SDFace::vnum()";
        return -1;
    }
    SDFace *nextFace(SDVertex *vert) { return f[vnum(vert)]; }  /**< 获取顶点的下一个邻面 */
    SDFace *prevFace(SDVertex *vert) { return f[PREV(vnum(vert))]; } /**< 获取顶点的上一个邻面 */
    SDVertex *nextVert(SDVertex *vert) { return v[NEXT(vnum(vert))]; } /**< 获取顺时针方向的下一个顶点 */
    SDVertex *prevVert(SDVertex *vert) { return v[PREV(vnum(vert))]; } /**< 获取逆时针方向的上一个顶点 */
    SDVertex *otherVert(SDVertex *v0, SDVertex *v1) { /**< 获取边(v0,v1)对面的第三个顶点 */
        for (int i = 0; i < 3; ++i)
            if (v[i] != v0 && v[i] != v1) return v[i];
        LOG(FATAL) << "Basic logic error in SDVertex::otherVert()";
        return nullptr;
    }
    SDVertex *v[3];      /**< 三个顶点 */
    SDFace *f[3];        /**< 三个相邻面 */
    SDFace *children[4]; /**< 细分后的四个子面(三个角面+一个中心面) */
};

/**
 * @brief 细分曲面边结构
 * 用于构建半边数据结构，存储边的两个端点和邻接面
 */
struct SDEdge {
    // SDEdge Constructor / 边构造函数
    SDEdge(SDVertex *v0 = nullptr, SDVertex *v1 = nullptr) {
        v[0] = std::min(v0, v1);  /**< 保证端点有序 */
        v[1] = std::max(v0, v1);
        f[0] = f[1] = nullptr;
        f0edgeNum = -1;
    }

    // SDEdge Comparison Function / 边比较函数(用于std::set/map排序)
    bool operator<(const SDEdge &e2) const {
        if (v[0] == e2.v[0]) return v[1] < e2.v[1];
        return v[0] < e2.v[0];
    }
    SDVertex *v[2];   /**< 边的两个端点 */
    SDFace *f[2];     /**< 边两侧的邻接面 */
    int f0edgeNum;    /**< 在第一个面中该边的索引 */
};

// LoopSubdiv Local Declarations / Loop细分局部函数声明
static Point3f weightOneRing(SDVertex *vert, Float beta);
static Point3f weightBoundary(SDVertex *vert, Float beta);

// LoopSubdiv Inline Functions / Loop细分内联函数
/**
 * @brief 计算顶点的价(valence)
 * 内部顶点: 邻接面数量即为价
 * 边界顶点: 邻接面数量+1(包含边界本身)
 */
inline int SDVertex::valence() {
    SDFace *f = startFace;
    if (!boundary) {
        // Compute valence of interior vertex
        int nf = 1;
        while ((f = f->nextFace(this)) != startFace) ++nf;
        return nf;
    } else {
        // Compute valence of boundary vertex
        int nf = 1;
        while ((f = f->nextFace(this)) != nullptr) ++nf;
        f = startFace;
        while ((f = f->prevFace(this)) != nullptr) ++nf;
        return nf + 1;
    }
}

/**
 * @brief 计算内部顶点的细分权重beta
 * Loop细分中偶数顶点的更新权重
 */
inline Float beta(int valence) {
    if (valence == 3)
        return 3.f / 16.f;
    else
        return 3.f / (8.f * valence);
}

/**
 * @brief 计算极限曲面推点权重gamma
 * 将细分曲面顶点推向极限位置时使用的权重
 */
inline Float loopGamma(int valence) {
    return 1.f / (valence + 3.f / (8.f * beta(valence)));
}

// LoopSubdiv Function Definitions / Loop细分函数实现
/**
 * @brief Loop细分主函数
 *
 * 执行指定层数的Loop细分，将粗网格细化为光滑曲面。
 * 主要步骤：
 * 1. 分配顶点和面片
 * 2. 设置邻接关系
 * 3. 执行多级细分(更新偶数顶点、创建奇数顶点、更新拓扑)
 * 4. 推点到极限曲面
 * 5. 输出三角形网格
 */
static std::vector<std::shared_ptr<Shape>> LoopSubdivide(
    const Transform *ObjectToWorld, const Transform *WorldToObject,
    bool reverseOrientation, int nLevels, int nIndices,
    const int *vertexIndices, int nVertices, const Point3f *p) {
    std::vector<SDVertex *> vertices;
    std::vector<SDFace *> faces;
    // Allocate _LoopSubdiv_ vertices and faces / 分配Loop细分所需的顶点和面片
    std::unique_ptr<SDVertex[]> verts(new SDVertex[nVertices]);
    for (int i = 0; i < nVertices; ++i) {
        verts[i] = SDVertex(p[i]);
        vertices.push_back(&verts[i]);
    }
    int nFaces = nIndices / 3;
    std::unique_ptr<SDFace[]> fs(new SDFace[nFaces]);
    for (int i = 0; i < nFaces; ++i) faces.push_back(&fs[i]);

    // Set face to vertex pointers / 设置面到顶点的指针
    const int *vp = vertexIndices;
    for (int i = 0; i < nFaces; ++i, vp += 3) {
        SDFace *f = faces[i];
        for (int j = 0; j < 3; ++j) {
            SDVertex *v = vertices[vp[j]];
            f->v[j] = v;
            v->startFace = f;
        }
    }

    // Set neighbor pointers in _faces_ / 设置面的邻接指针(通过边集建立邻接关系)
    std::set<SDEdge> edges;
    for (int i = 0; i < nFaces; ++i) {
        SDFace *f = faces[i];
        for (int edgeNum = 0; edgeNum < 3; ++edgeNum) {
            // Update neighbor pointer for _edgeNum_
            int v0 = edgeNum, v1 = NEXT(edgeNum);
            SDEdge e(f->v[v0], f->v[v1]);
            if (edges.find(e) == edges.end()) {
                // Handle new edge
                e.f[0] = f;
                e.f0edgeNum = edgeNum;
                edges.insert(e);
            } else {
                // Handle previously seen edge
                e = *edges.find(e);
                e.f[0]->f[e.f0edgeNum] = f;
                f->f[edgeNum] = e.f[0];
                edges.erase(e);
            }
        }
    }

    // Finish vertex initialization / 完成顶点初始化(判断边界和正则性)
    for (int i = 0; i < nVertices; ++i) {
        SDVertex *v = vertices[i];
        SDFace *f = v->startFace;
        do {
            f = f->nextFace(v);
        } while (f && f != v->startFace);
        v->boundary = (f == nullptr);
        if (!v->boundary && v->valence() == 6)
            v->regular = true;
        else if (v->boundary && v->valence() == 4)
            v->regular = true;
        else
            v->regular = false;
    }

    // Refine _LoopSubdiv_ into triangles / 循环执行细分细化
    std::vector<SDFace *> f = faces;
    std::vector<SDVertex *> v = vertices;
    MemoryArena arena;
    for (int i = 0; i < nLevels; ++i) {
        // Update _f_ and _v_ for next level of subdivision / 为下一级细分更新面和顶点集合
        std::vector<SDFace *> newFaces;
        std::vector<SDVertex *> newVertices;

        // Allocate next level of children in mesh tree / 在网格树中分配下一级子节点
        for (SDVertex *vertex : v) {
            vertex->child = arena.Alloc<SDVertex>();
            vertex->child->regular = vertex->regular;
            vertex->child->boundary = vertex->boundary;
            newVertices.push_back(vertex->child);
        }
        for (SDFace *face : f) {
            for (int k = 0; k < 4; ++k) {
                face->children[k] = arena.Alloc<SDFace>();
                newFaces.push_back(face->children[k]);
            }
        }

        // Update vertex positions and create new edge vertices / 更新顶点位置并创建新边顶点

        // Update vertex positions for even vertices / 更新偶数顶点位置
        for (SDVertex *vertex : v) {
            if (!vertex->boundary) {
                // Apply one-ring rule for even vertex / 应用一邻域规则更新内部偶数顶点
                if (vertex->regular)
                    vertex->child->p = weightOneRing(vertex, 1.f / 16.f);
                else
                    vertex->child->p =
                        weightOneRing(vertex, beta(vertex->valence()));
            } else {
                // Apply boundary rule for even vertex / 应用边界规则更新边界偶数顶点
                vertex->child->p = weightBoundary(vertex, 1.f / 8.f);
            }
        }

        // Compute new odd edge vertices / 计算新的奇数边顶点
        std::map<SDEdge, SDVertex *> edgeVerts;
        for (SDFace *face : f) {
            for (int k = 0; k < 3; ++k) {
                // Compute odd vertex on _k_th edge / 计算第k条边上的奇数顶点
                SDEdge edge(face->v[k], face->v[NEXT(k)]);
                SDVertex *vert = edgeVerts[edge];
                if (!vert) {
                    // Create and initialize new odd vertex / 创建并初始化新的奇数顶点
                    vert = arena.Alloc<SDVertex>();
                    newVertices.push_back(vert);
                    vert->regular = true;
                    vert->boundary = (face->f[k] == nullptr);
                    vert->startFace = face->children[3];

                    // Apply edge rules to compute new vertex position / 应用边规则计算新顶点位置
                    if (vert->boundary) {
                        vert->p = 0.5f * edge.v[0]->p;
                        vert->p += 0.5f * edge.v[1]->p;
                    } else {
                        vert->p = 3.f / 8.f * edge.v[0]->p;
                        vert->p += 3.f / 8.f * edge.v[1]->p;
                        vert->p += 1.f / 8.f *
                                   face->otherVert(edge.v[0], edge.v[1])->p;
                        vert->p +=
                            1.f / 8.f *
                            face->f[k]->otherVert(edge.v[0], edge.v[1])->p;
                    }
                    edgeVerts[edge] = vert;
                }
            }
        }

        // Update new mesh topology / 更新新网格拓扑

        // Update even vertex face pointers / 更新偶数顶点的面指针
        for (SDVertex *vertex : v) {
            int vertNum = vertex->startFace->vnum(vertex);
            vertex->child->startFace = vertex->startFace->children[vertNum];
        }

        // Update face neighbor pointers / 更新面片邻接指针
        for (SDFace *face : f) {
            for (int j = 0; j < 3; ++j) {
                // Update children _f_ pointers for siblings
                face->children[3]->f[j] = face->children[NEXT(j)];
                face->children[j]->f[NEXT(j)] = face->children[3];

                // Update children _f_ pointers for neighbor children
                SDFace *f2 = face->f[j];
                face->children[j]->f[j] =
                    f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
                f2 = face->f[PREV(j)];
                face->children[j]->f[PREV(j)] =
                    f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
            }
        }

        // Update face vertex pointers / 更新面片顶点指针
        for (SDFace *face : f) {
            for (int j = 0; j < 3; ++j) {
                // Update child vertex pointer to new even vertex
                face->children[j]->v[j] = face->v[j]->child;

                // Update child vertex pointer to new odd vertex
                SDVertex *vert =
                    edgeVerts[SDEdge(face->v[j], face->v[NEXT(j)])];
                face->children[j]->v[NEXT(j)] = vert;
                face->children[NEXT(j)]->v[j] = vert;
                face->children[3]->v[j] = vert;
            }
        }

        // Prepare for next level of subdivision / 准备下一级细分
        f = newFaces;
        v = newVertices;
    }

    // Push vertices to limit surface / 将顶点推向极限曲面
    std::unique_ptr<Point3f[]> pLimit(new Point3f[v.size()]);
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i]->boundary)
            pLimit[i] = weightBoundary(v[i], 1.f / 5.f);
        else
            pLimit[i] = weightOneRing(v[i], loopGamma(v[i]->valence()));
    }
    for (size_t i = 0; i < v.size(); ++i) v[i]->p = pLimit[i];

    // Compute vertex tangents on limit surface / 计算极限曲面上的顶点切向
    std::vector<Normal3f> Ns;
    Ns.reserve(v.size());
    std::vector<Point3f> pRing(16, Point3f());
    for (SDVertex *vertex : v) {
        Vector3f S(0, 0, 0), T(0, 0, 0);
        int valence = vertex->valence();
        if (valence > (int)pRing.size()) pRing.resize(valence);
        vertex->oneRing(&pRing[0]);
        if (!vertex->boundary) {
            // Compute tangents of interior face
            for (int j = 0; j < valence; ++j) {
                S += std::cos(2 * Pi * j / valence) * Vector3f(pRing[j]);
                T += std::sin(2 * Pi * j / valence) * Vector3f(pRing[j]);
            }
        } else {
            // Compute tangents of boundary face
            S = pRing[valence - 1] - pRing[0];
            if (valence == 2)
                T = Vector3f(pRing[0] + pRing[1] - 2 * vertex->p);
            else if (valence == 3)
                T = pRing[1] - vertex->p;
            else if (valence == 4)  // regular
                T = Vector3f(-1 * pRing[0] + 2 * pRing[1] + 2 * pRing[2] +
                             -1 * pRing[3] + -2 * vertex->p);
            else {
                Float theta = Pi / float(valence - 1);
                T = Vector3f(std::sin(theta) * (pRing[0] + pRing[valence - 1]));
                for (int k = 1; k < valence - 1; ++k) {
                    Float wt = (2 * std::cos(theta) - 2) * std::sin((k)*theta);
                    T += Vector3f(wt * pRing[k]);
                }
                T = -T;
            }
        }
        Ns.push_back(Normal3f(Cross(S, T)));
    }

    // Create triangle mesh from subdivision mesh / 从细分网格创建三角形网格
    {
        size_t ntris = f.size();
        std::unique_ptr<int[]> verts(new int[3 * ntris]);
        int *vp = verts.get();
        size_t totVerts = v.size();
        std::map<SDVertex *, int> usedVerts;
        for (size_t i = 0; i < totVerts; ++i) usedVerts[v[i]] = i;
        for (size_t i = 0; i < ntris; ++i) {
            for (int j = 0; j < 3; ++j) {
                *vp = usedVerts[f[i]->v[j]];
                ++vp;
            }
        }
        return CreateTriangleMesh(ObjectToWorld, WorldToObject,
                                  reverseOrientation, ntris, verts.get(),
                                  totVerts, pLimit.get(), nullptr, &Ns[0],
                                  nullptr, nullptr, nullptr);
    }
}

/**
 * @brief 创建Loop细分曲面形状的工厂函数
 * 从参数集中读取控制网格顶点、索引和细分层数
 */
std::vector<std::shared_ptr<Shape>> CreateLoopSubdiv(const Transform *o2w,
                                                     const Transform *w2o,
                                                     bool reverseOrientation,
                                                     const ParamSet &params) {
    int nLevels = params.FindOneInt("levels",
                                    params.FindOneInt("nlevels", 3));
    int nps, nIndices;
    const int *vertexIndices = params.FindInt("indices", &nIndices);
    const Point3f *P = params.FindPoint3f("P", &nps);
    if (!vertexIndices) {
        Error("Vertex indices \"indices\" not provided for LoopSubdiv shape.");
        return std::vector<std::shared_ptr<Shape>>();
    }
    if (!P) {
        Error("Vertex positions \"P\" not provided for LoopSubdiv shape.");
        return std::vector<std::shared_ptr<Shape>>();
    }

    // don't actually use this for now...
    std::string scheme = params.FindOneString("scheme", "loop");
    return LoopSubdivide(o2w, w2o, reverseOrientation, nLevels, nIndices,
                         vertexIndices, nps, P);
}

/**
 * @brief 计算内部顶点的加权一邻域位置
 * 使用beta权重对顶点及其一邻域顶点进行加权平均
 */
static Point3f weightOneRing(SDVertex *vert, Float beta) {
    // Put _vert_ one-ring in _pRing_
    int valence = vert->valence();
    Point3f *pRing = ALLOCA(Point3f, valence);
    vert->oneRing(pRing);
    Point3f p = (1 - valence * beta) * vert->p;
    for (int i = 0; i < valence; ++i) p += beta * pRing[i];
    return p;
}

/**
 * @brief 收集顶点的所有一邻域顶点
 * 对内部顶点，逆时针遍历周围面的邻接顶点
 * 对边界顶点，从边界开始沿两个方向遍历
 */
void SDVertex::oneRing(Point3f *p) {
    if (!boundary) {
        // Get one-ring vertices for interior vertex
        SDFace *face = startFace;
        do {
            *p++ = face->nextVert(this)->p;
            face = face->nextFace(this);
        } while (face != startFace);
    } else {
        // Get one-ring vertices for boundary vertex
        SDFace *face = startFace, *f2;
        while ((f2 = face->nextFace(this)) != nullptr) face = f2;
        *p++ = face->nextVert(this)->p;
        do {
            *p++ = face->prevVert(this)->p;
            face = face->prevFace(this);
        } while (face != nullptr);
    }
}

/**
 * @brief 计算边界顶点的加权位置
 * 边界顶点的更新只考虑边界上的两个邻点
 */
static Point3f weightBoundary(SDVertex *vert, Float beta) {
    // Put _vert_ one-ring in _pRing_
    int valence = vert->valence();
    Point3f *pRing = ALLOCA(Point3f, valence);
    vert->oneRing(pRing);
    Point3f p = (1 - 2 * beta) * vert->p;
    p += beta * pRing[0];
    p += beta * pRing[valence - 1];
    return p;
}

}  // namespace pbrt
