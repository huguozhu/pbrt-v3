
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

#ifndef PBRT_SHAPES_TRIANGLE_H
#define PBRT_SHAPES_TRIANGLE_H

// shapes/triangle.h*
/**
 * @file triangle.h
 * @brief 三角形网格(Triangle)几何体模块
 *
 * 定义了三角形网格几何体，是pbrt中最常用的几何表示方式。
 * 包括TriangleMesh(网格数据容器)和Triangle(单个三角形)。
 * 支持法线、切线、UV坐标、Alpha纹理和阴影Alpha纹理。
 */
#include "shape.h"
#include "stats.h"
#include <map>

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Triangle meshes", triMeshBytes);

// Triangle Declarations / 三角形声明

/**
 * @brief 三角形网格数据结构
 * 存储顶点位置、法线、切线、UV坐标以及面索引等网格数据
 */
struct TriangleMesh {
    // TriangleMesh Public Methods / 网格公有方法
    TriangleMesh(const Transform &ObjectToWorld, int nTriangles,
                 const int *vertexIndices, int nVertices, const Point3f *P,
                 const Vector3f *S, const Normal3f *N, const Point2f *uv,
                 const std::shared_ptr<Texture<Float>> &alphaMask,
                 const std::shared_ptr<Texture<Float>> &shadowAlphaMask,
                 const int *faceIndices);

    // TriangleMesh Data / 网格数据
    const int nTriangles, nVertices; /**< 三角形数量和顶点数量 */
    std::vector<int> vertexIndices;  /**< 顶点索引数组 */
    std::unique_ptr<Point3f[]> p;    /**< 顶点位置(世界空间) */
    std::unique_ptr<Normal3f[]> n;   /**< 顶点法线(可选) */
    std::unique_ptr<Vector3f[]> s;   /**< 顶点切线(可选) */
    std::unique_ptr<Point2f[]> uv;   /**< 顶点UV坐标(可选) */
    std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask; /**< 透明度纹理 */
    std::vector<int> faceIndices;    /**< 面片索引(可选) */
};

/**
 * @brief 单个三角形
 * 引用TriangleMesh中的数据，实现Shape接口。
 */
class Triangle : public Shape {
  public:
    // Triangle Public Methods / 三角形公有方法
    Triangle(const Transform *ObjectToWorld, const Transform *WorldToObject,
             bool reverseOrientation, const std::shared_ptr<TriangleMesh> &mesh,
             int triNumber)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation), mesh(mesh) {
        v = &mesh->vertexIndices[3 * triNumber];
        triMeshBytes += sizeof(*this);
        faceIndex = mesh->faceIndices.size() ? mesh->faceIndices[triNumber] : 0;
    }
    Bounds3f ObjectBound() const;
    Bounds3f WorldBound() const;
    bool Intersect(const Ray &ray, Float *tHit, SurfaceInteraction *isect,
                   bool testAlphaTexture = true) const;
    bool IntersectP(const Ray &ray, bool testAlphaTexture = true) const;
    Float Area() const;

    using Shape::Sample;  // Bring in the other Sample() overload. / 引入其它Sample重载版本
    Interaction Sample(const Point2f &u, Float *pdf) const;

    /**
     * @brief 计算三角形相对于某点的立体角
     * 使用Girard球面三角形面积定理计算
     */
    Float SolidAngle(const Point3f &p, int nSamples = 0) const;

  private:
    // Triangle Private Methods / 三角形私有方法
    /**
     * @brief 获取三角形的UV坐标(无UV时使用默认值)
     */
    void GetUVs(Point2f uv[3]) const {
        if (mesh->uv) {
            uv[0] = mesh->uv[v[0]];
            uv[1] = mesh->uv[v[1]];
            uv[2] = mesh->uv[v[2]];
        } else {
            uv[0] = Point2f(0, 0);
            uv[1] = Point2f(1, 0);
            uv[2] = Point2f(1, 1);
        }
    }

    // Triangle Private Data / 三角形私有数据
    std::shared_ptr<TriangleMesh> mesh; /**< 所属三角形网格 */
    const int *v;                       /**< 顶点索引指针 */
    int faceIndex;                      /**< 面片索引 */
};

std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    int nTriangles, const int *vertexIndices, int nVertices, const Point3f *p,
    const Vector3f *s, const Normal3f *n, const Point2f *uv,
    const std::shared_ptr<Texture<Float>> &alphaTexture,
    const std::shared_ptr<Texture<Float>> &shadowAlphaTexture,
    const int *faceIndices = nullptr);
std::vector<std::shared_ptr<Shape>> CreateTriangleMeshShape(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    const ParamSet &params,
    std::map<std::string, std::shared_ptr<Texture<Float>>> *floatTextures =
        nullptr);

bool WritePlyFile(const std::string &filename, int nTriangles,
                  const int *vertexIndices, int nVertices, const Point3f *P,
                  const Vector3f *S, const Normal3f *N, const Point2f *UV,
                  const int *faceIndices);

}  // namespace pbrt

#endif  // PBRT_SHAPES_TRIANGLE_H
