
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


// shapes/disk.cpp*
/**
 * @file disk.cpp
 * @brief 圆盘(Disk)几何体的实现
 *
 * 实现了圆盘的构造、光线求交、面积计算和采样功能。
 * 圆盘是一个平面圆环，位于z=height平面上。
 */
#include "shapes/disk.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// Disk Method Definitions / 圆盘方法实现
/**
 * @brief 计算圆盘在对象空间的包围盒
 */
Bounds3f Disk::ObjectBound() const {
    return Bounds3f(Point3f(-radius, -radius, height),
                    Point3f(radius, radius, height));
}

/**
 * @brief 光线-圆盘求交(完整求交)
 *
 * 计算光线与z=height平面的交点，然后测试该点是否在圆盘内外半径范围
 * 以及角度范围内。
 */
bool Disk::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                     bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersect);
    // Transform _Ray_ to object space / 将光线变换到对象空间
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute plane intersection for disk / 计算光线与圆盘所在平面的交点

    // Reject disk intersections for rays parallel to the disk's plane / 拒斥平行于盘面的光线
    if (ray.d.z == 0) return false;
    Float tShapeHit = (height - ray.o.z) / ray.d.z;
    if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

    // See if hit point is inside disk radii and $\phimax$ / 检查交点是否在圆盘半径范围内
    Point3f pHit = ray(tShapeHit);
    Float dist2 = pHit.x * pHit.x + pHit.y * pHit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
        return false;

    // Test disk $\phi$ value against $\phimax$ / 测试角度phi是否在phimax范围内
    Float phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;
    if (phi > phiMax) return false;

    // Find parametric representation of disk hit / 计算圆盘交点的参数化表示
    Float u = phi / phiMax;
    Float rHit = std::sqrt(dist2);
    Float v = (radius - rHit) / (radius - innerRadius);
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv =
        Vector3f(pHit.x, pHit.y, 0.) * (innerRadius - radius) / rHit;
    Normal3f dndu(0, 0, 0), dndv(0, 0, 0);

    // Refine disk intersection point / 精化圆盘交点(确保在正确平面上)
    pHit.z = height;

    // Compute error bounds for disk intersection / 计算圆盘交点误差边界
    Vector3f pError(0, 0, 0);

    // Initialize _SurfaceInteraction_ from parametric information
    *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
                                                 -ray.d, dpdu, dpdv, dndu, dndv,
                                                 ray.time, this));

    // Update _tHit_ for quadric intersection
    *tHit = (Float)tShapeHit;
    return true;
}

/**
 * @brief 光线-圆盘求交测试(仅判断是否相交)
 */
bool Disk::IntersectP(const Ray &r, bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersectP);
    // Transform _Ray_ to object space / 将光线变换到对象空间
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute plane intersection for disk / 计算光线与圆盘所在平面的交点

    // Reject disk intersections for rays parallel to the disk's plane / 拒斥平行于盘面的光线
    if (ray.d.z == 0) return false;
    Float tShapeHit = (height - ray.o.z) / ray.d.z;
    if (tShapeHit <= 0 || tShapeHit >= ray.tMax) return false;

    // See if hit point is inside disk radii and $\phimax$ / 检查交点是否在圆盘半径范围内
    Point3f pHit = ray(tShapeHit);
    Float dist2 = pHit.x * pHit.x + pHit.y * pHit.y;
    if (dist2 > radius * radius || dist2 < innerRadius * innerRadius)
        return false;

    // Test disk $\phi$ value against $\phimax$ / 测试角度phi是否在phimax范围内
    Float phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;
    if (phi > phiMax) return false;
    return true;
}

/**
 * @brief 计算圆盘面积(考虑内外半径和角度范围)
 */
Float Disk::Area() const {
    return phiMax * 0.5 * (radius * radius - innerRadius * innerRadius);
}

/**
 * @brief 在圆盘表面采样一个点
 * 使用同心圆盘采样(ConcentricSampleDisk)生成均匀分布
 */
Interaction Disk::Sample(const Point2f &u, Float *pdf) const {
    Point2f pd = ConcentricSampleDisk(u);
    Point3f pObj(pd.x * radius, pd.y * radius, height);
    Interaction it;
    it.n = Normalize((*ObjectToWorld)(Normal3f(0, 0, 1)));
    if (reverseOrientation) it.n *= -1;
    it.p = (*ObjectToWorld)(pObj, Vector3f(0, 0, 0), &it.pError);
    *pdf = 1 / Area();
    return it;
}

/**
 * @brief 创建圆盘形状的工厂函数
 */
std::shared_ptr<Disk> CreateDiskShape(const Transform *o2w,
                                      const Transform *w2o,
                                      bool reverseOrientation,
                                      const ParamSet &params) {
    Float height = params.FindOneFloat("height", 0.);
    Float radius = params.FindOneFloat("radius", 1);
    Float inner_radius = params.FindOneFloat("innerradius", 0);
    Float phimax = params.FindOneFloat("phimax", 360);
    return std::make_shared<Disk>(o2w, w2o, reverseOrientation, height, radius,
                                  inner_radius, phimax);
}

}  // namespace pbrt
