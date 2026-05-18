
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


// shapes/cylinder.cpp*
/**
 * @file cylinder.cpp
 * @brief 圆柱体(Cylinder)几何体的实现
 *
 * 实现了圆柱体的光线求交、面积计算和采样功能。
 * 圆柱体外表面方程为 x^2 + y^2 = radius^2, z在[zMin, zMax]之间。
 */
#include "shapes/cylinder.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

namespace pbrt {

// Cylinder Method Definitions / 圆柱体方法实现
/**
 * @brief 计算圆柱体在对象空间的包围盒
 */
Bounds3f Cylinder::ObjectBound() const {
    return Bounds3f(Point3f(-radius, -radius, zMin),
                    Point3f(radius, radius, zMax));
}

/**
 * @brief 光线-圆柱体求交(完整求交)
 *
 * 在对象空间求解光线与圆柱面的二次方程。
 * 圆柱隐式方程: x^2 + y^2 - radius^2 = 0, 同时z在[zMin, zMax]内
 */
bool Cylinder::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                         bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersect);
    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space / 将光线变换到对象空间
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic cylinder coefficients / 计算圆柱体二次型系数

    // Initialize _EFloat_ ray coordinate values / 用EFloat初始化光线坐标值
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    EFloat a = dx * dx + dy * dy;
    EFloat b = 2 * (dx * ox + dy * oy);
    EFloat c = ox * ox + oy * oy - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values / 求解二次方程得到t值
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection / 选择最近的合理交点
    if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) return false;
    }

    // Compute cylinder hit point and $\phi$ / 计算圆柱体交点和角度phi
    pHit = ray((Float)tShapeHit);

    // Refine cylinder intersection point / 精修圆柱体交点
    Float hitRad = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
    pHit.x *= radius / hitRad;
    pHit.y *= radius / hitRad;
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;

    // Test cylinder intersection against clipping parameters / 测试圆柱体交点是否在裁剪参数范围内
    if (pHit.z < zMin || pHit.z > zMax || phi > phiMax) {
        if (tShapeHit == t1) return false;
        tShapeHit = t1;
        if (t1.UpperBound() > ray.tMax) return false;
        // Compute cylinder hit point and $\phi$ / 计算圆柱体交点和角度phi
        pHit = ray((Float)tShapeHit);

        // Refine cylinder intersection point / 精修圆柱体交点
        Float hitRad = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
        pHit.x *= radius / hitRad;
        pHit.y *= radius / hitRad;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;
        if (pHit.z < zMin || pHit.z > zMax || phi > phiMax) return false;
    }

    // Find parametric representation of cylinder hit / 计算交点的参数化表示(u,v)
    Float u = phi / phiMax;
    Float v = (pHit.z - zMin) / (zMax - zMin);

    // Compute cylinder $\dpdu$ and $\dpdv$ / 计算圆柱体偏导数向量
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv(0, 0, zMax - zMin);

    // Compute cylinder $\dndu$ and $\dndv$ / 计算圆柱体法线偏导数
    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
    Vector3f d2Pduv(0, 0, 0), d2Pdvv(0, 0, 0);

    // Compute coefficients for fundamental forms
    Float E = Dot(dpdu, dpdu);
    Float F = Dot(dpdu, dpdv);
    Float G = Dot(dpdv, dpdv);
    Vector3f N = Normalize(Cross(dpdu, dpdv));
    Float e = Dot(N, d2Pduu);
    Float f = Dot(N, d2Pduv);
    Float g = Dot(N, d2Pdvv);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    Float invEGF2 = 1 / (E * G - F * F);
    Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu +
                             (e * F - f * E) * invEGF2 * dpdv);
    Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu +
                             (f * F - g * E) * invEGF2 * dpdv);

    // Compute error bounds for cylinder intersection / 计算圆柱体交点误差边界
    Vector3f pError = gamma(3) * Abs(Vector3f(pHit.x, pHit.y, 0));

    // Initialize _SurfaceInteraction_ from parametric information / 初始化表面交互信息
    *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
                                                 -ray.d, dpdu, dpdv, dndu, dndv,
                                                 ray.time, this));

    // Update _tHit_ for quadric intersection
    *tHit = (Float)tShapeHit;
    return true;
}

/**
 * @brief 光线-圆柱体求交测试(仅判断是否相交)
 */
bool Cylinder::IntersectP(const Ray &r, bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersectP);
    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space / 将光线变换到对象空间
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic cylinder coefficients / 计算圆柱体二次型系数

    // Initialize _EFloat_ ray coordinate values / 用EFloat初始化光线坐标值
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    EFloat a = dx * dx + dy * dy;
    EFloat b = 2 * (dx * ox + dy * oy);
    EFloat c = ox * ox + oy * oy - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values / 求解二次方程得到t值
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection / 选择最近的合理交点
    if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) return false;
    }

    // Compute cylinder hit point and $\phi$ / 计算圆柱体交点和角度phi
    pHit = ray((Float)tShapeHit);

    // Refine cylinder intersection point / 精修圆柱体交点
    Float hitRad = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
    pHit.x *= radius / hitRad;
    pHit.y *= radius / hitRad;
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;

    // Test cylinder intersection against clipping parameters / 测试圆柱体交点是否在裁剪参数范围内
    if (pHit.z < zMin || pHit.z > zMax || phi > phiMax) {
        if (tShapeHit == t1) return false;
        tShapeHit = t1;
        if (t1.UpperBound() > ray.tMax) return false;
        // Compute cylinder hit point and $\phi$ / 计算圆柱体交点和角度phi
        pHit = ray((Float)tShapeHit);

        // Refine cylinder intersection point / 精修圆柱体交点
        Float hitRad = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
        pHit.x *= radius / hitRad;
        pHit.y *= radius / hitRad;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;
        if (pHit.z < zMin || pHit.z > zMax || phi > phiMax) return false;
    }
    return true;
}

/**
 * @brief 计算圆柱体表面积
 */
Float Cylinder::Area() const { return (zMax - zMin) * radius * phiMax; }

/**
 * @brief 在圆柱体表面采样一个点
 * 通过z值和角度phi进行均匀采样
 */
Interaction Cylinder::Sample(const Point2f &u, Float *pdf) const {
    Float z = Lerp(u[0], zMin, zMax);
    Float phi = u[1] * phiMax;
    Point3f pObj = Point3f(radius * std::cos(phi), radius * std::sin(phi), z);
    Interaction it;
    it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, 0)));
    if (reverseOrientation) it.n *= -1;
    // Reproject _pObj_ to cylinder surface and compute _pObjError_
    Float hitRad = std::sqrt(pObj.x * pObj.x + pObj.y * pObj.y);
    pObj.x *= radius / hitRad;
    pObj.y *= radius / hitRad;
    Vector3f pObjError = gamma(3) * Abs(Vector3f(pObj.x, pObj.y, 0));
    it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
    *pdf = 1 / Area();
    return it;
}

/**
 * @brief 创建圆柱体形状的工厂函数
 */
std::shared_ptr<Cylinder> CreateCylinderShape(const Transform *o2w,
                                              const Transform *w2o,
                                              bool reverseOrientation,
                                              const ParamSet &params) {
    Float radius = params.FindOneFloat("radius", 1);
    Float zmin = params.FindOneFloat("zmin", -1);
    Float zmax = params.FindOneFloat("zmax", 1);
    Float phimax = params.FindOneFloat("phimax", 360);
    return std::make_shared<Cylinder>(o2w, w2o, reverseOrientation, radius,
                                      zmin, zmax, phimax);
}

}  // namespace pbrt
