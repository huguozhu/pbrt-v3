
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


// shapes/cone.cpp*
/**
 * @file cone.cpp
 * @brief 圆锥体(Cone)几何体的实现
 *
 * 实现了圆锥体的构造、光线求交、面积计算等功能。
 * 圆锥体是沿Z轴延伸的二次曲面，其方程由参数radius, height和phiMax决定。
 * 光线-圆锥求交通过求解二次方程实现。
 */
#include "shapes/cone.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

namespace pbrt {

// Cone Method Definitions / 圆锥体方法实现
/**
 * @brief 圆锥体构造函数
 * 初始化圆锥体的几何参数，将角度限制在0-360度范围内并转换为弧度
 */
Cone::Cone(const Transform *o2w, const Transform *w2o, bool ro, Float height,
           Float radius, Float phiMax)
    : Shape(o2w, w2o, ro),
      radius(radius),
      height(height),
      phiMax(Radians(Clamp(phiMax, 0, 360))) {}
/**
 * @brief 计算圆锥体在对象空间中的轴对齐包围盒
 * 包围盒从圆锥底部(0)延伸到顶部(height)，半径范围[-radius, radius]
 */
Bounds3f Cone::ObjectBound() const {
    Point3f p1 = Point3f(-radius, -radius, 0);
    Point3f p2 = Point3f(radius, radius, height);
    return Bounds3f(p1, p2);
}

/**
 * @brief 光线-圆锥体求交(完整求交)
 *
 * 通过求解二次方程计算光线与圆锥面的交点。
 * 圆锥体的隐式方程为 (x^2 + y^2) = (radius/height)^2 * (z - height)^2
 *
 * @param r 光线
 * @param tHit 返回交点沿光线的距离
 * @param isect 返回交点表面信息
 * @param testAlphaTexture 是否测试Alpha纹理
 * @return 是否相交
 */
bool Cone::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                     bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersect);
    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space / 将光线变换到对象空间
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic cone coefficients / 计算圆锥二次型系数

    // Initialize _EFloat_ ray coordinate values / 用EFloat初始化光线坐标值(高精度)
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    // 计算圆锥斜率因子k = (radius/height)^2
    EFloat k = EFloat(radius) / EFloat(height);
    k = k * k;
    // 二次方程系数: a*t^2 + b*t + c = 0
    EFloat a = dx * dx + dy * dy - k * dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy - k * dz * (oz - height));
    EFloat c = ox * ox + oy * oy - k * (oz - height) * (oz - height);

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

    // Compute cone inverse mapping / 计算圆锥逆映射
    pHit = ray((Float)tShapeHit);
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0.) phi += 2 * Pi;

    // Test cone intersection against clipping parameters / 测试圆锥交点是否在裁剪参数范围内
    if (pHit.z < 0 || pHit.z > height || phi > phiMax) {
        if (tShapeHit == t1) return false;
        tShapeHit = t1;
        if (t1.UpperBound() > ray.tMax) return false;
        // Compute cone inverse mapping / 使用第二个交点重新计算
        pHit = ray((Float)tShapeHit);
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0.) phi += 2 * Pi;
        if (pHit.z < 0 || pHit.z > height || phi > phiMax) return false;
    }

    // Find parametric representation of cone hit / 求交点的参数化表示(u,v)
    Float u = phi / phiMax;
    Float v = pHit.z / height;

    // Compute cone $\dpdu$ and $\dpdv$ / 计算圆锥的偏导数向量 dpdu 和 dpdv
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv(-pHit.x / (1.f - v), -pHit.y / (1.f - v), height);

    // Compute cone $\dndu$ and $\dndv$ / 计算法线偏导数 dndu 和 dndv
    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0.);
    Vector3f d2Pduv = phiMax / (1.f - v) * Vector3f(pHit.y, -pHit.x, 0.);
    Vector3f d2Pdvv(0, 0, 0);

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

    // Compute error bounds for cone intersection / 计算圆锥交点的误差边界

    // Compute error bounds for intersection computed with ray equation / 基于光线方程计算误差边界
    EFloat px = ox + tShapeHit * dx;
    EFloat py = oy + tShapeHit * dy;
    EFloat pz = oz + tShapeHit * dz;
    Vector3f pError = Vector3f(px.GetAbsoluteError(), py.GetAbsoluteError(),
                               pz.GetAbsoluteError());

    // Initialize _SurfaceInteraction_ from parametric information / 用参数信息初始化SurfaceInteraction
    *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
                                                 -ray.d, dpdu, dpdv, dndu, dndv,
                                                 ray.time, this));
    *tHit = (Float)tShapeHit;
    return true;
}

/**
 * @brief 光线-圆锥体求交测试(仅判断是否相交，不返回交点细节)
 * 用于阴影光线等只需要判断遮挡关系的场景
 */
bool Cone::IntersectP(const Ray &r, bool testAlphaTexture) const {
    ProfilePhase p(Prof::ShapeIntersectP);
    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space / 将光线变换到对象空间
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic cone coefficients / 计算圆锥二次型系数

    // Initialize _EFloat_ ray coordinate values / 用EFloat初始化光线坐标值
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    EFloat k = EFloat(radius) / EFloat(height);
    k = k * k;
    EFloat a = dx * dx + dy * dy - k * dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy - k * dz * (oz - height));
    EFloat c = ox * ox + oy * oy - k * (oz - height) * (oz - height);

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

    // Compute cone inverse mapping / 计算圆锥逆映射
    pHit = ray((Float)tShapeHit);
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0.) phi += 2 * Pi;

    // Test cone intersection against clipping parameters / 测试圆锥交点是否在裁剪参数范围内
    if (pHit.z < 0 || pHit.z > height || phi > phiMax) {
        if (tShapeHit == t1) return false;
        tShapeHit = t1;
        if (t1.UpperBound() > ray.tMax) return false;
        // Compute cone inverse mapping / 使用第二个交点重新计算
        pHit = ray((Float)tShapeHit);
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0.) phi += 2 * Pi;
        if (pHit.z < 0 || pHit.z > height || phi > phiMax) return false;
    }
    return true;
}

/**
 * @brief 计算圆锥体表面积
 * 圆锥侧面积公式: pi * r * sqrt(h^2 + r^2)
 */
Float Cone::Area() const {
    return radius * std::sqrt((height * height) + (radius * radius)) * phiMax /
           2;
}

/**
 * @brief 在圆锥体表面采样一个点(未实现)
 */
Interaction Cone::Sample(const Point2f &u, Float *pdf) const {
    LOG(FATAL) << "Cone::Sample not implemented.";
    return Interaction();
}

/**
 * @brief 创建圆锥体形状的工厂函数
 * 从参数集中读取圆锥参数：半径、高度和最大角度
 */
std::shared_ptr<Cone> CreateConeShape(const Transform *o2w,
                                      const Transform *w2o,
                                      bool reverseOrientation,
                                      const ParamSet &params) {
    Float radius = params.FindOneFloat("radius", 1);
    Float height = params.FindOneFloat("height", 1);
    Float phimax = params.FindOneFloat("phimax", 360);
    return std::make_shared<Cone>(o2w, w2o, reverseOrientation, height, radius,
                                  phimax);
}

}  // namespace pbrt
