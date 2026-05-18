
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

#ifndef PBRT_CORE_TRANSFORM_H
#define PBRT_CORE_TRANSFORM_H

// core/transform.h*
// 变换系统模块：定义4x4矩阵和齐次坐标变换的核心抽象。
// 支持平移、旋转、缩放、透视投影等仿射变换和投影变换，
// 以及支持时间域插值的动画变换（AnimatedTransform）。
// 提供带误差上界估计的安全变换操作。
#include "pbrt.h"
#include "stringprint.h"
#include "geometry.h"
#include "quaternion.h"

namespace pbrt {

// Matrix4x4 Declarations
// 矩阵类型声明

// Matrix4x4: 4x4矩阵结构体。
// 用于表示齐次坐标变换矩阵，支持矩阵乘法、转置、求逆等基本运算。
// 行优先存储，m[行][列]。
struct Matrix4x4 {
    // Matrix4x4 Public Methods
    // Matrix4x4 公有方法
    // 默认构造函数：初始化为4x4单位矩阵
    Matrix4x4() {
        m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.f;
        m[0][1] = m[0][2] = m[0][3] = m[1][0] = m[1][2] = m[1][3] = m[2][0] =
            m[2][1] = m[2][3] = m[3][0] = m[3][1] = m[3][2] = 0.f;
    }
    Matrix4x4(Float mat[4][4]);
    Matrix4x4(Float t00, Float t01, Float t02, Float t03, Float t10, Float t11,
              Float t12, Float t13, Float t20, Float t21, Float t22, Float t23,
              Float t30, Float t31, Float t32, Float t33);
    bool operator==(const Matrix4x4 &m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return false;
        return true;
    }
    bool operator!=(const Matrix4x4 &m2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (m[i][j] != m2.m[i][j]) return true;
        return false;
    }
    // Transpose: 矩阵转置（友元函数）
    friend Matrix4x4 Transpose(const Matrix4x4 &);
    // Print: 打印矩阵内容到文件流
    void Print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < 4; ++i) {
            fprintf(f, "  [ ");
            for (int j = 0; j < 4; ++j) {
                fprintf(f, "%f", m[i][j]);
                if (j != 3) fprintf(f, ", ");
            }
            fprintf(f, " ]\n");
        }
        fprintf(f, " ] ");
    }
    // Mul: 矩阵乘法（静态方法）
    static Matrix4x4 Mul(const Matrix4x4 &m1, const Matrix4x4 &m2) {
        Matrix4x4 r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
                            m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
        return r;
    }
    // Inverse: 矩阵求逆（友元函数）
    friend Matrix4x4 Inverse(const Matrix4x4 &);

    // operator<<: 矩阵输出到输出流
    friend std::ostream &operator<<(std::ostream &os, const Matrix4x4 &m) {
        // clang-format off
        os << StringPrintf("[ [ %f, %f, %f, %f ] "
                           "[ %f, %f, %f, %f ] "
                           "[ %f, %f, %f, %f ] "
                           "[ %f, %f, %f, %f ] ]",
                           m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3],
                           m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3],
                           m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3],
                           m.m[3][0], m.m[3][1], m.m[3][2], m.m[3][3]);
        // clang-format on
        return os;
    }

    Float m[4][4];  // 4x4矩阵数据，行优先存储
};

// Transform Declarations
// 变换类型声明
// Transform: 齐次坐标变换类。
// 封装了4x4变换矩阵及其逆矩阵，提供对点、向量、法线、光线等
// 几何实体的变换操作。支持变换组合、手性判断和误差分析。
class Transform {
  public:
    // Transform Public Methods
    // Transform 公有方法
    // 默认构造函数：初始化为单位变换
    Transform() {}
    // 从4x4矩阵构造变换（自动计算逆矩阵）
    Transform(const Float mat[4][4]) {
        m = Matrix4x4(mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0],
                      mat[1][1], mat[1][2], mat[1][3], mat[2][0], mat[2][1],
                      mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2],
                      mat[3][3]);
        mInv = Inverse(m);
    }
    // 从Matrix4x4构造变换（自动计算逆矩阵）
    Transform(const Matrix4x4 &m) : m(m), mInv(Inverse(m)) {}
    // 从变换矩阵和逆矩阵构造变换
    Transform(const Matrix4x4 &m, const Matrix4x4 &mInv) : m(m), mInv(mInv) {}
    // Print: 打印变换矩阵到文件流
    void Print(FILE *f) const;
    // Inverse: 求变换的逆变换（通过交换m和mInv实现）
    friend Transform Inverse(const Transform &t) {
        return Transform(t.mInv, t.m);
    }
    // Transpose: 求变换的转置变换
    friend Transform Transpose(const Transform &t) {
        return Transform(Transpose(t.m), Transpose(t.mInv));
    }
    bool operator==(const Transform &t) const {
        return t.m == m && t.mInv == mInv;
    }
    bool operator!=(const Transform &t) const {
        return t.m != m || t.mInv != mInv;
    }
    bool operator<(const Transform &t2) const {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                if (m.m[i][j] < t2.m.m[i][j]) return true;
                if (m.m[i][j] > t2.m.m[i][j]) return false;
            }
        return false;
    }
    // IsIdentity: 判断变换是否为单位变换
    bool IsIdentity() const {
        return (m.m[0][0] == 1.f && m.m[0][1] == 0.f && m.m[0][2] == 0.f &&
                m.m[0][3] == 0.f && m.m[1][0] == 0.f && m.m[1][1] == 1.f &&
                m.m[1][2] == 0.f && m.m[1][3] == 0.f && m.m[2][0] == 0.f &&
                m.m[2][1] == 0.f && m.m[2][2] == 1.f && m.m[2][3] == 0.f &&
                m.m[3][0] == 0.f && m.m[3][1] == 0.f && m.m[3][2] == 0.f &&
                m.m[3][3] == 1.f);
    }
    // GetMatrix: 获取变换矩阵
    const Matrix4x4 &GetMatrix() const { return m; }
    // GetInverseMatrix: 获取逆变换矩阵
    const Matrix4x4 &GetInverseMatrix() const { return mInv; }
    // HasScale: 检查变换是否包含缩放（单位向量变换后长度是否改变）
    bool HasScale() const {
        Float la2 = (*this)(Vector3f(1, 0, 0)).LengthSquared();
        Float lb2 = (*this)(Vector3f(0, 1, 0)).LengthSquared();
        Float lc2 = (*this)(Vector3f(0, 0, 1)).LengthSquared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
        return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
    }
    // operator(): 变换点（使用齐次坐标，自动处理w分量）
    template <typename T>
    inline Point3<T> operator()(const Point3<T> &p) const;
    // operator(): 变换向量（忽略平移部分，仅使用左上3x3子矩阵）
    template <typename T>
    inline Vector3<T> operator()(const Vector3<T> &v) const;
    // operator(): 变换法线（使用逆矩阵的转置，保证法线方向正确）
    template <typename T>
    inline Normal3<T> operator()(const Normal3<T> &) const;
    // operator(): 变换光线（包括原点偏移以避免自相交）
    inline Ray operator()(const Ray &r) const;
    // operator(): 变换差分光线
    inline RayDifferential operator()(const RayDifferential &r) const;
    // operator(): 变换包围盒
    Bounds3f operator()(const Bounds3f &b) const;
    // operator*: 变换组合（左乘）
    Transform operator*(const Transform &t2) const;
    // SwapsHandedness: 判断变换是否改变了坐标系手性（左右手系互换）
    bool SwapsHandedness() const;
    // operator(): 变换表面交互点
    SurfaceInteraction operator()(const SurfaceInteraction &si) const;
    // operator(): 变换点并计算绝对误差上界
    template <typename T>
    inline Point3<T> operator()(const Point3<T> &pt,
                                Vector3<T> *absError) const;
    // operator(): 变换带误差的点，并计算传播后的误差
    template <typename T>
    inline Point3<T> operator()(const Point3<T> &p, const Vector3<T> &pError,
                                Vector3<T> *pTransError) const;
    // operator(): 变换向量并计算绝对误差上界
    template <typename T>
    inline Vector3<T> operator()(const Vector3<T> &v,
                                 Vector3<T> *vTransError) const;
    // operator(): 变换带误差的向量，并计算传播后的误差
    template <typename T>
    inline Vector3<T> operator()(const Vector3<T> &v, const Vector3<T> &vError,
                                 Vector3<T> *vTransError) const;
    // operator(): 变换光线，计算原点偏移误差
    inline Ray operator()(const Ray &r, Vector3f *oError,
                          Vector3f *dError) const;
    inline Ray operator()(const Ray &r, const Vector3f &oErrorIn,
                          const Vector3f &dErrorIn, Vector3f *oErrorOut,
                          Vector3f *dErrorOut) const;

    friend std::ostream &operator<<(std::ostream &os, const Transform &t) {
        os << "t=" << t.m << ", inv=" << t.mInv;
        return os;
    }

  private:
    // Transform Private Data
    // Transform 私有数据
    Matrix4x4 m;     // 变换矩阵
    Matrix4x4 mInv;  // 逆变换矩阵
    friend class AnimatedTransform;
    friend struct Quaternion;
};

// Translate: 创建平移变换
Transform Translate(const Vector3f &delta);
// Scale: 创建缩放变换
Transform Scale(Float x, Float y, Float z);
// RotateX: 创建绕X轴旋转变换（角度以弧度为单位）
Transform RotateX(Float theta);
// RotateY: 创建绕Y轴旋转变换（角度以弧度为单位）
Transform RotateY(Float theta);
// RotateZ: 创建绕Z轴旋转变换（角度以弧度为单位）
Transform RotateZ(Float theta);
// Rotate: 创建绕任意轴旋转变换
Transform Rotate(Float theta, const Vector3f &axis);
// LookAt: 创建观察变换（从pos位置看向look目标点）
Transform LookAt(const Point3f &pos, const Point3f &look, const Vector3f &up);
// Orthographic: 创建正交投影变换
Transform Orthographic(Float znear, Float zfar);
// Perspective: 创建透视投影变换
Transform Perspective(Float fov, Float znear, Float zfar);
bool SolveLinearSystem2x2(const Float A[2][2], const Float B[2], Float *x0,
                          Float *x1);

// Transform Inline Functions
// 变换内联函数实现

// operator(): 变换点（齐次坐标变换，自动除以w分量）
template <typename T>
inline Point3<T> Transform::operator()(const Point3<T> &p) const {
    T x = p.x, y = p.y, z = p.z;
    T xp = m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z + m.m[0][3];
    T yp = m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z + m.m[1][3];
    T zp = m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z + m.m[2][3];
    T wp = m.m[3][0] * x + m.m[3][1] * y + m.m[3][2] * z + m.m[3][3];
    CHECK_NE(wp, 0);
    if (wp == 1)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

// operator(): 变换向量（仅使用3x3子矩阵，忽略平移）
template <typename T>
inline Vector3<T> Transform::operator()(const Vector3<T> &v) const {
    T x = v.x, y = v.y, z = v.z;
    return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                      m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                      m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

// operator(): 变换法线（使用逆矩阵转置，保证法线与切线的正交性）
template <typename T>
inline Normal3<T> Transform::operator()(const Normal3<T> &n) const {
    T x = n.x, y = n.y, z = n.z;
    return Normal3<T>(mInv.m[0][0] * x + mInv.m[1][0] * y + mInv.m[2][0] * z,
                      mInv.m[0][1] * x + mInv.m[1][1] * y + mInv.m[2][1] * z,
                      mInv.m[0][2] * x + mInv.m[1][2] * y + mInv.m[2][2] * z);
}

// operator(): 变换光线（包含原点误差偏移以防止自相交）
inline Ray Transform::operator()(const Ray &r) const {
    Vector3f oError;
    Point3f o = (*this)(r.o, &oError);
    Vector3f d = (*this)(r.d);
    // Offset ray origin to edge of error bounds and compute _tMax_
    Float lengthSquared = d.LengthSquared();
    Float tMax = r.tMax;
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), oError) / lengthSquared;
        o += d * dt;
        tMax -= dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

// operator(): 变换差分光线（包含主光线和辅助光线的变换）
inline RayDifferential Transform::operator()(const RayDifferential &r) const {
    Ray tr = (*this)(Ray(r));
    RayDifferential ret(tr.o, tr.d, tr.tMax, tr.time, tr.medium);
    ret.hasDifferentials = r.hasDifferentials;
    ret.rxOrigin = (*this)(r.rxOrigin);
    ret.ryOrigin = (*this)(r.ryOrigin);
    ret.rxDirection = (*this)(r.rxDirection);
    ret.ryDirection = (*this)(r.ryDirection);
    return ret;
}

// operator(): 变换点并计算绝对误差上界（使用gamma(3)误差分析）
template <typename T>
inline Point3<T> Transform::operator()(const Point3<T> &p,
                                       Vector3<T> *pError) const {
    T x = p.x, y = p.y, z = p.z;
    // Compute transformed coordinates from point _pt_
    T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);

    // Compute absolute error for transformed point
    T xAbsSum = (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
                 std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    T yAbsSum = (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
                 std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    T zAbsSum = (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
                 std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    *pError = gamma(3) * Vector3<T>(xAbsSum, yAbsSum, zAbsSum);
    CHECK_NE(wp, 0);
    if (wp == 1)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

// operator(): 变换带输入误差的点，计算误差传播后的新误差上界
template <typename T>
inline Point3<T> Transform::operator()(const Point3<T> &pt,
                                       const Vector3<T> &ptError,
                                       Vector3<T> *absError) const {
    T x = pt.x, y = pt.y, z = pt.z;
    T xp = (m.m[0][0] * x + m.m[0][1] * y) + (m.m[0][2] * z + m.m[0][3]);
    T yp = (m.m[1][0] * x + m.m[1][1] * y) + (m.m[1][2] * z + m.m[1][3]);
    T zp = (m.m[2][0] * x + m.m[2][1] * y) + (m.m[2][2] * z + m.m[2][3]);
    T wp = (m.m[3][0] * x + m.m[3][1] * y) + (m.m[3][2] * z + m.m[3][3]);
    absError->x =
        (gamma(3) + (T)1) *
            (std::abs(m.m[0][0]) * ptError.x + std::abs(m.m[0][1]) * ptError.y +
             std::abs(m.m[0][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[0][0] * x) + std::abs(m.m[0][1] * y) +
                    std::abs(m.m[0][2] * z) + std::abs(m.m[0][3]));
    absError->y =
        (gamma(3) + (T)1) *
            (std::abs(m.m[1][0]) * ptError.x + std::abs(m.m[1][1]) * ptError.y +
             std::abs(m.m[1][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[1][0] * x) + std::abs(m.m[1][1] * y) +
                    std::abs(m.m[1][2] * z) + std::abs(m.m[1][3]));
    absError->z =
        (gamma(3) + (T)1) *
            (std::abs(m.m[2][0]) * ptError.x + std::abs(m.m[2][1]) * ptError.y +
             std::abs(m.m[2][2]) * ptError.z) +
        gamma(3) * (std::abs(m.m[2][0] * x) + std::abs(m.m[2][1] * y) +
                    std::abs(m.m[2][2] * z) + std::abs(m.m[2][3]));
    CHECK_NE(wp, 0);
    if (wp == 1.)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

// operator(): 变换向量并计算绝对误差上界
template <typename T>
inline Vector3<T> Transform::operator()(const Vector3<T> &v,
                                        Vector3<T> *absError) const {
    T x = v.x, y = v.y, z = v.z;
    absError->x =
        gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) +
                    std::abs(m.m[0][2] * v.z));
    absError->y =
        gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) +
                    std::abs(m.m[1][2] * v.z));
    absError->z =
        gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) +
                    std::abs(m.m[2][2] * v.z));
    return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                      m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                      m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

// operator(): 变换带输入误差的向量，计算误差传播后的新误差上界
template <typename T>
inline Vector3<T> Transform::operator()(const Vector3<T> &v,
                                        const Vector3<T> &vError,
                                        Vector3<T> *absError) const {
    T x = v.x, y = v.y, z = v.z;
    absError->x =
        (gamma(3) + (T)1) *
            (std::abs(m.m[0][0]) * vError.x + std::abs(m.m[0][1]) * vError.y +
             std::abs(m.m[0][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[0][0] * v.x) + std::abs(m.m[0][1] * v.y) +
                    std::abs(m.m[0][2] * v.z));
    absError->y =
        (gamma(3) + (T)1) *
            (std::abs(m.m[1][0]) * vError.x + std::abs(m.m[1][1]) * vError.y +
             std::abs(m.m[1][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[1][0] * v.x) + std::abs(m.m[1][1] * v.y) +
                    std::abs(m.m[1][2] * v.z));
    absError->z =
        (gamma(3) + (T)1) *
            (std::abs(m.m[2][0]) * vError.x + std::abs(m.m[2][1]) * vError.y +
             std::abs(m.m[2][2]) * vError.z) +
        gamma(3) * (std::abs(m.m[2][0] * v.x) + std::abs(m.m[2][1] * v.y) +
                    std::abs(m.m[2][2] * v.z));
    return Vector3<T>(m.m[0][0] * x + m.m[0][1] * y + m.m[0][2] * z,
                      m.m[1][0] * x + m.m[1][1] * y + m.m[1][2] * z,
                      m.m[2][0] * x + m.m[2][1] * y + m.m[2][2] * z);
}

// operator(): 变换光线并计算原点和方向误差
inline Ray Transform::operator()(const Ray &r, Vector3f *oError,
                                 Vector3f *dError) const {
    Point3f o = (*this)(r.o, oError);
    Vector3f d = (*this)(r.d, dError);
    Float tMax = r.tMax;
    Float lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), *oError) / lengthSquared;
        o += d * dt;
        //        tMax -= dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

// operator(): 变换带输入误差的光线，计算输出误差
inline Ray Transform::operator()(const Ray &r, const Vector3f &oErrorIn,
                                 const Vector3f &dErrorIn, Vector3f *oErrorOut,
                                 Vector3f *dErrorOut) const {
    Point3f o = (*this)(r.o, oErrorIn, oErrorOut);
    Vector3f d = (*this)(r.d, dErrorIn, dErrorOut);
    Float tMax = r.tMax;
    Float lengthSquared = d.LengthSquared();
    if (lengthSquared > 0) {
        Float dt = Dot(Abs(d), *oErrorOut) / lengthSquared;
        o += d * dt;
        //        tMax -= dt;
    }
    return Ray(o, d, tMax, r.time, r.medium);
}

// AnimatedTransform Declarations
// 动画变换类型声明

// AnimatedTransform: 动画变换类。
// 表示在时间区间[startTime, endTime]内从startTransform到endTransform
// 的插值变换。通过矩阵分解（平移、旋转、缩放）实现平滑插值，
// 用于运动模糊渲染。
class AnimatedTransform {
  public:
    // AnimatedTransform Public Methods
    // AnimatedTransform 公有方法
    // 构造函数：指定起始/结束变换及其对应时间点
    AnimatedTransform(const Transform *startTransform, Float startTime,
                      const Transform *endTransform, Float endTime);
    // Decompose: 将变换矩阵分解为平移(T)、旋转(R)和缩放(S)分量
    static void Decompose(const Matrix4x4 &m, Vector3f *T, Quaternion *R,
                          Matrix4x4 *S);
    // Interpolate: 在指定时间点插值得到中间变换
    void Interpolate(Float time, Transform *t) const;
    // operator(): 在插值时间点变换光线
    Ray operator()(const Ray &r) const;
    // operator(): 在插值时间点变换差分光线
    RayDifferential operator()(const RayDifferential &r) const;
    // operator(): 在指定时间点变换点
    Point3f operator()(Float time, const Point3f &p) const;
    // operator(): 在指定时间点变换向量
    Vector3f operator()(Float time, const Vector3f &v) const;
    // HasScale: 检查动画变换是否包含缩放
    bool HasScale() const {
        return startTransform->HasScale() || endTransform->HasScale();
    }
    // MotionBounds: 计算运动模糊包围盒
    Bounds3f MotionBounds(const Bounds3f &b) const;
    // BoundPointMotion: 计算单点运动轨迹的包围盒
    Bounds3f BoundPointMotion(const Point3f &p) const;

  private:
    // AnimatedTransform Private Data
    // AnimatedTransform 私有数据
    const Transform *startTransform, *endTransform;  // 起始和结束变换指针
    const Float startTime, endTime;          // 起始和结束时间
    const bool actuallyAnimated;              // 是否真正有动画效果
    Vector3f T[2];         // 平移分量（两个时间点）
    Quaternion R[2];       // 旋转分量（四元数，两个时间点）
    Matrix4x4 S[2];        // 缩放分量（两个时间点）
    bool hasRotation;      // 是否包含旋转分量
    // DerivativeTerm: 用于运动边界计算的导数项
    struct DerivativeTerm {
        DerivativeTerm() {}
        DerivativeTerm(Float c, Float x, Float y, Float z)
            : kc(c), kx(x), ky(y), kz(z) {}
        Float kc, kx, ky, kz;
        Float Eval(const Point3f &p) const {
            return kc + kx * p.x + ky * p.y + kz * p.z;
        }
    };
    DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];  // 运动包围盒计算系数数组
};

}  // namespace pbrt

#endif  // PBRT_CORE_TRANSFORM_H
