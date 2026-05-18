
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

#ifndef PBRT_CORE_QUATERNION_H
#define PBRT_CORE_QUATERNION_H

// core/quaternion.h*
// 四元数: 提供四元数表示和运算，用于平滑的旋转插值（球面线性插值Slerp），
// 在动画相机和变换中实现平滑旋转
#include "pbrt.h"
#include "stringprint.h"
#include "geometry.h"

namespace pbrt {

// Quaternion Declarations
// Quaternion: 四元数结构体，用于表示旋转和进行平滑旋转插值
struct Quaternion {
    // Quaternion Public Methods
    Quaternion() : v(0, 0, 0), w(1) {}  // 默认构造：单位四元数（无旋转）
    // 四元数加法
    Quaternion &operator+=(const Quaternion &q) {
        v += q.v;
        w += q.w;
        return *this;
    }
    friend Quaternion operator+(const Quaternion &q1, const Quaternion &q2) {
        Quaternion ret = q1;
        return ret += q2;
    }
    Quaternion &operator-=(const Quaternion &q) {
        v -= q.v;
        w -= q.w;
        return *this;
    }
    Quaternion operator-() const {
        Quaternion ret;
        ret.v = -v;
        ret.w = -w;
        return ret;
    }
    friend Quaternion operator-(const Quaternion &q1, const Quaternion &q2) {
        Quaternion ret = q1;
        return ret -= q2;
    }
    Quaternion &operator*=(Float f) {
        v *= f;
        w *= f;
        return *this;
    }
    Quaternion operator*(Float f) const {
        Quaternion ret = *this;
        ret.v *= f;
        ret.w *= f;
        return ret;
    }
    Quaternion &operator/=(Float f) {
        v /= f;
        w /= f;
        return *this;
    }
    Quaternion operator/(Float f) const {
        Quaternion ret = *this;
        ret.v /= f;
        ret.w /= f;
        return ret;
    }
    // ToTransform: 将四元数转换为对应的旋转变换矩阵
    Transform ToTransform() const;
    // 从变换矩阵构造四元数
    Quaternion(const Transform &t);

    friend std::ostream &operator<<(std::ostream &os, const Quaternion &q) {
        os << StringPrintf("[ %f, %f, %f, %f ]", q.v.x, q.v.y, q.v.z,
                           q.w);
        return os;
    }

    // Quaternion Public Data
    Vector3f v;  // 虚部向量 (x, y, z)
    Float w;     // 实部 (标量部分)
};

// Slerp: 球面线性插值，在两个四元数之间进行平滑旋转插值
Quaternion Slerp(Float t, const Quaternion &q1, const Quaternion &q2);

// Quaternion Inline Functions
// 标量乘以四元数
inline Quaternion operator*(Float f, const Quaternion &q) { return q * f; }

// Dot: 计算两个四元数的点积
inline Float Dot(const Quaternion &q1, const Quaternion &q2) {
    return Dot(q1.v, q2.v) + q1.w * q2.w;
}

// Normalize: 四元数归一化
inline Quaternion Normalize(const Quaternion &q) {
    return q / std::sqrt(Dot(q, q));
}

}  // namespace pbrt

#endif  // PBRT_CORE_QUATERNION_H
