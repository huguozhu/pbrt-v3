
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

#ifndef PBRT_CORE_TEXTURE_H
#define PBRT_CORE_TEXTURE_H

// core/texture.h*
// 纹理系统模块：定义纹理映射和纹理求值的核心抽象基类。
// 支持2D纹理映射（UV映射、球面映射、柱面映射、平面映射）
// 和3D纹理映射，以及噪声函数（Value Noise、FBm、Turbulence）。
// 具体纹理实现（ConstantTexture、ScaleTexture、MixTexture等）位于src/textures/目录下。
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
#include "memory.h"

namespace pbrt {

// Texture Declarations
// 纹理类型声明

// TextureMapping2D: 2D纹理映射的抽象基类。
// 将表面交点映射到2D纹理坐标(s,t)，
// 同时计算纹理坐标在屏幕空间中的偏导数。
class TextureMapping2D {
  public:
    // TextureMapping2D Interface
    // TextureMapping2D 接口：将SurfaceInteraction映射到2D纹理坐标(s,t)
    virtual ~TextureMapping2D();
    virtual Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                        Vector2f *dstdy) const = 0;
};

// UVMapping2D: UV纹理映射。
// 直接使用表面的UV参数坐标，支持缩放和平移变换。
class UVMapping2D : public TextureMapping2D {
  public:
    // UVMapping2D Public Methods
    // UVMapping2D 公有方法
    UVMapping2D(Float su = 1, Float sv = 1, Float du = 0, Float dv = 0);
    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                Vector2f *dstdy) const;

  private:
	// �û�����Ķ�����ֵ�������ź�ƫ��
    const Float su, sv, du, dv;
};

// SphericalMapping2D: 球面纹理映射。
// 将3D点映射到球面上的(s,t)纹理坐标，常用于环境贴图。
class SphericalMapping2D : public TextureMapping2D {
  public:
    // SphericalMapping2D Public Methods
    // SphericalMapping2D 公有方法
    SphericalMapping2D(const Transform &WorldToTexture)
        : WorldToTexture(WorldToTexture) {}
    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                Vector2f *dstdy) const;

  private:
	// ��Բ������pת��ΪԲ��st��������
    Point2f sphere(const Point3f &P) const;
    const Transform WorldToTexture;
};

// CylindricalMapping2D: 柱面纹理映射。
// 将3D点映射到圆柱面上的(s,t)纹理坐标。
class CylindricalMapping2D : public TextureMapping2D {
  public:
    // CylindricalMapping2D Public Methods
    // CylindricalMapping2D 公有方法
    CylindricalMapping2D(const Transform &WorldToTexture)
        : WorldToTexture(WorldToTexture) {}
    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                Vector2f *dstdy) const;

  private:
    // CylindricalMapping2D Private Methods
    Point2f cylinder(const Point3f &p) const {
        Vector3f vec = Normalize(WorldToTexture(p) - Point3f(0, 0, 0));
        return Point2f((Pi + std::atan2(vec.y, vec.x)) * Inv2Pi, vec.z);
    }
    const Transform WorldToTexture;
};

// PlanarMapping2D: 平面纹理映射。
// 将3D点投影到平面上的(s,t)纹理坐标，使用两个向量定义投影平面。
class PlanarMapping2D : public TextureMapping2D {
  public:
    // PlanarMapping2D Public Methods
    // PlanarMapping2D 公有方法
    Point2f Map(const SurfaceInteraction &si, Vector2f *dstdx,
                Vector2f *dstdy) const;
    PlanarMapping2D(const Vector3f &vs, const Vector3f &vt, Float ds = 0,
                    Float dt = 0)
        : vs(vs), vt(vt), ds(ds), dt(dt) {}

  private:
    const Vector3f vs, vt;
    const Float ds, dt;
};

// TextureMapping3D: 3D纹理映射的抽象基类。
// 将表面交点映射到3D纹理坐标，用于3D噪声纹理和体纹理。
class TextureMapping3D {
  public:
    // TextureMapping3D Interface
    // TextureMapping3D 接口：将SurfaceInteraction映射到3D纹理坐标(p)
    virtual ~TextureMapping3D();
    virtual Point3f Map(const SurfaceInteraction &si, Vector3f *dpdx,
                        Vector3f *dpdy) const = 0;
};

// IdentityMapping3D: 3D标识纹理映射。
// 直接将世界空间或物体空间的3D坐标作为纹理坐标，
// 常被ConstantTexture(常量纹理)、ScaleTexture(缩放纹理)、
// MixTexture(混合纹理)等间接使用。
class IdentityMapping3D : public TextureMapping3D {
  public:
    // IdentityMapping3D Public Methods
    // IdentityMapping3D 公有方法
    IdentityMapping3D(const Transform &WorldToTexture)
        : WorldToTexture(WorldToTexture) {}
    Point3f Map(const SurfaceInteraction &si, Vector3f *dpdx,
                Vector3f *dpdy) const;

  private:
    const Transform WorldToTexture;
};

// Texture: 纹理抽象基类模板。
// 定义了所有纹理类型的核心接口，所有具体纹理实现
//（如ConstantTexture常量纹理、ScaleTexture缩放纹理、
//  MixTexture混合纹理、ImageTexture图像纹理等）
// 均继承此类并实现Evaluate方法。
template <typename T>
class Texture {
  public:
    // Texture Interface
    // Texture 接口：在给定的表面交点上求取纹理值
    virtual T Evaluate(const SurfaceInteraction &) const = 0;
    virtual ~Texture() {}
};

// Lanczos: Lanczos sinc函数重采样滤波器，用于纹理重建和图像缩放
Float Lanczos(Float, Float tau = 2);
// Noise: Perlin噪声函数（值噪声），生成连续的随机噪声值
Float Noise(Float x, Float y = .5f, Float z = .5f);
Float Noise(const Point3f &p);
// FBm: 分形布朗运动噪声，通过叠加多个频率的噪声生成自然纹理
Float FBm(const Point3f &p, const Vector3f &dpdx, const Vector3f &dpdy,
          Float omega, int octaves);
// Turbulence: 湍流噪声，取噪声音频分量的绝对值之和，产生湍流效果
Float Turbulence(const Point3f &p, const Vector3f &dpdx, const Vector3f &dpdy,
                 Float omega, int octaves);

}  // namespace pbrt

#endif  // PBRT_CORE_TEXTURE_H
