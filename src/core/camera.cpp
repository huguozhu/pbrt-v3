
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


// core/camera.cpp*
//
// 此文件实现了相机（Camera）基类的主要方法。
// Camera 是 pbrt 中所有相机类型的抽象基类，定义了从图像样本生成光线的接口，
// 以及相机对光场重要性采样的方法。具体相机类型（透视、正交、环境等）在子类中实现。
//

#include "camera.h"
#include "sampling.h"
#include "sampler.h"

namespace pbrt {

// Camera Method Definitions

// Camera 析构函数 —— 释放关联的 film 对象。
Camera::~Camera() { delete film; }

// Camera 构造函数。
// CameraToWorld: 从相机空间到世界空间的动画变换。
// shutterOpen/shutterClose: 快门开启和关闭时间。
// film: 与此相机关联的胶片对象。
// medium: 相机所在介质。
Camera::Camera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
               Float shutterClose, Film *film, const Medium *medium)
    : CameraToWorld(CameraToWorld),
      shutterOpen(shutterOpen),
      shutterClose(shutterClose),
      film(film),
      medium(medium) {
    if (CameraToWorld.HasScale())
        Warning(
            "Scaling detected in world-to-camera transformation!\n"
            "The system has numerous assumptions, implicit and explicit,\n"
            "that this transform will have no scale factors in it.\n"
            "Proceed at your own risk; your image may have errors or\n"
            "the system may crash as a result of this.");
}

// 生成具有差分信息的相机光线，用于计算纹理滤波的屏幕空间导数。
// 通过在 x 和 y 方向上各偏移半个像素来产生辅助光线，
// 从而估算光线原点和方向相对于像素位置的偏导数。
Float Camera::GenerateRayDifferential(const CameraSample &sample,
                                      RayDifferential *rd) const {
    Float wt = GenerateRay(sample, rd);
    if (wt == 0) return 0;

    // Find camera ray after shifting a fraction of a pixel in the $x$ direction
    Float wtx;
    for (Float eps : { .05, -.05 }) {
        CameraSample sshift = sample;
        sshift.pFilm.x += eps;
        Ray rx;
        wtx = GenerateRay(sshift, &rx);
        rd->rxOrigin = rd->o + (rx.o - rd->o) / eps;
        rd->rxDirection = rd->d + (rx.d - rd->d) / eps;
        if (wtx != 0)
            break;
    }
    if (wtx == 0)
        return 0;

    // Find camera ray after shifting a fraction of a pixel in the $y$ direction
    Float wty;
    for (Float eps : { .05, -.05 }) {
        CameraSample sshift = sample;
        sshift.pFilm.y += eps;
        Ray ry;
        wty = GenerateRay(sshift, &ry);
        rd->ryOrigin = rd->o + (ry.o - rd->o) / eps;
        rd->ryDirection = rd->d + (ry.d - rd->d) / eps;
        if (wty != 0)
            break;
    }
    if (wty == 0)
        return 0;

    rd->hasDifferentials = true;
    return wt;
}

// 计算相机沿给定方向的重要性（importance）We，即相机对场景辐射度的贡献权重。
// 此方法在基类中未实现，具体子类需要覆盖。
Spectrum Camera::We(const Ray &ray, Point2f *raster) const {
    LOG(FATAL) << "Camera::We() is not implemented!";
    return Spectrum(0.f);
}

// 计算相机重要性函数 We 对应的概率密度。
// pdfPos: 位置的概率密度, pdfDir: 方向的概率密度。
void Camera::Pdf_We(const Ray &ray, Float *pdfPos, Float *pdfDir) const {
    LOG(FATAL) << "Camera::Pdf_We() is not implemented!";
}

// 从参考点 ref 对相机方向进行重要性采样，用于双向路径追踪等算法。
// 采样来自相机的入射方向 wi，并计算该采样对应的概率密度和可见性测试对象。
Spectrum Camera::Sample_Wi(const Interaction &ref, const Point2f &u,
                           Vector3f *wi, Float *pdf, Point2f *pRaster,
                           VisibilityTester *vis) const {
    LOG(FATAL) << "Camera::Sample_Wi() is not implemented!";
    return Spectrum(0.f);
}

}  // namespace pbrt
