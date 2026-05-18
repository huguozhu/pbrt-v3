
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

#ifndef PBRT_CAMERAS_PERSPECTIVE_H
#define PBRT_CAMERAS_PERSPECTIVE_H

// cameras/perspective.h*
// PerspectiveCamera: 透视相机，模拟标准针孔相机模型，
// 支持景深效果（透镜光圈）和透视投影，是最常用的相机类型
#include "pbrt.h"
#include "camera.h"
#include "film.h"

namespace pbrt {

// PerspectiveCamera Declarations
class PerspectiveCamera : public ProjectiveCamera {
  public:
    // PerspectiveCamera Public Methods
    // 构造函数：初始化透视投影矩阵、视场角(fov)等参数
    PerspectiveCamera(const AnimatedTransform &CameraToWorld,
                      const Bounds2f &screenWindow, Float shutterOpen,
                      Float shutterClose, Float lensRadius, Float focalDistance,
                      Float fov, Film *film, const Medium *medium);
    // GenerateRay: 根据采样点生成透视光线，从视点穿过像素位置
    Float GenerateRay(const CameraSample &sample, Ray *) const;
    // GenerateRayDifferential: 生成包含微分信息的光线，用于纹理抗锯齿
    Float GenerateRayDifferential(const CameraSample &sample,
                                  RayDifferential *ray) const;
    // We: 计算给定光线方向在相机上的重要性权重（用于路径追踪中的相机采样）
    Spectrum We(const Ray &ray, Point2f *pRaster2 = nullptr) const;
    // Pdf_We: 计算相机采样光线方向和位置的概率密度
    void Pdf_We(const Ray &ray, Float *pdfPos, Float *pdfDir) const;
    // Sample_Wi: 从相机对给定参考点采样入射方向（用于双向路径追踪）
    Spectrum Sample_Wi(const Interaction &ref, const Point2f &sample,
                       Vector3f *wi, Float *pdf, Point2f *pRaster,
                       VisibilityTester *vis) const;

  private:
    // PerspectiveCamera Private Data
    Vector3f dxCamera, dyCamera;		// ��x���y�᷽�����һ�����ص�λ�ã���Ӧ������ռ�ľ���
    Float A;
};

PerspectiveCamera *CreatePerspectiveCamera(const ParamSet &params,
                                           const AnimatedTransform &cam2world,
                                           Film *film, const Medium *medium);

}  // namespace pbrt

#endif  // PBRT_CAMERAS_PERSPECTIVE_H
