
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

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

// cameras/realistic.h*
// RealisticCamera: 真实相机模型，使用多片透镜组模拟真实物理相机，
// 支持光圈直径、对焦距离和透镜系统的光线追踪，实现逼真的光学效果
#include "pbrt.h"
#include "camera.h"
#include "film.h"

namespace pbrt {

// RealisticCamera Declarations
class RealisticCamera : public Camera {
  public:
    // RealisticCamera Public Methods
    // 构造函数：初始化透镜组、光圈直径、对焦距离等参数
    RealisticCamera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
                    Float shutterClose, Float apertureDiameter,
                    Float focusDistance, bool simpleWeighting,
                    std::vector<Float> &lensData, Film *film,
                    const Medium *medium);
    // GenerateRay: 对给定采样点，通过透镜系统追踪生成光线
    Float GenerateRay(const CameraSample &sample, Ray *) const;

  private:
    // RealisticCamera Private Declarations
    // LensElementInterface: 透镜组中单个透镜元件的接口参数
    struct LensElementInterface {
        Float curvatureRadius;  // 透镜曲率半径
        Float thickness;        // 透镜厚度
        Float eta;              // 折射率
        Float apertureRadius;   // 光圈半径
    };

    // RealisticCamera Private Data
    const bool simpleWeighting;                              // 是否使用简单重要性权重
    std::vector<LensElementInterface> elementInterfaces;       // 透镜组元件列表
    std::vector<Bounds2f> exitPupilBounds;                     // 出射光瞳边界

    // RealisticCamera Private Methods
    // LensRearZ: 计算透镜组背面（靠近胶片）的Z坐标
    Float LensRearZ() const { return elementInterfaces.back().thickness; }
    // LensFrontZ: 计算透镜组正面（靠近场景）的Z坐标
    Float LensFrontZ() const {
        Float zSum = 0;
        for (const LensElementInterface &element : elementInterfaces)
            zSum += element.thickness;
        return zSum;
    }
    // RearElementRadius: 获取最后一个透镜元件的光圈半径
    Float RearElementRadius() const {
        return elementInterfaces.back().apertureRadius;
    }
    // TraceLensesFromFilm: 从胶片侧向场景侧追踪光线通过透镜系统
    bool TraceLensesFromFilm(const Ray &ray, Ray *rOut) const;
    // IntersectSphericalElement: 计算光线与球面透镜元件的交点
    static bool IntersectSphericalElement(Float radius, Float zCenter,
                                          const Ray &ray, Float *t,
                                          Normal3f *n);
    // TraceLensesFromScene: 从场景侧向胶片侧追踪光线通过透镜系统
    bool TraceLensesFromScene(const Ray &rCamera, Ray *rOut) const;
    // 以下是透镜系统可视化调试方法
    void DrawLensSystem() const;
    void DrawRayPathFromFilm(const Ray &r, bool arrow,
                             bool toOpticalIntercept) const;
    void DrawRayPathFromScene(const Ray &r, bool arrow,
                              bool toOpticalIntercept) const;
    // ComputeCardinalPoints: 计算透镜组的基点（主点和焦点）
    static void ComputeCardinalPoints(const Ray &rIn, const Ray &rOut, Float *p,
                                      Float *f);
    // ComputeThickLensApproximation: 计算厚透镜近似参数
    void ComputeThickLensApproximation(Float pz[2], Float f[2]) const;
    // 以下是对焦方法
    Float FocusThickLens(Float focusDistance);       // 使用厚透镜公式自动对焦
    Float FocusBinarySearch(Float focusDistance);    // 使用二分搜索自动对焦
    Float FocusDistance(Float filmDist);             // 计算特定胶片距离对应的对焦距离
    // 出射光瞳相关方法
    Bounds2f BoundExitPupil(Float pFilmX0, Float pFilmX1) const;  // 计算出射光瞳边界
    void RenderExitPupil(Float sx, Float sy, const char *filename) const;  // 渲染出射光瞳图像
    Point3f SampleExitPupil(const Point2f &pFilm, const Point2f &lensSample,
                            Float *sampleBoundsArea) const;  // 采样出射光瞳位置
    void TestExitPupilBounds() const;  // 测试出射光瞳边界计算是否正确
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
                                       const AnimatedTransform &cam2world,
                                       Film *film, const Medium *medium);

}  // namespace pbrt

#endif  // PBRT_CAMERAS_REALISTIC_H
