
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


// core/shape.cpp*
// 本文件实现了形状(Shape)基类的默认方法。
// Shape是pbrt中所有几何形状的抽象基类，定义了形状的基本接口，
// 包括求交、包围盒计算、面积采样、PDF计算等操作。
// 具体的形状实现(如球体、三角形网格等)继承此类并实现其虚方法。
#include "shape.h"
#include "stats.h"
#include "lowdiscrepancy.h"

namespace pbrt {

// Shape Method Definitions
// Shape析构函数（虚函数），默认实现为空
Shape::~Shape() {}

STAT_COUNTER("Scene/Shapes created", nShapesCreated);
// Shape构造函数：初始化ObjectToWorld和WorldToObject变换矩阵，
// 并记录形状的创建统计信息
Shape::Shape(const Transform *ObjectToWorld, const Transform *WorldToObject,
             bool reverseOrientation)
    : ObjectToWorld(ObjectToWorld),
      WorldToObject(WorldToObject),
      reverseOrientation(reverseOrientation),
      transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {
    ++nShapesCreated;
}

// WorldBound：计算形状的世界空间包围盒，
// 将对象空间的包围盒通过ObjectToWorld变换到世界空间
Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }

// Sample（面积光源重要性采样）：根据参考点ref对形状进行采样，
// 返回采样得到的交点和对应的概率密度。将面积度量转换为立体角度量。
Interaction Shape::Sample(const Interaction &ref, const Point2f &u,
                          Float *pdf) const {
    // 调用形状的面积采样方法，得到采样点intr和面积度量下的PDF
    Interaction intr = Sample(u, pdf);
    // 计算从参考点到采样点的方向向量
    Vector3f wi = intr.p - ref.p;
    // 如果参考点与采样点重合，则PDF为0（立体角为零）
    if (wi.LengthSquared() == 0)
        *pdf = 0;
    else {
        wi = Normalize(wi);
        // 将PDF从面积度量转换为立体角度量
        // 转换公式：solid_angle_pdf = area_pdf * distance^2 / |dot(n, -wi)|
        *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
        // 处理可能的浮点数溢出情况
        if (std::isinf(*pdf)) *pdf = 0.f;
    }
    return intr;
}

// Pdf：计算从参考点ref沿方向wi采样该形状的概率密度值。
// 通过从ref发射射线与形状求交，将PDF从面积度量转换为立体角度量。
Float Shape::Pdf(const Interaction &ref, const Vector3f &wi) const {
    // 从参考点沿方向wi发射射线，与形状求交
    Ray ray = ref.SpawnRay(wi);
    Float tHit;
    SurfaceInteraction isectLight;
    // 忽略用于裁剪形状的alpha纹理（针对"San Miguel"场景的hack，
    // 用于创建不可见面光源）
    if (!Intersect(ray, &tHit, &isectLight, false)) return 0;

    // 将面积度量下的PDF转换为立体角度量
    // 公式：pdf_solid = distance^2 / (|dot(n, -wi)| * area)
    Float pdf = DistanceSquared(ref.p, isectLight.p) /
                (AbsDot(isectLight.n, -wi) * Area());
    if (std::isinf(pdf)) pdf = 0.f;
    return pdf;
}

// SolidAngle：计算形状相对于点p的立体角。
// 使用nSamples个样本进行蒙特卡洛积分估计，
// 通过采样形状表面并检测可见性来计算立体角大小。
Float Shape::SolidAngle(const Point3f &p, int nSamples) const {
    // 构建一个临时的参考交互点，位置为p，法线和方向均为零/默认值
    Interaction ref(p, Normal3f(), Vector3f(), Vector3f(0, 0, 1), 0,
                    MediumInterface{});
    double solidAngle = 0;
    // 使用低差异序列进行蒙特卡洛采样
    for (int i = 0; i < nSamples; ++i) {
        // 生成2D低差异采样点（使用RadicalInverse）
        Point2f u{RadicalInverse(0, i), RadicalInverse(1, i)};
        Float pdf;
        // 对形状进行采样
        Interaction pShape = Sample(ref, u, &pdf);
        // 如果采样有效且点p到采样点之间没有遮挡，则累加贡献
        if (pdf > 0 && !IntersectP(Ray(p, pShape.p - p, .999f))) {
            solidAngle += 1 / pdf;
        }
    }
    return solidAngle / nSamples;
}

}  // namespace pbrt
