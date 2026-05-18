
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


// core/texture.cpp*
// 本文件实现了纹理映射和程序化纹理生成的基础设施。
// 包含2D/3D纹理坐标映射方法（UV/球形/圆柱/平面映射），
// 以及Perlin噪声、FBm分形布朗运动、湍流(Turbulence)
// 和Lanczos滤波器等程序化纹理生成函数。
#include "texture.h"
#include "shape.h"

namespace pbrt {

// Texture Inline Functions
// SmoothStep：平滑阶跃函数，在[min, max]区间内进行三次Hermite插值。
// 返回值在0到1之间平滑过渡，用于噪声函数中的平滑插值。
inline Float SmoothStep(Float min, Float max, Float value) {
    Float v = Clamp((value - min) / (max - min), 0, 1);
    return v * v * (-2 * v + 3);
}

// Texture Forward Declarations
// Grad：使用Permutation Table计算Perlin噪声在整数坐标(x,y,z)处的梯度值
inline Float Grad(int x, int y, int z, Float dx, Float dy, Float dz);
// NoiseWeight：Perlin噪声的插值权重函数，使用6t^5-15t^4+10t^3的平滑曲线
inline Float NoiseWeight(Float t);

// Perlin Noise Data
// Perlin噪声的排列表，大小为256的排列数组，通过重复两倍实现无边界环绕查找
static PBRT_CONSTEXPR int NoisePermSize = 256;
static int NoisePerm[2 * NoisePermSize] = {
    151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140,
    36, 103, 30, 69, 142,
    // Remainder of the noise permutation table
    8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62,
    94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174,
    20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77,
    146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55,
    46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76,
    132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100,
    109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147,
    118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28,
    42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101,
    155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
    178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12,
    191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31,
    181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
    138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66,
    215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194,
    233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6,
    148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32,
    57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74,
    165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60,
    211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25,
    63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135,
    130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226,
    250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59,
    227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2,
    44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19,
    98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251,
    34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249,
    14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115,
    121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72,
    243, 141, 128, 195, 78, 66, 215, 61, 156, 180};

// Texture Method Definitions
// 2D纹理映射基类虚析构函数
TextureMapping2D::~TextureMapping2D() { }
// 3D纹理映射基类虚析构函数
TextureMapping3D::~TextureMapping3D() { }

// UVMapping2D构造函数：设置UV缩放和平移参数
UVMapping2D::UVMapping2D(Float su, Float sv, Float du, Float dv)
    : su(su), sv(sv), du(du), dv(dv) {}
// UVMapping2D::Map：执行2D UV纹理映射，返回纹理坐标并计算纹理坐标的屏幕空间微分
Point2f UVMapping2D::Map(const SurfaceInteraction &si, Vector2f *dstdx,
                         Vector2f *dstdy) const {
    // Compute texture differentials for 2D identity mapping
    *dstdx = Vector2f(su * si.dudx, sv * si.dvdx);
    *dstdy = Vector2f(su * si.dudy, sv * si.dvdy);
    return Point2f(su * si.uv[0] + du, sv * si.uv[1] + dv);
}

// SphericalMapping2D::Map：执行球体纹理映射，基于表面位置计算球面纹理坐标，
// 并通过有限差分计算纹理坐标的屏幕空间微分
Point2f SphericalMapping2D::Map(const SurfaceInteraction &si, Vector2f *dstdx,
                                Vector2f *dstdy) const {
    Point2f st = sphere(si.p);
    // 使用有限差分法计算球体映射的纹理坐标微分
    const Float delta = .1f;
    Point2f stDeltaX = sphere(si.p + delta * si.dpdx);
    *dstdx = (stDeltaX - st) / delta;
    Point2f stDeltaY = sphere(si.p + delta * si.dpdy);
    *dstdy = (stDeltaY - st) / delta;

    // 处理球体映射中v方向在phi=0/2Pi处的纹理坐标微分不连续性
    if ((*dstdx)[1] > .5)
        (*dstdx)[1] = 1 - (*dstdx)[1];
    else if ((*dstdx)[1] < -.5f)
        (*dstdx)[1] = -((*dstdx)[1] + 1);
    if ((*dstdy)[1] > .5)
        (*dstdy)[1] = 1 - (*dstdy)[1];
    else if ((*dstdy)[1] < -.5f)
        (*dstdy)[1] = -((*dstdy)[1] + 1);
    return st;
}

// SphericalMapping2D::sphere：将3D点p映射到球面坐标(theta, phi)，
// 返回归一化到[0,1]范围的纹理坐标
Point2f SphericalMapping2D::sphere(const Point3f &p) const {
    Vector3f vec = Normalize(WorldToTexture(p) - Point3f(0, 0, 0));
    Float theta = SphericalTheta(vec), phi = SphericalPhi(vec);
    return Point2f(theta * InvPi, phi * Inv2Pi);
}

// CylindricalMapping2D::Map：执行圆柱纹理映射，基于表面位置计算柱面纹理坐标，
// 并通过有限差分计算纹理坐标的屏幕空间微分
Point2f CylindricalMapping2D::Map(const SurfaceInteraction &si, Vector2f *dstdx,
                                  Vector2f *dstdy) const {
    Point2f st = cylinder(si.p);
    // Compute texture coordinate differentials for cylinder $(u,v)$ mapping
    const Float delta = .01f;
    Point2f stDeltaX = cylinder(si.p + delta * si.dpdx);
    *dstdx = (stDeltaX - st) / delta;
    if ((*dstdx)[1] > .5)
        (*dstdx)[1] = 1.f - (*dstdx)[1];
    else if ((*dstdx)[1] < -.5f)
        (*dstdx)[1] = -((*dstdx)[1] + 1);
    Point2f stDeltaY = cylinder(si.p + delta * si.dpdy);
    *dstdy = (stDeltaY - st) / delta;
    if ((*dstdy)[1] > .5)
        (*dstdy)[1] = 1.f - (*dstdy)[1];
    else if ((*dstdy)[1] < -.5f)
        (*dstdy)[1] = -((*dstdy)[1] + 1);
    return st;
}

// PlanarMapping2D::Map：执行平面纹理映射，使用两个正交向量(vs, vt)
// 将3D点投影到平面上得到2D纹理坐标，并计算屏幕空间微分
Point2f PlanarMapping2D::Map(const SurfaceInteraction &si, Vector2f *dstdx,
                             Vector2f *dstdy) const {
    Vector3f vec(si.p);
    *dstdx = Vector2f(Dot(si.dpdx, vs), Dot(si.dpdx, vt));
    *dstdy = Vector2f(Dot(si.dpdy, vs), Dot(si.dpdy, vt));
    return Point2f(ds + Dot(vec, vs), dt + Dot(vec, vt));
}

// IdentityMapping3D::Map：执行3D恒等纹理映射，
// 将表面交互点和其微分从世界空间变换到纹理空间
Point3f IdentityMapping3D::Map(const SurfaceInteraction &si, Vector3f *dpdx,
                               Vector3f *dpdy) const {
    *dpdx = WorldToTexture(si.dpdx);
    *dpdy = WorldToTexture(si.dpdy);
    return WorldToTexture(si.p);
}

// Noise：计算3D Perlin噪声值。
// 通过对整数格点上的随机梯度进行三线性插值生成连续的噪声函数。
Float Noise(Float x, Float y, Float z) {
    // 计算噪声格点坐标（整数部分）和格内偏移（小数部分）
    int ix = std::floor(x), iy = std::floor(y), iz = std::floor(z);
    Float dx = x - ix, dy = y - iy, dz = z - iz;

    // 使用排列表计算8个格点角的梯度权重
    ix &= NoisePermSize - 1;
    iy &= NoisePermSize - 1;
    iz &= NoisePermSize - 1;
    Float w000 = Grad(ix, iy, iz, dx, dy, dz);
    Float w100 = Grad(ix + 1, iy, iz, dx - 1, dy, dz);
    Float w010 = Grad(ix, iy + 1, iz, dx, dy - 1, dz);
    Float w110 = Grad(ix + 1, iy + 1, iz, dx - 1, dy - 1, dz);
    Float w001 = Grad(ix, iy, iz + 1, dx, dy, dz - 1);
    Float w101 = Grad(ix + 1, iy, iz + 1, dx - 1, dy, dz - 1);
    Float w011 = Grad(ix, iy + 1, iz + 1, dx, dy - 1, dz - 1);
    Float w111 = Grad(ix + 1, iy + 1, iz + 1, dx - 1, dy - 1, dz - 1);

    // 对8个格点角的梯度权重执行三线性插值，
    // 使用NoiseWeight生成平滑的插值系数
    Float wx = NoiseWeight(dx), wy = NoiseWeight(dy), wz = NoiseWeight(dz);
    Float x00 = Lerp(wx, w000, w100);
    Float x10 = Lerp(wx, w010, w110);
    Float x01 = Lerp(wx, w001, w101);
    Float x11 = Lerp(wx, w011, w111);
    Float y0 = Lerp(wy, x00, x10);
    Float y1 = Lerp(wy, x01, x11);
    return Lerp(wz, y0, y1);
}

// Noise：基于Point3f重载的Perlin噪声函数
Float Noise(const Point3f &p) { return Noise(p.x, p.y, p.z); }
// Grad：使用Perlin排列表计算格点(x,y,z)处对应于偏移量(dx,dy,dz)的梯度值。
// 通过排列表查找到16种预定义梯度方向之一，计算梯度与偏移向量的点积。
inline Float Grad(int x, int y, int z, Float dx, Float dy, Float dz) {
    int h = NoisePerm[NoisePerm[NoisePerm[x] + y] + z];
    h &= 15;
    Float u = h < 8 || h == 12 || h == 13 ? dx : dy;
    Float v = h < 4 || h == 12 || h == 13 ? dy : dz;
    return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
}

// NoiseWeight：Perlin噪声的6次多项式插值权重函数，
// 公式为 6t^5 - 15t^4 + 10t^3，其一阶和二阶导数在t=0和t=1处均为零，
// 保证了噪声函数在格点处的平滑连续性。
inline Float NoiseWeight(Float t) {
    Float t3 = t * t * t;
    Float t4 = t3 * t;
    return 6 * t4 * t - 15 * t4 + 10 * t3;
}

// FBm：分形布朗运动(Fractional Brownian Motion)纹理函数。
// 叠加多个频率和幅度递减的噪声八度(octave)来模拟自然的自相似纹理。
// 使用屏幕空间微分自动计算合适的八度数以实现抗锯齿。
Float FBm(const Point3f &p, const Vector3f &dpdx, const Vector3f &dpdy,
          Float omega, int maxOctaves) {
    // Compute number of octaves for antialiased FBm
    Float len2 = std::max(dpdx.LengthSquared(), dpdy.LengthSquared());
    Float n = Clamp(-1 - .5f * Log2(len2), 0, maxOctaves);
    int nInt = std::floor(n);

    // 累加多个噪声八度的贡献，频率逐步倍增(lambda *= 1.99)，
    // 幅度按omega因子递减
    Float sum = 0, lambda = 1, o = 1;
    for (int i = 0; i < nInt; ++i) {
        sum += o * Noise(lambda * p);
        lambda *= 1.99f;
        o *= omega;
    }
    // 处理小数八度部分：使用SmoothStep平滑过渡以避免视觉突变
    Float nPartial = n - nInt;
    sum += o * SmoothStep(.3f, .7f, nPartial) * Noise(lambda * p);
    return sum;
}

// Turbulence：湍流噪声函数，类似于FBm但对每个噪声八度取绝对值，
// 产生更尖锐、更不规则的纹理效果，适合模拟火焰、云层等自然现象。
Float Turbulence(const Point3f &p, const Vector3f &dpdx, const Vector3f &dpdy,
                 Float omega, int maxOctaves) {
    // Compute number of octaves for antialiased FBm
    Float len2 = std::max(dpdx.LengthSquared(), dpdy.LengthSquared());
    Float n = Clamp(-1 - .5f * Log2(len2), 0, maxOctaves);
    int nInt = std::floor(n);

    // 累加多个噪声八度的绝对值（产生尖锐的湍流效果）
    Float sum = 0, lambda = 1, o = 1;
    for (int i = 0; i < nInt; ++i) {
        sum += o * std::abs(Noise(lambda * p));
        lambda *= 1.99f;
        o *= omega;
    }

    // 处理被截断的高频八度：对小数部分使用平滑过渡，
    // 对完全被截断的八度使用其预期贡献的期望值0.2
    // Account for contributions of clamped octaves in turbulence
    Float nPartial = n - nInt;
    sum += o * Lerp(SmoothStep(.3f, .7f, nPartial), 0.2,
                    std::abs(Noise(lambda * p)));
    for (int i = nInt; i < maxOctaves; ++i) {
        sum += o * 0.2f;
        o *= omega;
    }
    return sum;
}

// Texture Function Definitions
// Lanczos：Lanczos sinc重采样滤波器。
// tau控制滤波器的宽度（通常为2或3），
// 在图像重采样和纹理滤波中用于平衡振铃伪影和细节保留。
Float Lanczos(Float x, Float tau) {
    x = std::abs(x);
    if (x < 1e-5f) return 1;
    if (x > 1.f) return 0;
    x *= Pi;
    Float s = std::sin(x * tau) / (x * tau);
    Float lanczos = std::sin(x) / x;
    return s * lanczos;
}

}  // namespace pbrt
