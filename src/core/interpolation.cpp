
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

// core/interpolation.cpp*
//
// 本模块提供了 pbrt 中使用的插值函数，包括：
// - Catmull-Rom 样条插值及其相关操作（权重计算、采样、积分、求逆）
// - Fourier 级数插值及采样
// Catmull-Rom 样条是一种 C1 连续的插值样条，具有局部控制特性，
// 在 pbrt 中广泛用于纹理映射、灯光分布等需要平滑插值的场景。
//
#include "interpolation.h"

namespace pbrt {

// Spline Interpolation Definitions

// CatmullRom: 计算 Catmull-Rom 样条在位置 x 处的插值结果
// size - 节点/值数组长度
// nodes - 节点位置数组（单调递增）
// values - 节点处的函数值数组
// x - 需要插值的自变量位置
// 返回值：插值结果，若 x 超出节点范围则返回 0
Float CatmullRom(int size, const Float *nodes, const Float *values, Float x) {
    // 检查 x 是否在节点范围内，超出则返回 0
    if (!(x >= nodes[0] && x <= nodes[size - 1])) return 0;
    // 二分查找包含 x 的区间
    int idx = FindInterval(size, [&](int i) { return nodes[i] <= x; });
    Float x0 = nodes[idx], x1 = nodes[idx + 1];
    Float f0 = values[idx], f1 = values[idx + 1];
    Float width = x1 - x0;
    // 使用有限差分近似导数（Catmull-Rom 样条的标准方法）
    Float d0, d1;
    if (idx > 0)
        d0 = width * (f1 - values[idx - 1]) / (x1 - nodes[idx - 1]);
    else
        d0 = f1 - f0;

    if (idx + 2 < size)
        d1 = width * (values[idx + 2] - f0) / (nodes[idx + 2] - x0);
    else
        d1 = f1 - f0;

    // 计算局部参数 t 并求三次 Hermite 插值
    Float t = (x - x0) / (x1 - x0), t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * f0 + (-2 * t3 + 3 * t2) * f1 +
           (t3 - 2 * t2 + t) * d0 + (t3 - t2) * d1;
}

// CatmullRomWeights: 计算 Catmull-Rom 样条插值的权重和偏移量
// 与 CatmullRom 不同，此函数返回四个相邻节点的权重，
// 便于调用者自行组合多个样条结果（如 2D/3D 插值）
// size - 节点数组长度
// nodes - 节点位置数组
// x - 插值位置
// offset - 输出参数，返回四个权重对应的起始节点索引减 1
// weights - 输出参数，长度为 4 的权重数组
// 返回值：若 x 在范围内则返回 true，否则返回 false
bool CatmullRomWeights(int size, const Float *nodes, Float x, int *offset,
                       Float *weights) {
    // Return _false_ if _x_ is out of bounds
    if (!(x >= nodes[0] && x <= nodes[size - 1])) return false;

    // Search for the interval _idx_ containing _x_
    int idx = FindInterval(size, [&](int i) { return nodes[i] <= x; });
    *offset = idx - 1;
    Float x0 = nodes[idx], x1 = nodes[idx + 1];

    // Compute the $t$ parameter and powers
    Float t = (x - x0) / (x1 - x0), t2 = t * t, t3 = t2 * t;

    // Compute initial node weights $w_1$ and $w_2$
    weights[1] = 2 * t3 - 3 * t2 + 1;
    weights[2] = -2 * t3 + 3 * t2;

    // Compute first node weight $w_0$
    if (idx > 0) {
        Float w0 = (t3 - 2 * t2 + t) * (x1 - x0) / (x1 - nodes[idx - 1]);
        weights[0] = -w0;
        weights[2] += w0;
    } else {
        Float w0 = t3 - 2 * t2 + t;
        weights[0] = 0;
        weights[1] -= w0;
        weights[2] += w0;
    }

    // Compute last node weight $w_3$
    if (idx + 2 < size) {
        Float w3 = (t3 - t2) * (x1 - x0) / (nodes[idx + 2] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        Float w3 = t3 - t2;
        weights[1] -= w3;
        weights[2] += w3;
        weights[3] = 0;
    }
    return true;
}

// SampleCatmullRom: 从 Catmull-Rom 样条表示的概率分布中采样
// 通过对 CDF 求逆实现重要性采样，使用 Newton-Bisection 方法求解
// n - 样本点数量
// x - 节点位置数组
// f - 节点处的函数值（未归一化的 PDF）
// F - 累积分布函数（CDF）值
// u - [0,1) 均匀随机数
// fval - 输出参数，返回采样点处的 PDF 值
// pdf - 输出参数，返回归一化后的概率密度
// 返回值：采样位置
Float SampleCatmullRom(int n, const Float *x, const Float *f, const Float *F,
                       Float u, Float *fval, Float *pdf) {
    // Map _u_ to a spline interval by inverting _F_
    // 通过求逆 CDF 将 u 映射到样条区间
    u *= F[n - 1];
    int i = FindInterval(n, [&](int i) { return F[i] <= u; });

    // Look up $x_i$ and function values of spline segment _i_
    // 查找样条段 i 的节点位置和函数值
    Float x0 = x[i], x1 = x[i + 1];
    Float f0 = f[i], f1 = f[i + 1];
    Float width = x1 - x0;

    // Approximate derivatives using finite differences
    // 使用有限差分近似导数
    Float d0, d1;
    if (i > 0)
        d0 = width * (f1 - f[i - 1]) / (x1 - x[i - 1]);
    else
        d0 = f1 - f0;
    if (i + 2 < n)
        d1 = width * (f[i + 2] - f0) / (x[i + 2] - x0);
    else
        d1 = f1 - f0;

    // Re-scale _u_ for continous spline sampling step
    // 重新缩放 u 以用于样条段内的连续采样
    u = (u - F[i]) / width;

    // Invert definite integral over spline segment and return solution
    // 求解样条段上的定积分逆变换

    // Set initial guess for $t$ by importance sampling a linear interpolant
    // 使用线性插值的重要性采样设置初始猜测 t
    Float t;
    if (f0 != f1)
        t = (f0 - std::sqrt(std::max((Float)0, f0 * f0 + 2 * u * (f1 - f0)))) /
            (f0 - f1);
    else
        t = u / f0;
    Float a = 0, b = 1, Fhat, fhat;
    while (true) {
        // Fall back to a bisection step when _t_ is out of bounds
        // 当 t 超出范围时回退到二分法
        if (!(t > a && t < b)) t = 0.5f * (a + b);

        // Evaluate target function and its derivative in Horner form
        // 使用 Horner 形式计算目标函数及其导数
        Fhat = t * (f0 +
                    t * (.5f * d0 +
                         t * ((1.f / 3.f) * (-2 * d0 - d1) + f1 - f0 +
                              t * (.25f * (d0 + d1) + .5f * (f0 - f1)))));
        fhat = f0 +
               t * (d0 +
                    t * (-2 * d0 - d1 + 3 * (f1 - f0) +
                         t * (d0 + d1 + 2 * (f0 - f1))));

        // Stop the iteration if converged
        // 收敛判定
        if (std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) break;

        // Update bisection bounds using updated _t_
        // 更新二分法边界
        if (Fhat - u < 0)
            a = t;
        else
            b = t;

        // Perform a Newton step
        // 执行 Newton 迭代步
        t -= (Fhat - u) / fhat;
    }

    // Return the sample position and function value
    // 返回采样位置和函数值
    if (fval) *fval = fhat;
    if (pdf) *pdf = fhat / F[n - 1];
    return x0 + width * t;
}

// SampleCatmullRom2D: 从二维 Catmull-Rom 样条分布中采样
// 利用已计算的一维 CDF，对给定 alpha（第一维参数）进行条件采样
// size1, size2 - 两个维度的样本点数量
// nodes1, nodes2 - 两个维度的节点位置数组
// values - 二维函数值表（大小为 size1 * size2）
// cdf - 第二维的 CDF 值表（大小为 size1 * size2）
// alpha - 第一维的参数值
// u - [0,1) 均匀随机数
// fval, pdf - 输出参数，PDF 值和概率密度
// 返回值：第二维的采样位置
Float SampleCatmullRom2D(int size1, int size2, const Float *nodes1,
                         const Float *nodes2, const Float *values,
                         const Float *cdf, Float alpha, Float u, Float *fval,
                         Float *pdf) {
    // Determine offset and coefficients for the _alpha_ parameter
    // 确定 alpha 参数对应的 Catmull-Rom 权重和偏移量
    int offset;
    Float weights[4];
    if (!CatmullRomWeights(size1, nodes1, alpha, &offset, weights)) return 0;

    // Define a lambda function to interpolate table entries
    // 定义 lambda 函数，用于在二维表中沿第二维插值
    auto interpolate = [&](const Float *array, int idx) {
        Float value = 0;
        for (int i = 0; i < 4; ++i)
            if (weights[i] != 0)
                value += array[(offset + i) * size2 + idx] * weights[i];
        return value;
    };

    // Map _u_ to a spline interval by inverting the interpolated _cdf_
    // 通过求逆插值后的 CDF 将 u 映射到样条区间
    Float maximum = interpolate(cdf, size2 - 1);
    u *= maximum;
    int idx =
        FindInterval(size2, [&](int i) { return interpolate(cdf, i) <= u; });

    // Look up node positions and interpolated function values
    // 查找节点位置和插值后的函数值
    Float f0 = interpolate(values, idx), f1 = interpolate(values, idx + 1);
    Float x0 = nodes2[idx], x1 = nodes2[idx + 1];
    Float width = x1 - x0;
    Float d0, d1;

    // Re-scale _u_ using the interpolated _cdf_
    // 使用插值后的 CDF 重新缩放 u
    u = (u - interpolate(cdf, idx)) / width;

    // Approximate derivatives using finite differences of the interpolant
    // 使用插值函数的有限差分近似导数
    if (idx > 0)
        d0 = width * (f1 - interpolate(values, idx - 1)) /
             (x1 - nodes2[idx - 1]);
    else
        d0 = f1 - f0;
    if (idx + 2 < size2)
        d1 = width * (interpolate(values, idx + 2) - f0) /
             (nodes2[idx + 2] - x0);
    else
        d1 = f1 - f0;

    // Invert definite integral over spline segment and return solution
    // 求解样条段上的定积分逆变换（Newton-Bisection）

    // Set initial guess for $t$ by importance sampling a linear interpolant
    Float t;
    if (f0 != f1)
        t = (f0 - std::sqrt(std::max((Float)0, f0 * f0 + 2 * u * (f1 - f0)))) /
            (f0 - f1);
    else
        t = u / f0;
    Float a = 0, b = 1, Fhat, fhat;
    while (true) {
        // Fall back to a bisection step when _t_ is out of bounds
        if (!(t >= a && t <= b)) t = 0.5f * (a + b);

        // Evaluate target function and its derivative in Horner form
        Fhat = t * (f0 +
                    t * (.5f * d0 +
                         t * ((1.f / 3.f) * (-2 * d0 - d1) + f1 - f0 +
                              t * (.25f * (d0 + d1) + .5f * (f0 - f1)))));
        fhat = f0 +
               t * (d0 +
                    t * (-2 * d0 - d1 + 3 * (f1 - f0) +
                         t * (d0 + d1 + 2 * (f0 - f1))));

        // Stop the iteration if converged
        if (std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) break;

        // Update bisection bounds using updated _t_
        if (Fhat - u < 0)
            a = t;
        else
            b = t;

        // Perform a Newton step
        t -= (Fhat - u) / fhat;
    }

    // Return the sample position and function value
    if (fval) *fval = fhat;
    if (pdf) *pdf = fhat / maximum;
    return x0 + width * t;
}

// IntegrateCatmullRom: 对 Catmull-Rom 样条进行数值积分并构建 CDF
// n - 样本点数量
// x - 节点位置数组
// values - 节点处的函数值数组
// cdf - 输出参数，累积分布函数数组（长度为 n）
// 返回值：积分总和（归一化常数）
Float IntegrateCatmullRom(int n, const Float *x, const Float *values,
                          Float *cdf) {
    Float sum = 0;
    cdf[0] = 0;
    for (int i = 0; i < n - 1; ++i) {
        // Look up $x_i$ and function values of spline segment _i_
        // 查找样条段 i 的节点位置和函数值
        Float x0 = x[i], x1 = x[i + 1];
        Float f0 = values[i], f1 = values[i + 1];
        Float width = x1 - x0;

        // Approximate derivatives using finite differences
        // 使用有限差分近似导数
        Float d0, d1;
        if (i > 0)
            d0 = width * (f1 - values[i - 1]) / (x1 - x[i - 1]);
        else
            d0 = f1 - f0;
        if (i + 2 < n)
            d1 = width * (values[i + 2] - f0) / (x[i + 2] - x0);
        else
            d1 = f1 - f0;

        // Keep a running sum and build a cumulative distribution function
        // 计算样条段上的解析积分并构建 CDF
        sum += ((d0 - d1) * (1.f / 12.f) + (f0 + f1) * .5f) * width;
        cdf[i + 1] = sum;
    }
    return sum;
}

// InvertCatmullRom: 求 Catmull-Rom 样条的逆函数
// 给定函数值 u，返回对应的自变量 x（假设样条单调）
// 使用 Newton-Bisection 混合方法求解
// n - 样本点数量
// x - 节点位置数组
// values - 节点处的函数值数组
// u - 需要求逆的函数值
// 返回值：对应的自变量位置
Float InvertCatmullRom(int n, const Float *x, const Float *values, Float u) {
    // Stop when _u_ is out of bounds
    if (!(u > values[0]))
        return x[0];
    else if (!(u < values[n - 1]))
        return x[n - 1];

    // Map _u_ to a spline interval by inverting _values_
    int i = FindInterval(n, [&](int i) { return values[i] <= u; });

    // Look up $x_i$ and function values of spline segment _i_
    Float x0 = x[i], x1 = x[i + 1];
    Float f0 = values[i], f1 = values[i + 1];
    Float width = x1 - x0;

    // Approximate derivatives using finite differences
    Float d0, d1;
    if (i > 0)
        d0 = width * (f1 - values[i - 1]) / (x1 - x[i - 1]);
    else
        d0 = f1 - f0;
    if (i + 2 < n)
        d1 = width * (values[i + 2] - f0) / (x[i + 2] - x0);
    else
        d1 = f1 - f0;

    // Invert the spline interpolant using Newton-Bisection
    Float a = 0, b = 1, t = .5f;
    Float Fhat, fhat;
    while (true) {
        // Fall back to a bisection step when _t_ is out of bounds
        if (!(t > a && t < b)) t = 0.5f * (a + b);

        // Compute powers of _t_
        Float t2 = t * t, t3 = t2 * t;

        // Set _Fhat_ using Equation (8.27)
        Fhat = (2 * t3 - 3 * t2 + 1) * f0 + (-2 * t3 + 3 * t2) * f1 +
               (t3 - 2 * t2 + t) * d0 + (t3 - t2) * d1;

        // Set _fhat_ using Equation (not present)
        fhat = (6 * t2 - 6 * t) * f0 + (-6 * t2 + 6 * t) * f1 +
               (3 * t2 - 4 * t + 1) * d0 + (3 * t2 - 2 * t) * d1;

        // Stop the iteration if converged
        if (std::abs(Fhat - u) < 1e-6f || b - a < 1e-6f) break;

        // Update bisection bounds using updated _t_
        if (Fhat - u < 0)
            a = t;
        else
            b = t;

        // Perform a Newton step
        t -= (Fhat - u) / fhat;
    }
    return x0 + t * width;
}

// Fourier Interpolation Definitions

// Fourier: 使用余弦级数计算 Fourier 插值
// a - Fourier 级数系数数组
// m - 级数项数
// cosPhi - cos(phi) 的值
// 返回值：Fourier 级数在角度 phi 处的求值结果
// 使用 Chebyshev 递推关系高效计算余弦项
Float Fourier(const Float *a, int m, double cosPhi) {
    double value = 0.0;
    // Initialize cosine iterates
    // 初始化余弦递推：cos(0*phi)=1, cos(1*phi)=cosPhi
    double cosKMinusOnePhi = cosPhi;
    double cosKPhi = 1;
    for (int k = 0; k < m; ++k) {
        // Add the current summand and update the cosine iterates
        // 使用 Chebyshev 递推公式：cos((k+1)phi) = 2*cosPhi*cos(k*phi) - cos((k-1)phi)
        value += a[k] * cosKPhi;
        double cosKPlusOnePhi = 2 * cosPhi * cosKPhi - cosKMinusOnePhi;
        cosKMinusOnePhi = cosKPhi;
        cosKPhi = cosKPlusOnePhi;
    }
    return value;
}

// SampleFourier: 从 Fourier 级数表示的各向异性相位函数中采样方向
// 通过对累积分布函数求逆来采样角度 phi
// ak - Fourier 级数的余弦系数
// recip - 系数 k 的倒数数组（recip[k] = 1/k）
// m - 级数项数
// u - [0,1) 均匀随机数
// pdf - 输出参数，返回采样方向的概率密度
// phiPtr - 输出参数，返回采样的方位角
// 返回值：采样点处的相位函数值
Float SampleFourier(const Float *ak, const Float *recip, int m, Float u,
                    Float *pdf, Float *phiPtr) {
    // Pick a side and declare bisection variables
    // 利用对称性处理 [0, 2pi) 区间，映射到 [0, pi)
    bool flip = (u >= 0.5);
    if (flip)
        u = 1 - 2 * (u - .5f);
    else
        u *= 2;
    double a = 0, b = Pi, phi = 0.5 * Pi;
    double F, f;
    while (true) {
        // Evaluate $F(\phi)$ and its derivative $f(\phi)$
        // 计算 Fourier 级数在 phi 处的 CDF 值 F 和 PDF 值 f

        // Initialize sine and cosine iterates
        // 初始化正弦和余弦递推
        double cosPhi = std::cos(phi);
        double sinPhi = std::sqrt(std::max(0., 1 - cosPhi * cosPhi));
        double cosPhiPrev = cosPhi, cosPhiCur = 1;
        double sinPhiPrev = -sinPhi, sinPhiCur = 0;

        // Initialize _F_ and _f_ with the first series term
        // 用级数第一项初始化 F 和 f
        F = ak[0] * phi;
        f = ak[0];
        for (int k = 1; k < m; ++k) {
            // Compute next sine and cosine iterates
            // 使用 Chebyshev 递推计算下一项的正弦和余弦
            double sinPhiNext = 2 * cosPhi * sinPhiCur - sinPhiPrev;
            double cosPhiNext = 2 * cosPhi * cosPhiCur - cosPhiPrev;
            sinPhiPrev = sinPhiCur;
            sinPhiCur = sinPhiNext;
            cosPhiPrev = cosPhiCur;
            cosPhiCur = cosPhiNext;

            // Add the next series term to _F_ and _f_
            // 累加下一项到 F（CDF）和 f（PDF）
            F += ak[k] * recip[k] * sinPhiNext;
            f += ak[k] * cosPhiNext;
        }
        F -= u * ak[0] * Pi;

        // Update bisection bounds using updated $\phi$
        // 根据 F 的符号更新二分法边界
        if (F > 0)
            b = phi;
        else
            a = phi;

        // Stop the Fourier bisection iteration if converged
        // 收敛判定
        if (std::abs(F) < 1e-6f || b - a < 1e-6f) break;

        // Perform a Newton step given $f(\phi)$ and $F(\phi)$
        // 执行 Newton 迭代步
        phi -= F / f;

        // Fall back to a bisection step when $\phi$ is out of bounds
        // 当 phi 超出范围时回退到二分法
        if (!(phi > a && phi < b)) phi = 0.5f * (a + b);
    }
    // Potentially flip $\phi$ and return the result
    // 如果之前翻转了，将 phi 映射回 [0, 2pi)
    if (flip) phi = 2 * Pi - phi;
    *pdf = (Float)(Inv2Pi * f / ak[0]);
    *phiPtr = (Float)phi;
    return f;
}

}  // namespace pbrt
