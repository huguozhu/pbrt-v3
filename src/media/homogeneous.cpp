
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


// media/homogeneous.cpp*
// 文件描述: 均匀参与介质的实现。提供空间均匀介质的透射率计算和散射事件采样。
#include "media/homogeneous.h"
#include "sampler.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"

namespace pbrt {

// HomogeneousMedium Method Definitions
// 计算光线穿过均匀介质的透射率: Tr = exp(-sigma_t * distance)
// 对于均匀介质，透射率有解析解，即比尔-朗伯定律的指数衰减形式。
Spectrum HomogeneousMedium::Tr(const Ray &ray, Sampler &sampler) const {
    ProfilePhase _(Prof::MediumTr);
    return Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxFloat));
}

// 采样均匀介质中的散射事件: 随机选择一个光谱通道，沿光线方向采样距离，
// 如果采样点在介质内则创建介质交互，并计算相应的透射率和采样密度。
Spectrum HomogeneousMedium::Sample(const Ray &ray, Sampler &sampler,
                                   MemoryArena &arena,
                                   MediumInteraction *mi) const {
    ProfilePhase _(Prof::MediumSample);
    // 随机选择一个光谱通道，并沿光线采样传播距离
    // 这是基于比尔-朗伯定律的重要性采样
    int channel = std::min((int)(sampler.Get1D() * Spectrum::nSamples),
                           Spectrum::nSamples - 1);
    Float dist = -std::log(1 - sampler.Get1D()) / sigma_t[channel];
    Float t = std::min(dist / ray.d.Length(), ray.tMax);
    bool sampledMedium = t < ray.tMax;
    if (sampledMedium)
        *mi = MediumInteraction(ray(t), -ray.d, ray.time, this,
                                ARENA_ALLOC(arena, HenyeyGreenstein)(g));

    // 计算透射率和采样密度
    Spectrum Tr = Exp(-sigma_t * std::min(t, MaxFloat) * ray.d.Length());

    // 返回均匀介质的散射权重因子
    Spectrum density = sampledMedium ? (sigma_t * Tr) : Tr;
    Float pdf = 0;
    for (int i = 0; i < Spectrum::nSamples; ++i) pdf += density[i];
    pdf *= 1 / (Float)Spectrum::nSamples;
    if (pdf == 0) {
        CHECK(Tr.IsBlack());
        pdf = 1;
    }
    return sampledMedium ? (Tr * sigma_s / pdf) : (Tr / pdf);
}

}  // namespace pbrt
