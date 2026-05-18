
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


// textures/fbm.cpp*
// 文件功能：分形布朗运动纹理（FBmTexture）的创建函数实现
#include "textures/fbm.h"

namespace pbrt {

// FBmTexture Method Definitions
// 创建浮点型FBM纹理：从参数中读取八度音阶数和粗糙度
FBmTexture<Float> *CreateFBmFloatTexture(const Transform &tex2world,
                                         const TextureParams &tp) {
    // Initialize 3D texture mapping _map_ from _tp_
    std::unique_ptr<TextureMapping3D> map(new IdentityMapping3D(tex2world));
    return new FBmTexture<Float>(std::move(map), tp.FindInt("octaves", 8),
                                 tp.FindFloat("roughness", .5f));
}

FBmTexture<Spectrum> *CreateFBmSpectrumTexture(const Transform &tex2world,
                                               const TextureParams &tp) {
    // Initialize 3D texture mapping _map_ from _tp_
    // 创建光谱型FBM纹理：使用3D恒等映射，默认8个八度，粗糙度0.5
    std::unique_ptr<TextureMapping3D> map(new IdentityMapping3D(tex2world));
    return new FBmTexture<Spectrum>(std::move(map), tp.FindInt("octaves", 8),
                                    tp.FindFloat("roughness", .5f));
}

}  // namespace pbrt
