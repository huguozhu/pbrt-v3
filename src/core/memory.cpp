
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


// core/memory.cpp*
//
// 本模块提供了对齐内存分配和释放的基础设施。
// AllocAligned 分配缓存行对齐的内存，FreeAligned 释放之。
// 对齐分配对于 SIMD 向量化和避免伪共享（false sharing）
// 具有重要意义，可提升多线程渲染性能。
//
#include "memory.h"

namespace pbrt {

// Memory Allocation Functions

// AllocAligned: 分配指定大小、缓存行对齐的内存
// 在 Windows 上使用 _aligned_malloc，在 POSIX 上使用 posix_memalign
// 对齐到 PBRT_L1_CACHE_LINE_SIZE（通常为 64 字节）
void *AllocAligned(size_t size) {
#if defined(PBRT_HAVE__ALIGNED_MALLOC)
    return _aligned_malloc(size, PBRT_L1_CACHE_LINE_SIZE);
#elif defined(PBRT_HAVE_POSIX_MEMALIGN)
    void *ptr;
    if (posix_memalign(&ptr, PBRT_L1_CACHE_LINE_SIZE, size) != 0) ptr = nullptr;
    return ptr;
#else
    return memalign(PBRT_L1_CACHE_LINE_SIZE, size);
#endif
}

// FreeAligned: 释放由 AllocAligned 分配的对齐内存
// 在 Windows 上使用 _aligned_free，其他平台使用标准 free
void FreeAligned(void *ptr) {
    if (!ptr) return;
#if defined(PBRT_HAVE__ALIGNED_MALLOC)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

}  // namespace pbrt
