
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

#ifndef PBRT_CORE_FILEUTIL_H
#define PBRT_CORE_FILEUTIL_H

// core/fileutil.h*
// 文件工具: 提供平台无关的文件路径操作函数，包括路径解析、目录提取、扩展名检查等
#include "pbrt.h"
#include <string>
#include <cctype>
#include <string.h>

namespace pbrt {

// Platform independent filename-handling functions.
// 平台无关的文件名处理函数
bool IsAbsolutePath(const std::string &filename);       // 判断是否为绝对路径
std::string AbsolutePath(const std::string &filename);  // 获取绝对路径
std::string ResolveFilename(const std::string &filename);  // 解析文件名
std::string DirectoryContaining(const std::string &filename);  // 获取文件所在目录
void SetSearchDirectory(const std::string &dirname);  // 设置搜索目录

// HasExtension: 检查文件名是否具有指定扩展名（大小写不敏感）
inline bool HasExtension(const std::string &value, const std::string &ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(
        ending.rbegin(), ending.rend(), value.rbegin(),
        [](char a, char b) { return std::tolower(a) == std::tolower(b); });
}

}  // namespace pbrt

#endif  // PBRT_CORE_FILEUTIL_H
