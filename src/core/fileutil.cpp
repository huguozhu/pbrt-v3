
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


// core/fileutil.cpp*
//
// 此文件实现了文件路径相关的工具函数。
// 包括判断路径是否为绝对路径、获取绝对路径、解析文件名、
// 获取文件所在目录以及设置搜索目录等功能。
// 在 Windows 和 POSIX 系统下分别使用不同的实现。
//

#include "fileutil.h"
#include <cstdlib>
#include <climits>
#ifndef PBRT_IS_WINDOWS
#include <libgen.h>
#endif

namespace pbrt {

static std::string searchDirectory;

#ifdef PBRT_IS_WINDOWS
// 判断文件名是否为绝对路径（Windows 版本）。
// 以 '\'、'/' 开头或包含 ':'（如 "C:\"）均视为绝对路径。
bool IsAbsolutePath(const std::string &filename) {
    if (filename.empty()) return false;
    return (filename[0] == '\\' || filename[0] == '/' ||
            filename.find(':') != std::string::npos);
}

// 将相对路径转换为绝对路径（Windows 版本）。
// 使用 _fullpath 系统调用，失败时返回原文件名。
std::string AbsolutePath(const std::string &filename) {
    char full[_MAX_PATH];
    if (_fullpath(full, filename.c_str(), _MAX_PATH))
        return std::string(full);
    else
        return filename;
}

// 解析并补全文件名：如果是相对路径则在搜索目录前添加前缀（Windows 版本）。
std::string ResolveFilename(const std::string &filename) {
    if (searchDirectory.empty() || filename.empty())
        return filename;
    else if (IsAbsolutePath(filename))
        return filename;

    char searchDirectoryEnd = searchDirectory[searchDirectory.size() - 1];
    if (searchDirectoryEnd == '\\' || searchDirectoryEnd == '/')
        return searchDirectory + filename;
    else
        return searchDirectory + "\\" + filename;
}

// 返回文件所在目录路径（Windows 版本）。
// 使用 _splitpath_s 和 _makepath_s 分割和重组路径组件。
std::string DirectoryContaining(const std::string &filename) {
    // This code isn't tested but I believe it should work. Might need to add
    // some const_casts to make it compile though.
    char drive[_MAX_DRIVE];
    char dir[_MAX_DIR];
    char ext[_MAX_EXT];

    errno_t err = _splitpath_s(filename.c_str(), drive, _MAX_DRIVE, dir,
                               _MAX_DIR, nullptr, 0, ext, _MAX_EXT);
    if (err == 0) {
        char fullDir[_MAX_PATH];
        err = _makepath_s(fullDir, _MAX_PATH, drive, dir, nullptr, nullptr);
        if (err == 0) return std::string(fullDir);
    }
    return filename;
}

#else

// 判断文件名是否为绝对路径（POSIX 版本）。以 '/' 开头为绝对路径。
bool IsAbsolutePath(const std::string &filename) {
    return (filename.size() > 0) && filename[0] == '/';
}

// 将相对路径转换为绝对路径（POSIX 版本）。使用 realpath 系统调用。
std::string AbsolutePath(const std::string &filename) {
    char full[PATH_MAX];
    if (realpath(filename.c_str(), full))
        return std::string(full);
    else
        return filename;
}

// 解析并补全文件名：如果是相对路径则在搜索目录前添加前缀（POSIX 版本）。
std::string ResolveFilename(const std::string &filename) {
    if (searchDirectory.empty() || filename.empty())
        return filename;
    else if (IsAbsolutePath(filename))
        return filename;
    else if (searchDirectory[searchDirectory.size() - 1] == '/')
        return searchDirectory + filename;
    else
        return searchDirectory + "/" + filename;
}

// 返回文件所在目录路径（POSIX 版本）。使用 dirname() 系统函数。
std::string DirectoryContaining(const std::string &filename) {
    char *t = strdup(filename.c_str());
    std::string result = dirname(t);
    free(t);
    return result;
}

#endif

// 设置全局搜索目录，用于 ResolveFilename 补全相对路径。
void SetSearchDirectory(const std::string &dirname) {
    searchDirectory = dirname;
}

}  // namespace pbrt
