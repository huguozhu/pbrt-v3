
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

// core/error.cpp*
//
// 此文件实现了 pbrt 的报错和警告功能。
// 提供 Warning() 和 Error() 函数，支持格式化字符串输出，
// 并包含行号/列号信息用于定位场景描述文件中的错误位置。
// 使用互斥锁确保多线程环境下错误消息不会交错输出。
//

#include "error.h"
#include "stringprint.h"
#include "parallel.h"
#include "progressreporter.h"
#include "parser.h"

#include <mutex>

// Error Reporting Includes
#include <stdarg.h>

namespace pbrt {

// Error Reporting Functions

// 可变参数格式化辅助函数：将格式化字符串与 va_list 参数拼接为 std::string。
// 内部使用 vsnprintf 两次调用来确定所需缓冲区大小。
template <typename... Args>
static std::string StringVaprintf(const std::string &fmt, va_list args) {
    // Figure out how much space we need to allocate; add an extra
    // character for '\0'.
    va_list argsCopy;
    va_copy(argsCopy, args);
    size_t size = vsnprintf(nullptr, 0, fmt.c_str(), args) + 1;
    std::string str;
    str.resize(size);
    vsnprintf(&str[0], size, fmt.c_str(), argsCopy);
    str.pop_back();  // remove trailing NUL
    return str;
}

// 内部错误处理函数：构建完整的错误消息字符串并输出到 stderr。
// loc: 场景文件位置（可为 nullptr），用于精确定位错误行。
// format/args: 格式化字符串和参数。
// errorType: 错误类型前缀（如 "Warning" 或 "Error"）。
// 使用静态变量去重，避免相同错误消息重复输出。
static void processError(Loc *loc, const char *format, va_list args,
                         const char *errorType) {
    // Build up an entire formatted error string and print it all at once;
    // this way, if multiple threads are printing messages at once, they
    // don't get jumbled up...
    std::string errorString;

    // Print line and position in input file, if available
    if (loc)
        errorString = StringPrintf("%s:%d:%d: ", loc->filename.c_str(),
                                   loc->line, loc->column);

    errorString += errorType;
    errorString += ": ";
    errorString += StringVaprintf(format, args);

    // Print the error message (but not more than one time).
    static std::string lastError;
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);
    if (errorString != lastError) {
        LOG(INFO) << errorString;
        fprintf(stderr, "%s\n", errorString.c_str());
        lastError = errorString;
    }
}

// 输出警告消息。在安静模式（quiet）下不输出。
// 使用 processError 构建消息，错误类型为 "Warning"。
void Warning(const char *format, ...) {
    if (PbrtOptions.quiet) return;
    va_list args;
    va_start(args, format);
    processError(parserLoc, format, args, "Warning");
    va_end(args);
}

// 输出错误消息。与 Warning 类似，但错误类型为 "Error"。
void Error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(parserLoc, format, args, "Error");
    va_end(args);
}

}  // namespace pbrt
