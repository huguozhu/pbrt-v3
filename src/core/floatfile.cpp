
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

// core/floatfile.cpp*
//
// 本模块提供了从文本文件中读取浮点数数组的功能。
// ReadFloatFile 函数解析包含空格分隔的浮点数字的文本文件，
// 支持注释行（以#开头）、换行符分隔等常见文本格式。
//
#include "floatfile.h"
#include <ctype.h>
#include <stdlib.h>

namespace pbrt {

// ReadFloatFile: 从文本文件中读取浮点数列表
// 文件名由 filename 指定，解析结果以 Float 类型存入 values 向量中。
// 文件格式：以空白字符分隔的浮点数，支持单行注释（# 开头）。
// 返回值表示是否成功读取。
bool ReadFloatFile(const char *filename, std::vector<Float> *values) {
    // 打开文件，若失败则报告错误并返回 false
    FILE *f = fopen(filename, "r");
    if (!f) {
        Error("Unable to open file \"%s\"", filename);
        return false;
    }

    // 逐字符扫描文件，解析浮点数
    int c;
    bool inNumber = false;       // 是否正在解析一个数字
    char curNumber[32];          // 当前数字的字符缓冲区
    int curNumberPos = 0;        // 当前数字缓冲区写入位置
    int lineNumber = 1;          // 当前行号（用于错误报告）
    while ((c = getc(f)) != EOF) {
        if (c == '\n') ++lineNumber;
        if (inNumber) {
            // 正在解析数字：收集数字、小数点、指数标记等字符
            CHECK_LT(curNumberPos, (int)sizeof(curNumber))
                << "Overflowed buffer for parsing number in file: " << filename
                << ", at line " << lineNumber;
            if (isdigit(c) || c == '.' || c == 'e' || c == '-' || c == '+')
                curNumber[curNumberPos++] = c;
            else {
                curNumber[curNumberPos++] = '\0';
                values->push_back(atof(curNumber));
                inNumber = false;
                curNumberPos = 0;
            }
        } else {
            // 非数字状态：判断新字符是否为数字起始字符
            if (isdigit(c) || c == '.' || c == '-' || c == '+') {
                inNumber = true;
                curNumber[curNumberPos++] = c;
            } else if (c == '#') {
                // # 开头表示注释，跳过整行
                while ((c = getc(f)) != '\n' && c != EOF)
                    ;
                ++lineNumber;
            } else if (!isspace(c)) {
                // 遇到非空白字符且非数字起始，发出警告
                Warning("Unexpected text found at line %d of float file \"%s\"",
                        lineNumber, filename);
            }
        }
    }
    fclose(f);
    return true;
}

}  // namespace pbrt
