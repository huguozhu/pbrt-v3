
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

#ifndef PBRT_CORE_PROGRESSREPORTER_H
#define PBRT_CORE_PROGRESSREPORTER_H

// core/progressreporter.h*
// ProgressReporter: 渲染进度报告器，在控制台显示渲染进度条和预计剩余时间，
// 使用独立线程定期更新显示
#include "pbrt.h"
#include <atomic>
#include <chrono>
#include <thread>

namespace pbrt {

// ProgressReporter Declarations
// ProgressReporter: 渲染进度报告器，显示进度百分比、已用时间和预计剩余时间
class ProgressReporter {
  public:
    // ProgressReporter Public Methods
    // 构造函数：totalWork为总工作量，title为任务名称
    ProgressReporter(int64_t totalWork, const std::string &title);
    ~ProgressReporter();
    // Update: 增加已完成的工作量（默认为1）
    void Update(int64_t num = 1) {
        if (num == 0 || PbrtOptions.quiet) return;
        workDone += num;
    }
    Float ElapsedMS() const {
        // 返回从开始到现在的毫秒数
        std::chrono::system_clock::time_point now =
            std::chrono::system_clock::now();
        int64_t elapsedMS =
            std::chrono::duration_cast<std::chrono::milliseconds>(now -
                                                                  startTime)
                .count();
        return (Float)elapsedMS;
    }
    void Done();

  private:
    // ProgressReporter Private Methods
    void PrintBar();  // 打印进度条

    // ProgressReporter Private Data
    const int64_t totalWork;              // 总工作量
    const std::string title;              // 任务标题
    const std::chrono::system_clock::time_point startTime;  // 开始时间
    std::atomic<int64_t> workDone;        // 已完成工作量（原子操作，线程安全）
    std::atomic<bool> exitThread;         // 退出线程标记
    std::thread updateThread;             // 更新显示的后台线程
};

}  // namespace pbrt

#endif  // PBRT_CORE_PROGRESSREPORTER_H
