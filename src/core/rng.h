
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

#ifndef PBRT_CORE_RNG_H
#define PBRT_CORE_RNG_H

// core/rng.h*
// RNG: 基于PCG(permuted congruential generator)算法的伪随机数生成器，
// 提供高质量、高效的随机数生成，支持序列跳跃和状态比较
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
static const double DoubleOneMinusEpsilon = 0.99999999999999989;
static const float FloatOneMinusEpsilon = 0.99999994;
#else
static const double DoubleOneMinusEpsilon = 0x1.fffffffffffffp-1;
static const float FloatOneMinusEpsilon = 0x1.fffffep-1;
#endif

#ifdef PBRT_FLOAT_AS_DOUBLE
static const Float OneMinusEpsilon = DoubleOneMinusEpsilon;
#else
static const Float OneMinusEpsilon = FloatOneMinusEpsilon;
#endif

#define PCG32_DEFAULT_STATE 0x853c49e6748fea9bULL
#define PCG32_DEFAULT_STREAM 0xda3e39cb94b95bdbULL
#define PCG32_MULT 0x5851f42d4c957f2dULL
// RNG: PCG伪随机数生成器，提供高质量32位随机数和序列控制
class RNG {
  public:
    // RNG Public Methods
    RNG();  // 默认构造函数，使用预设状态初始化
    RNG(uint64_t sequenceIndex) { SetSequence(sequenceIndex); }  // 使用序列索引初始化
    void SetSequence(uint64_t sequenceIndex);  // 设置随机序列
    uint32_t UniformUInt32();  // 生成32位均匀分布无符号整数
    uint32_t UniformUInt32(uint32_t b) {  // 生成[0,b)范围的均匀分布随机数
        uint32_t threshold = (~b + 1u) % b;
        while (true) {
            uint32_t r = UniformUInt32();
            if (r >= threshold) return r % b;
        }
    }
    Float UniformFloat() {  // 生成[0,1)范围的均匀分布浮点随机数
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
        return std::min(OneMinusEpsilon,
                        Float(UniformUInt32() * 2.3283064365386963e-10f));
#else
        return std::min(OneMinusEpsilon, Float(UniformUInt32() * 0x1p-32f));
#endif
    }
    template <typename Iterator>
    void Shuffle(Iterator begin, Iterator end) {
        for (Iterator it = end - 1; it > begin; --it)
            std::iter_swap(it,
                           begin + UniformUInt32((uint32_t)(it - begin + 1)));
    }
    void Advance(int64_t idelta) {  // 将随机数生成器向前跳跃idelta步
        uint64_t cur_mult = PCG32_MULT, cur_plus = inc, acc_mult = 1u,
                 acc_plus = 0u, delta = (uint64_t)idelta;
        while (delta > 0) {
            if (delta & 1) {
                acc_mult *= cur_mult;
                acc_plus = acc_plus * cur_mult + cur_plus;
            }
            cur_plus = (cur_mult + 1) * cur_plus;
            cur_mult *= cur_mult;
            delta /= 2;
        }
        state = acc_mult * state + acc_plus;
    }
    int64_t operator-(const RNG &other) const {  // 计算两个RNG状态之间的距离（步数）
        CHECK_EQ(inc, other.inc);
        uint64_t cur_mult = PCG32_MULT, cur_plus = inc, cur_state = other.state,
                 the_bit = 1u, distance = 0u;
        while (state != cur_state) {
            if ((state & the_bit) != (cur_state & the_bit)) {
                cur_state = cur_state * cur_mult + cur_plus;
                distance |= the_bit;
            }
            CHECK_EQ(state & the_bit, cur_state & the_bit);
            the_bit <<= 1;
            cur_plus = (cur_mult + 1ULL) * cur_plus;
            cur_mult *= cur_mult;
        }
        return (int64_t)distance;
    }

  private:
    // RNG Private Data
    uint64_t state, inc;  // RNG内部状态：state为当前状态值，inc为序列增量
};

// RNG Inline Method Definitions
inline RNG::RNG() : state(PCG32_DEFAULT_STATE), inc(PCG32_DEFAULT_STREAM) {}
inline void RNG::SetSequence(uint64_t initseq) {
    state = 0u;
    inc = (initseq << 1u) | 1u;
    UniformUInt32();
    state += PCG32_DEFAULT_STATE;
    UniformUInt32();
}

inline uint32_t RNG::UniformUInt32() {
    uint64_t oldstate = state;
    state = oldstate * PCG32_MULT + inc;
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
}

}  // namespace pbrt

#endif  // PBRT_CORE_RNG_H
