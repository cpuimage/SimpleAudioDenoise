// fast fourier transform library (suitable for power of 2 and non-power of 2) Public domain. See "unlicense" statement at the end of this file.
// stb_fft - v0.11 - 2018-12-28
//
// ZhiHan Gao - 200759103@qq.com
// USAGE
//
// This is a single-file library. To use it, do something like the following in one .c file.
//     #define STB_FFT_IMPLEMENTAION
//     #include "stb_fft.h"
//
//
// You can then #include this file in other parts of the program as you would with any other header file.
//
// The default real type is float.
// To define a macro named USE_DOUBLE_TYPE to change it to double.
// You can read the following function: STB_FFT,STB_IFFT,STB_FFT_R2C,STB_IFFT_C2R for easily usage
//
//

#ifndef _stb_fft_h_
#define _stb_fft_h_
#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_DOUBLE_TYPE
#define stb_real_t double
#else
#define stb_real_t float
#endif

typedef struct {
    stb_real_t real;
    stb_real_t imag;
} cmplx;

typedef struct {
    int count;
    int* radix;
    int* remainder;
    int* offsets;
} stb_fft_stages;

typedef struct {
    int N;
    cmplx* twiddles;
    cmplx* twiddles_ordered;
    stb_fft_stages stages;
} stb_fft_plan;

typedef struct {
    stb_fft_plan* half_plan;
    cmplx* buffer;
    cmplx* twiddles;
} stb_fft_real_plan;

stb_fft_plan* stb_fft_plan_dft_1d(int N);

stb_fft_real_plan* stb_fft_real_plan_dft_1d(int N);

void stb_fft_r2c_exec(stb_fft_real_plan* plan, const stb_real_t* in, cmplx* out);

void stb_fft_c2r_exec(stb_fft_real_plan* plan, const cmplx* in, stb_real_t* out);

void stb_fft_exec(const stb_fft_plan* plan, cmplx* in, cmplx* out);

void stb_ifft_exec(const stb_fft_plan* plan, cmplx* in, cmplx* out);

// for easily usage
void STB_FFT(cmplx* input, cmplx* output, int n);

void STB_IFFT(cmplx* input, cmplx* output, int n);

void STB_FFT_R2C(stb_real_t* input, cmplx* output, int n);

void STB_IFFT_C2R(cmplx* input, stb_real_t* output, int n);

#ifdef __cplusplus
}
#endif
#endif

#ifdef STB_FFT_IMPLEMENTAION

#ifndef STB_KP5
#define STB_KP5 ((stb_real_t)0.5)
#endif

#ifndef STB_KP25
#define STB_KP25 ((stb_real_t)0.25)
#endif

#ifndef STB_KP951056516
#define STB_KP951056516 ((stb_real_t)0.951056516295153572116439333379382143405698634)
#endif

#ifndef STB_KP587785252
#define STB_KP587785252 ((stb_real_t)0.587785252292473129168705954639072768597652438)
#endif

#ifndef STB_KP559016994
#define STB_KP559016994 ((stb_real_t)0.559016994374947424102293417182819058860154590)
#endif

#ifndef STB_KP866025403
#define STB_KP866025403 ((stb_real_t)0.866025403784438646763723170752936183471402627)
#endif

#ifndef STB_KP900968867
#define STB_KP900968867 ((stb_real_t)0.900968867902419126236102319507445051165919162)
#endif

#ifndef STB_KP222520933
#define STB_KP222520933 ((stb_real_t)0.222520933956314404288902564496794759466355569)
#endif

#ifndef STB_KP623489801
#define STB_KP623489801 ((stb_real_t)0.623489801858733530525004884004239810632274731)
#endif

#ifndef STB_KP781831482
#define STB_KP781831482 ((stb_real_t)0.781831482468029808708444526674057750232334519)
#endif

#ifndef STB_KP974927912
#define STB_KP974927912 ((stb_real_t)0.974927912181823607018131682993931217232785801)
#endif

#ifndef STB_KP433883739
#define STB_KP433883739 ((stb_real_t)0.433883739117558120475768332848358754609990728)
#endif

#ifndef STB_KP707106781
#define STB_KP707106781 ((stb_real_t)0.707106781186547524400844362104849039284835938)
#endif

#ifndef STB_TWOPI
#define STB_TWOPI ((stb_real_t)6.283185307179586476925286766559005768394338798750211641949889)
#endif

#ifndef STB_PI
#define STB_PI ((stb_real_t)3.141592653589793238462643383279502884197169399375105820974944)
#endif

void stbSinCos(double x, stb_real_t* pSin, stb_real_t* pCos)
{
    *pSin = (stb_real_t)sin(x);
    *pCos = (stb_real_t)cos(x);
}

void stb_make_twiddles(int n, int count, cmplx* twiddles)
{
    double w_pi = STB_TWOPI / n;
    double t = 0;
    for (int i = 0; i < count; ++i) {
        stbSinCos(t, &twiddles[i].imag, &twiddles[i].real);
        t += w_pi;
    }
};

int stb_make_twiddles_sequential(int n, cmplx* twiddles, stb_fft_stages* stages)
{
    int size = 0;
    {
        int offset = 0;
        for (int s = 0; s < stages->count; s++) {
            int r = stages->radix[s];
            int count = stages->remainder[s];
            int amount = (r <= 8) ? (r - 1) * (count - 1) : 0;
            stages->offsets[s] = amount;
            offset += amount;
        }
        size = offset;
        for (int s = 0; s < stages->count; s++) {
            int count = stages->offsets[s];
            offset -= count;
            stages->offsets[s] = offset;
        }
    }
    if (twiddles) {
        int w = 1;
        for (int s = 0; s < stages->count; s++) {
            double w_pi = STB_TWOPI * w / n;
            const int r = stages->radix[s];
            const int count = stages->remainder[s];
            int offset = stages->offsets[s];
            if (r <= 8) {
                for (int i = 1; i < count; i++) {
                    for (int j = 1; j < r; j++) {
                        stbSinCos(w_pi * i * j, &twiddles[offset].imag, &twiddles[offset].real);
                        offset++;
                    }
                }
            }
            w *= r;
        }
    }
    return size;
};

typedef struct {
    int calc_twiddles;
    int stage_count;
    int twiddles_size;
} stb_stage_info;

stb_stage_info stb_calculate_stages(int n, stb_fft_plan* plan)
{
    int calc_twiddles = 0;
    int stage = 0;
    int twiddles_size = 0;
    while (n > 1) {
        int i = 8;
        for (; i > 1; i--) {
            if (!(n % i)) {
                twiddles_size += ((i - 1) * (n - 1));
                break;
            }
        }
        if (i == 1) {
            calc_twiddles = 1;
            i = 7;
            for (; i <= n; i++) {
                if (!(n % i)) {
                    break;
                }
            }
        }
        n /= i;
        if (plan) {
            plan->stages.radix[stage] = i;
            plan->stages.remainder[stage] = n;
        }
        stage++;
    }
    stb_stage_info result = {
        calc_twiddles, stage, twiddles_size
    };
    return result;
}

int stb_make_fft_plan(int N, stb_fft_plan* plan)
{
    if (N == 0)
        return 0;
    stb_stage_info info = stb_calculate_stages(N, NULL);
    const int size_plan = sizeof(stb_fft_plan);
    const int size_radix = info.stage_count * sizeof(int);
    const int size_remainder = info.stage_count * sizeof(int);
    const int size_offsets = info.stage_count * sizeof(int);
    const int twiddles_ordered_size = info.twiddles_size * sizeof(cmplx);
    const int twiddles_size = N * sizeof(cmplx) * info.calc_twiddles;
    const int total_size = size_plan
        + twiddles_size
        + size_radix + size_remainder + size_offsets + twiddles_ordered_size;
    if (plan) {
        uint8_t* data = (uint8_t*)(plan);
        uint8_t* data_twiddles = data + size_plan;
        uint8_t* data_radix = data_twiddles + twiddles_size;
        uint8_t* data_remainder = data_radix + size_radix;
        uint8_t* data_offsets = data_remainder + size_offsets;
        uint8_t* data_twiddles_ordered = data_offsets + size_remainder;
        plan->twiddles = (info.calc_twiddles) ? (cmplx*)(data_twiddles) : NULL;
        plan->stages.radix = (int*)(data_radix);
        plan->stages.remainder = (int*)(data_remainder);
        plan->stages.offsets = (int*)(data_offsets);
        plan->stages.count = info.stage_count;
        plan->twiddles_ordered = (cmplx*)(data_twiddles_ordered);
        plan->N = N;
        if (plan->twiddles) {
            stb_make_twiddles(N, N, plan->twiddles);
        }
        stb_calculate_stages(N, plan);
        stb_make_twiddles_sequential(N, plan->twiddles_ordered, &plan->stages);
    }
    return total_size;
}

stb_fft_plan* stb_fft_plan_dft_1d(int N)
{
    if (N == 0)
        return 0;
    stb_stage_info info = stb_calculate_stages(N, NULL);
    const int size_plan = sizeof(stb_fft_plan);
    const int size_radix = info.stage_count * sizeof(int);
    const int size_remainder = info.stage_count * sizeof(int);
    const int size_offsets = info.stage_count * sizeof(int);
    const int twiddles_ordered_size = info.twiddles_size * sizeof(cmplx);
    const int twiddles_size = N * sizeof(cmplx) * info.calc_twiddles;
    const int total_size = size_plan
        + twiddles_size
        + size_radix + size_remainder + size_offsets + twiddles_ordered_size;
    stb_fft_plan* plan = (stb_fft_plan*)calloc(total_size, 1);
    if (plan) {
        uint8_t* data = (uint8_t*)(plan);
        uint8_t* data_twiddles = data + size_plan;
        uint8_t* data_radix = data_twiddles + twiddles_size;
        uint8_t* data_remainder = data_radix + size_radix;
        uint8_t* data_offsets = data_remainder + size_offsets;
        uint8_t* data_twiddles_ordered = data_offsets + size_remainder;
        plan->twiddles = (info.calc_twiddles) ? (cmplx*)(data_twiddles) : NULL;
        plan->stages.radix = (int*)(data_radix);
        plan->stages.remainder = (int*)(data_remainder);
        plan->stages.offsets = (int*)(data_offsets);
        plan->stages.count = info.stage_count;
        plan->twiddles_ordered = (cmplx*)(data_twiddles_ordered);
        plan->N = N;
        if (plan->twiddles) {
            stb_make_twiddles(N, N, plan->twiddles);
        }
        stb_calculate_stages(N, plan);
        stb_make_twiddles_sequential(N, plan->twiddles_ordered, &plan->stages);
    }
    return plan;
}

stb_fft_real_plan* stb_fft_real_plan_dft_1d(int N)
{
    if (N == 0)
        return 0;
    if (N & 1) {
        fprintf(stderr, "Real FFT must be even.\n");
        return 0;
    }
    N >>= 1;
    int c_plan_size = stb_make_fft_plan(N, NULL);
    int mem_needed = sizeof(stb_fft_real_plan) + c_plan_size + sizeof(cmplx) * (N * 3 / 2);
    stb_fft_real_plan* plan = (stb_fft_real_plan*)calloc(mem_needed, 1);
    if (plan) {
        plan->half_plan = (stb_fft_plan*)(plan + 1);
        plan->buffer = (cmplx*)(((char*)plan->half_plan) + c_plan_size);
        plan->twiddles = plan->buffer + N;
        stb_make_fft_plan(N, plan->half_plan);
        for (int i = 0; i < N / 2; ++i) {
            stbSinCos(-STB_PI * ((double)(i + 1) / N + 0.5), &plan->twiddles[i].imag, &plan->twiddles[i].real);
        }
    }
    return plan;
}

void stb_general_dit(cmplx* twiddles, cmplx* out, int count,
    int cur_sign, int radix, int N, int reverse)
{
    cmplx* scratch = (cmplx*)calloc(radix, sizeof(cmplx));
    if (scratch == 0)
        return;
    if (reverse) {
        for (int butterfly = 0; butterfly < count; ++butterfly) {
            int idx = butterfly;
            for (int r = 0; r < radix; ++r) {
                scratch[r] = out[idx];
                idx += count;
            }
            idx = butterfly;
            for (int r = 0; r < radix; ++r) {
                int tw_idx = 0;
                out[idx] = *scratch;
                for (int cr = 1; cr < radix; ++cr) {
                    tw_idx += cur_sign * idx;
                    if (tw_idx >= N)
                        tw_idx -= N;
                    out[idx].real += twiddles[tw_idx].real * scratch[cr].real
                        - twiddles[tw_idx].imag * scratch[cr].imag;
                    out[idx].imag += twiddles[tw_idx].real * scratch[cr].imag
                        + scratch[cr].real * twiddles[tw_idx].imag;
                }
                idx += count;
            }
        }
    } else {
        for (int butterfly = 0; butterfly < count; ++butterfly) {
            int idx = butterfly;
            for (int r = 0; r < radix; ++r) {
                scratch[r] = out[idx];
                idx += count;
            }
            idx = butterfly;
            for (int r = 0; r < radix; ++r) {
                int tw_idx = 0;
                out[idx] = *scratch;
                for (int cr = 1; cr < radix; ++cr) {
                    tw_idx += cur_sign * idx;
                    if (tw_idx >= N)
                        tw_idx -= N;
                    out[idx].real += twiddles[tw_idx].real * scratch[cr].real
                        + twiddles[tw_idx].imag * scratch[cr].imag;
                    out[idx].imag += twiddles[tw_idx].real * scratch[cr].imag
                        - scratch[cr].real * twiddles[tw_idx].imag;
                }
                idx += count;
            }
        }
    }
    free(scratch);
}

void stb_radix_2_dit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T1 = out[0].real;
    stb_real_t T2 = out[count].real;
    out[count].real = out[0].real - T2;
    out[0].real = T1 + T2;
    stb_real_t T3 = out[0].imag;
    stb_real_t T4 = out[count].imag;
    out[count].imag = T3 - T4;
    out[0].imag = T3 + T4;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T8 = output[0].imag;
        stb_real_t T3_0 = output[count].real;
        stb_real_t T5 = output[count].imag;
        stb_real_t T4_0 = tw[0].imag;
        stb_real_t T6 = T4_0 * T5 + tw[0].real * T3_0;
        stb_real_t T7 = tw[0].real * T5 - T4_0 * T3_0;
        output[count].real = output[0].real - T6;
        output[count].imag = T8 - T7;
        output[0].real = T1_0 + T6;
        output[0].imag = T7 + T8;
        ++m;
        ++output;
        ++tw;
    }
}

void stb_radix_3_dit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T1 = out[0].real;
    stb_real_t T10 = out[0].imag;
    stb_real_t T2 = out[count].real;
    stb_real_t T3 = out[2 * count].real;
    stb_real_t T9 = (T3 - T2) * STB_KP866025403;
    stb_real_t T6 = out[count].imag;
    stb_real_t T7 = out[2 * count].imag;
    stb_real_t T8 = (T6 - T7) * STB_KP866025403;
    out[0].real = out[0].real + T2 + T3;
    out[0].imag = T10 + T6 + T7;
    stb_real_t T5 = T1 - STB_KP5 * (T2 + T3);
    out[2 * count].real = T5 - T8;
    out[count].real = T5 + T8;
    stb_real_t T12 = T10 - STB_KP5 * (T6 + T7);
    out[count].imag = T9 + T12;
    out[2 * count].imag = T12 - T9;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T18 = output[0].imag;
        stb_real_t T3_0 = output[count].real;
        stb_real_t T5_0 = output[count].imag;
        stb_real_t T4_0 = tw[0].imag;
        stb_real_t T6_0 = T4_0 * T5_0 + tw[0].real * T3_0;
        stb_real_t T14 = tw[0].real * T5_0 - T4_0 * T3_0;
        stb_real_t T8_0 = output[2 * count].real;
        stb_real_t T10_0 = output[2 * count].imag;
        stb_real_t T7_0 = tw[1].real;
        stb_real_t T9_0 = tw[1].imag;
        stb_real_t T11_0 = T9_0 * T10_0 + T7_0 * T8_0;
        stb_real_t T15 = T7_0 * T10_0 - T9_0 * T8_0;
        stb_real_t T12_0 = T6_0 + T11_0;
        stb_real_t T17 = T14 + T15;
        output[0].real = output[0].real + T12_0;
        output[0].imag = T17 + T18;
        stb_real_t T16 = (T14 - T15) * STB_KP866025403;
        output[2 * count].real = T1_0 - STB_KP5 * T12_0 - T16;
        output[count].real = T1_0 - STB_KP5 * T12_0 + T16;
        stb_real_t T20 = T18 - STB_KP5 * T17;
        output[count].imag = (T11_0 - T6_0) * STB_KP866025403 + T20;
        output[2 * count].imag = T20 - (T11_0 - T6_0) * STB_KP866025403;
        ++m;
        ++output;
        tw += 2;
    }
}

void stb_radix_4_dit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T2 = out[2 * count].real;
    stb_real_t T3 = out[0].real + T2;
    stb_real_t T11 = out[0].real - T2;
    stb_real_t T7 = out[0].imag;
    stb_real_t T8 = out[2 * count].imag;
    stb_real_t T4 = out[count].real;
    stb_real_t T5 = out[3 * count].real;
    stb_real_t T12 = out[count].imag;
    stb_real_t T13 = out[3 * count].imag;
    out[2 * count].real = T3 - (T4 + T5);
    out[2 * count].imag = T7 + T8 - (T12 + T13);
    out[0].real = T3 + T4 + T5;
    out[0].imag = T7 + T8 + T12 + T13;
    out[count].imag = T7 - T8 - (T4 - T5);
    out[count].real = T11 + T12 - T13;
    out[3 * count].imag = T4 - T5 + T7 - T8;
    out[3 * count].real = T11 - (T12 - T13);
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T25 = output[0].imag;
        stb_real_t T3_0 = output[2 * count].real;
        stb_real_t T5_0 = output[2 * count].imag;
        stb_real_t T2_0 = tw[1].real;
        stb_real_t T4_0 = tw[1].imag;
        stb_real_t T6_0 = T4_0 * T5_0 + T2_0 * T3_0;
        stb_real_t T24 = T2_0 * T5_0 - T4_0 * T3_0;
        stb_real_t T9_0 = output[count].real;
        stb_real_t T11_0 = output[count].imag;
        stb_real_t T10_0 = tw[0].imag;
        stb_real_t T12_0 = T10_0 * T11_0 + tw[0].real * T9_0;
        stb_real_t T20 = tw[0].real * T11_0 - T10_0 * T9_0;
        stb_real_t T14_0 = output[3 * count].real;
        stb_real_t T16_0 = output[3 * count].imag;
        stb_real_t T13_0 = tw[2].real;
        stb_real_t T15_0 = tw[2].imag;
        stb_real_t T17 = T15_0 * T16_0 + T13_0 * T14_0;
        stb_real_t T21 = T13_0 * T16_0 - T15_0 * T14_0;
        stb_real_t T18 = T12_0 + T17;
        output[2 * count].real = output[0].real + T6_0 - T18;
        output[0].real = T1_0 + T6_0 + T18;
        output[0].imag = T20 + T21 + T24 + T25;
        output[2 * count].imag = T24 + T25 - (T20 + T21);
        output[3 * count].real = T1_0 - T6_0 - (T20 - T21);
        output[count].real = T1_0 - T6_0 + T20 - T21;
        output[count].imag = T25 - T24 - (T12_0 - T17);
        output[3 * count].imag = T12_0 - T17 + T25 - T24;
        ++m;
        ++output;
        tw += 3;
    }
}

void stb_radix_5_dit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T1 = out[0].real;
    stb_real_t T24 = out[0].imag;
    stb_real_t T2 = out[count].real;
    stb_real_t T3 = out[4 * count].real;
    stb_real_t T5 = out[2 * count].real;
    stb_real_t T6 = out[3 * count].real;
    stb_real_t T8 = T2 + T3 + T5 + T6;
    stb_real_t T9 = (T2 + T3 - (T5 + T6)) * STB_KP559016994;
    stb_real_t T12 = out[count].imag;
    stb_real_t T13 = out[4 * count].imag;
    stb_real_t T15 = out[2 * count].imag;
    stb_real_t T16 = out[3 * count].imag;
    stb_real_t T25 = T12 + T13 + T15 + T16;
    stb_real_t T23 = (T12 + T13 - (T15 + T16)) * STB_KP559016994;
    out[0].real = out[0].real + T8;
    out[0].imag = T24 + T25;
    stb_real_t T18 = STB_KP587785252 * (T15 - T16) + STB_KP951056516 * (T12 - T13);
    stb_real_t T20 = STB_KP951056516 * (T15 - T16) - STB_KP587785252 * (T12 - T13);
    stb_real_t T11 = T9 + T1 - STB_KP25 * T8;
    stb_real_t T19 = T1 - STB_KP25 * T8 - T9;
    out[4 * count].real = T11 - T18;
    out[3 * count].real = T19 + T20;
    out[count].real = T11 + T18;
    out[2 * count].real = T19 - T20;
    stb_real_t T30 = STB_KP587785252 * (T5 - T6) + STB_KP951056516 * (T2 - T3);
    stb_real_t T31 = STB_KP951056516 * (T5 - T6) - STB_KP587785252 * (T2 - T3);
    stb_real_t T27 = T23 + T24 - STB_KP25 * T25;
    stb_real_t T32 = T24 - STB_KP25 * T25 - T23;
    out[count].imag = T27 - T30;
    out[3 * count].imag = T32 - T31;
    out[4 * count].imag = T30 + T27;
    out[2 * count].imag = T31 + T32;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T40 = output[0].imag;
        stb_real_t T3_0 = output[count].real;
        stb_real_t T5_0 = output[count].imag;
        stb_real_t T4_0 = tw[0].imag;
        stb_real_t T19_0 = output[3 * count].real;
        stb_real_t T21_0 = output[3 * count].imag;
        stb_real_t T18_0 = tw[2].real;
        stb_real_t T20_0 = tw[2].imag;
        stb_real_t T22_0 = T20_0 * T21_0 + T18_0 * T19_0;
        stb_real_t T32_0 = T18_0 * T21_0 - T20_0 * T19_0;
        stb_real_t T8_0 = output[4 * count].real;
        stb_real_t T10_0 = output[4 * count].imag;
        stb_real_t T7_0 = tw[3].real;
        stb_real_t T9_0 = tw[3].imag;
        stb_real_t T11_0 = T9_0 * T10_0 + T7_0 * T8_0;
        stb_real_t T29_0 = T7_0 * T10_0 - T9_0 * T8_0;
        stb_real_t T14_0 = output[2 * count].real;
        stb_real_t T16_0 = output[2 * count].imag;
        stb_real_t T13_0 = tw[1].real;
        stb_real_t T15_0 = tw[1].imag;
        stb_real_t T30_0 = tw[0].real * T5_0 - T4_0 * T3_0 - T29_0;
        stb_real_t T33 = T13_0 * T16_0 - T15_0 * T14_0 - T32_0;
        stb_real_t T45 = T15_0 * T16_0 + T13_0 * T14_0 - T22_0;
        stb_real_t T44 = T4_0 * T5_0 + tw[0].real * T3_0 - T11_0;
        stb_real_t T37 = tw[0].real * T5_0 - T4_0 * T3_0 + T29_0;
        stb_real_t T38 = T13_0 * T16_0 - T15_0 * T14_0 + T32_0;
        stb_real_t T39 = T37 + T38;
        stb_real_t T12_0 = T4_0 * T5_0 + tw[0].real * T3_0 + T11_0;
        stb_real_t T23_0 = T15_0 * T16_0 + T13_0 * T14_0 + T22_0;
        stb_real_t T24_0 = T12_0 + T23_0;
        output[0].real = output[0].real + T24_0;
        output[0].imag = T39 + T40;
        stb_real_t T34 = STB_KP587785252 * T33 + STB_KP951056516 * T30_0;
        stb_real_t T36 = STB_KP951056516 * T33 - STB_KP587785252 * T30_0;
        stb_real_t T27_0 = (T12_0 - T23_0) * STB_KP559016994 + T1_0 - STB_KP25 * T24_0;
        stb_real_t T35 = T1_0 - STB_KP25 * T24_0 - (T12_0 - T23_0) * STB_KP559016994;
        output[4 * count].real = T27_0 - T34;
        output[3 * count].real = T35 + T36;
        output[count].real = T27_0 + T34;
        output[2 * count].real = T35 - T36;
        stb_real_t T46 = STB_KP587785252 * T45 + STB_KP951056516 * T44;
        stb_real_t T47 = STB_KP951056516 * T45 - STB_KP587785252 * T44;
        stb_real_t T43 = (T37 - T38) * STB_KP559016994 + T40 - STB_KP25 * T39;
        stb_real_t T48 = T40 - STB_KP25 * T39 - (T37 - T38) * STB_KP559016994;
        output[count].imag = T43 - T46;
        output[3 * count].imag = T48 - T47;
        output[4 * count].imag = T46 + T43;
        output[2 * count].imag = T47 + T48;
        ++m;
        ++output;
        tw += 4;
    }
}

void stb_radix_6_dit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T2 = out[3 * count].real;
    stb_real_t T3 = out[0].real - T2;
    stb_real_t T11 = out[0].real + T2;
    stb_real_t T24 = out[0].imag;
    stb_real_t T25 = out[3 * count].imag;
    stb_real_t T4 = out[2 * count].real;
    stb_real_t T5 = out[5 * count].real;
    stb_real_t T7 = out[4 * count].real;
    stb_real_t T8 = out[count].real;
    stb_real_t T10 = T4 - T5 + T7 - T8;
    stb_real_t T14 = T4 + T5 + T7 + T8;
    stb_real_t T16 = out[2 * count].imag;
    stb_real_t T17 = out[5 * count].imag;
    stb_real_t T19 = out[4 * count].imag;
    stb_real_t T20 = out[count].imag;
    stb_real_t T27 = T16 - T17 + T19 - T20;
    stb_real_t T34 = T16 + T17 + T19 + T20;
    out[3 * count].real = T3 + T10;
    out[3 * count].imag = T24 - T25 + T27;
    out[0].real = T11 + T14;
    out[0].imag = T24 + T25 + T34;
    stb_real_t T22 = (T16 - T17 - (T19 - T20)) * STB_KP866025403;
    out[5 * count].real = T3 - STB_KP5 * T10 - T22;
    out[count].real = T3 - STB_KP5 * T10 + T22;
    stb_real_t T23 = (T7 - T8 - (T4 - T5)) * STB_KP866025403;
    stb_real_t T28 = T24 - T25 - STB_KP5 * T27;
    out[count].imag = T23 + T28;
    out[5 * count].imag = T28 - T23;
    stb_real_t T32 = (T16 + T17 - (T19 + T20)) * STB_KP866025403;
    out[2 * count].real = T11 - STB_KP5 * T14 - T32;
    out[4 * count].real = T11 - STB_KP5 * T14 + T32;
    stb_real_t T35 = T24 + T25 - STB_KP5 * T34;
    stb_real_t T36 = (T7 + T8 - (T4 + T5)) * STB_KP866025403;
    out[2 * count].imag = T35 - T36;
    out[4 * count].imag = T36 + T35;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T49 = output[0].imag;
        stb_real_t T3_0 = output[3 * count].real;
        stb_real_t T5_0 = output[3 * count].imag;
        stb_real_t T2_0 = tw[2].real;
        stb_real_t T4_0 = tw[2].imag;
        stb_real_t T6_0 = T4_0 * T5_0 + T2_0 * T3_0;
        stb_real_t T54 = T49 - (T2_0 * T5_0 - T4_0 * T3_0);
        stb_real_t T31_0 = output[0].real + T6_0;
        stb_real_t T50 = T2_0 * T5_0 - T4_0 * T3_0 + T49;
        stb_real_t T20_0 = output[4 * count].real;
        stb_real_t T22_0 = output[4 * count].imag;
        stb_real_t T19_0 = tw[3].real;
        stb_real_t T21_0 = tw[3].imag;
        stb_real_t T25_0 = output[count].real;
        stb_real_t T27_0 = output[count].imag;
        stb_real_t T26_0 = tw[0].imag;
        stb_real_t T28_0 = T26_0 * T27_0 + tw[0].real * T25_0;
        stb_real_t T40 = tw[0].real * T27_0 - T26_0 * T25_0;
        stb_real_t T29_0 = T21_0 * T22_0 + T19_0 * T20_0 - T28_0;
        stb_real_t T45 = T19_0 * T22_0 - T21_0 * T20_0 + T40;
        stb_real_t T33_0 = T21_0 * T22_0 + T19_0 * T20_0 + T28_0;
        stb_real_t T41 = T19_0 * T22_0 - T21_0 * T20_0 - T40;
        stb_real_t T9_0 = output[2 * count].real;
        stb_real_t T11_0 = output[2 * count].imag;
        stb_real_t T8_0 = tw[1].real;
        stb_real_t T10_0 = tw[1].imag;
        stb_real_t T14_0 = output[5 * count].real;
        stb_real_t T16_0 = output[5 * count].imag;
        stb_real_t T13_0 = tw[4].real;
        stb_real_t T15_0 = tw[4].imag;
        stb_real_t T17_0 = T15_0 * T16_0 + T13_0 * T14_0;
        stb_real_t T37 = T13_0 * T16_0 - T15_0 * T14_0;
        stb_real_t T44 = T8_0 * T11_0 - T10_0 * T9_0 + T37;
        stb_real_t T32_0 = T10_0 * T11_0 + T8_0 * T9_0 + T17_0;
        stb_real_t T42 = (T8_0 * T11_0 - T10_0 * T9_0 - T37 - T41) * STB_KP866025403;
        stb_real_t T30_0 = T10_0 * T11_0 + T8_0 * T9_0 - T17_0 + T29_0;
        stb_real_t T35_0 = output[0].real - T6_0 - STB_KP5 * T30_0;
        output[3 * count].real = output[0].real - T6_0 + T30_0;
        output[count].real = T35_0 + T42;
        output[5 * count].real = T35_0 - T42;
        stb_real_t T53 = (T29_0 - (T10_0 * T11_0 + T8_0 * T9_0 - T17_0)) * STB_KP866025403;
        stb_real_t T55 = T8_0 * T11_0 - T10_0 * T9_0 - T37 + T41;
        stb_real_t T56 = T54 - STB_KP5 * T55;
        output[count].imag = T53 + T56;
        output[3 * count].imag = T55 + T54;
        output[5 * count].imag = T56 - T53;
        stb_real_t T46 = (T44 - T45) * STB_KP866025403;
        stb_real_t T43 = T31_0 - STB_KP5 * (T32_0 + T33_0);
        output[0].real = T31_0 + T32_0 + T33_0;
        output[4 * count].real = T43 + T46;
        output[2 * count].real = T43 - T46;
        stb_real_t T52 = (T33_0 - T32_0) * STB_KP866025403;
        stb_real_t T51 = T50 - STB_KP5 * (T44 + T45);
        output[0].imag = T44 + T45 + T50;
        output[4 * count].imag = T52 + T51;
        output[2 * count].imag = T51 - T52;
        ++m;
        ++output;
        tw += 5;
    }
}

void stb_radix_7_dit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T1 = out[0].real;
    stb_real_t T30 = out[0].imag;
    stb_real_t T2 = out[count].real;
    stb_real_t T3 = out[6 * count].real;
    stb_real_t T12 = out[count].imag;
    stb_real_t T13 = out[6 * count].imag;
    stb_real_t T5 = out[2 * count].real;
    stb_real_t T6 = out[5 * count].real;
    stb_real_t T18 = out[2 * count].imag;
    stb_real_t T19 = out[5 * count].imag;
    stb_real_t T8 = out[3 * count].real;
    stb_real_t T9 = out[4 * count].real;
    stb_real_t T15 = out[3 * count].imag;
    stb_real_t T16 = out[4 * count].imag;
    out[0].real = out[0].real + T2 + T3 + T5 + T6 + T8 + T9;
    out[0].imag = T30 + T12 + T13 + T18 + T19 + T15 + T16;
    stb_real_t T21 = STB_KP974927912 * (T12 - T13) - STB_KP781831482 * (T15 - T16) - STB_KP433883739 * (T18 - T19);
    out[5 * count].real = STB_KP623489801 * (T8 + T9)
        + T1
        - (STB_KP222520933 * (T2 + T3)
              + STB_KP900968867 * (T5 + T6))
        - T21;
    out[2 * count].real = STB_KP623489801 * (T8 + T9)
        + T1
        - (STB_KP222520933 * (T2 + T3)
              + STB_KP900968867 * (T5 + T6))
        + T21;
    stb_real_t T38 = STB_KP623489801 * (T15 + T16) + T30 - (STB_KP222520933 * (T12 + T13) + STB_KP900968867 * (T18 + T19));
    out[2 * count].imag = STB_KP974927912 * (T3 - T2)
        - STB_KP781831482 * (T9 - T8)
        - STB_KP433883739 * (T6 - T5)
        + T38;
    out[5 * count].imag = T38
        - (STB_KP974927912 * (T3 - T2)
              - STB_KP781831482 * (T9 - T8)
              - STB_KP433883739 * (T6 - T5));
    stb_real_t T23 = STB_KP433883739 * (T15 - T16) + STB_KP781831482 * (T12 - T13) + STB_KP974927912 * (T18 - T19);
    out[6 * count].real = STB_KP623489801 * (T2 + T3)
        + T1
        - (STB_KP222520933 * (T5 + T6)
              + STB_KP900968867 * (T8 + T9))
        - T23;
    out[count].real = STB_KP623489801 * (T2 + T3)
        + T1
        - (STB_KP222520933 * (T5 + T6)
              + STB_KP900968867 * (T8 + T9))
        + T23;
    stb_real_t T36 = STB_KP623489801 * (T12 + T13) + T30 - (STB_KP222520933 * (T18 + T19) + STB_KP900968867 * (T15 + T16));
    out[count].imag = STB_KP433883739 * (T9 - T8)
        + STB_KP781831482 * (T3 - T2)
        + STB_KP974927912 * (T6 - T5)
        + T36;
    out[6 * count].imag = T36
        - (STB_KP433883739 * (T9 - T8)
              + STB_KP781831482 * (T3 - T2)
              + STB_KP974927912 * (T6 - T5));
    stb_real_t T25 = STB_KP974927912 * (T15 - T16) + STB_KP433883739 * (T12 - T13) - STB_KP781831482 * (T18 - T19);
    stb_real_t T24 = STB_KP623489801 * (T5 + T6) + T1 - (STB_KP900968867 * (T2 + T3) + STB_KP222520933 * (T8 + T9));
    out[4 * count].real = T24 - T25;
    out[3 * count].real = T24 + T25;
    stb_real_t T29 = STB_KP974927912 * (T9 - T8) + STB_KP433883739 * (T3 - T2) - STB_KP781831482 * (T6 - T5);
    stb_real_t T34 = STB_KP623489801 * (T18 + T19) + T30 - (STB_KP900968867 * (T12 + T13) + STB_KP222520933 * (T15 + T16));
    out[3 * count].imag = T29 + T34;
    out[4 * count].imag = T34 - T29;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T53 = output[0].imag;
        stb_real_t T3_0 = output[count].real;
        stb_real_t T5_0 = output[count].imag;
        stb_real_t T4_0 = tw[0].imag;
        stb_real_t T8_0 = output[6 * count].real;
        stb_real_t T10_0 = output[6 * count].imag;
        stb_real_t T7_0 = tw[5].real;
        stb_real_t T9_0 = tw[5].imag;
        stb_real_t T37_0 = T7_0 * T10_0 - T9_0 * T8_0;
        stb_real_t T12_0 = T4_0 * T5_0 + tw[0].real * T3_0 + T9_0 * T10_0 + T7_0 * T8_0;
        stb_real_t T54 = T9_0 * T10_0 + T7_0 * T8_0 - (T4_0 * T5_0 + tw[0].real * T3_0);
        stb_real_t T38_0 = tw[0].real * T5_0 - T4_0 * T3_0 - T37_0;
        stb_real_t T50 = tw[0].real * T5_0 - T4_0 * T3_0 + T37_0;
        stb_real_t T14_0 = output[2 * count].real;
        stb_real_t T16_0 = output[2 * count].imag;
        stb_real_t T13_0 = tw[1].real;
        stb_real_t T15_0 = tw[1].imag;
        stb_real_t T19_0 = output[5 * count].real;
        stb_real_t T21_0 = output[5 * count].imag;
        stb_real_t T18_0 = tw[4].real;
        stb_real_t T20_0 = tw[4].imag;
        stb_real_t T43 = T18_0 * T21_0 - T20_0 * T19_0;
        stb_real_t T23_0 = T15_0 * T16_0 + T13_0 * T14_0 + T20_0 * T21_0 + T18_0 * T19_0;
        stb_real_t T55 = T20_0 * T21_0 + T18_0 * T19_0 - (T15_0 * T16_0 + T13_0 * T14_0);
        stb_real_t T44 = T13_0 * T16_0 - T15_0 * T14_0 - T43;
        stb_real_t T51 = T13_0 * T16_0 - T15_0 * T14_0 + T43;
        stb_real_t T25_0 = output[3 * count].real;
        stb_real_t T27_0 = output[3 * count].imag;
        stb_real_t T24_0 = tw[2].real;
        stb_real_t T26_0 = tw[2].imag;
        stb_real_t T30_0 = output[4 * count].real;
        stb_real_t T32_0 = output[4 * count].imag;
        stb_real_t T29_0 = tw[3].real;
        stb_real_t T31_0 = tw[3].imag;
        stb_real_t T40 = T29_0 * T32_0 - T31_0 * T30_0;
        stb_real_t T34_0 = T26_0 * T27_0 + T24_0 * T25_0 + T31_0 * T32_0 + T29_0 * T30_0;
        stb_real_t T56 = T31_0 * T32_0 + T29_0 * T30_0 - (T26_0 * T27_0 + T24_0 * T25_0);
        stb_real_t T41 = T24_0 * T27_0 - T26_0 * T25_0 - T40;
        stb_real_t T52 = T24_0 * T27_0 - T26_0 * T25_0 + T40;
        output[0].real = output[0].real + T12_0 + T23_0 + T34_0;
        output[0].imag = T50 + T51 + T52 + T53;
        stb_real_t T45 = STB_KP974927912 * T38_0 - STB_KP781831482 * T41 - STB_KP433883739 * T44;
        output[5 * count].real = STB_KP623489801 * T34_0
            + T1_0
            - (STB_KP222520933 * T12_0
                  + STB_KP900968867 * T23_0)
            - T45;
        output[2 * count].real = STB_KP623489801 * T34_0
            + T1_0
            - (STB_KP222520933 * T12_0
                  + STB_KP900968867 * T23_0)
            + T45;
        stb_real_t T60 = STB_KP623489801 * T52 + T53 - (STB_KP222520933 * T50 + STB_KP900968867 * T51);
        output[2 * count].imag = STB_KP974927912 * T54 - STB_KP781831482 * T56 - STB_KP433883739 * T55 + T60;
        output[5 * count].imag = T60 - (STB_KP974927912 * T54 - STB_KP781831482 * T56 - STB_KP433883739 * T55);
        stb_real_t T47 = STB_KP433883739 * T41 + STB_KP781831482 * T38_0 + STB_KP974927912 * T44;
        output[6 * count].real = STB_KP623489801 * T12_0
            + T1_0
            - (STB_KP222520933 * T23_0
                  + STB_KP900968867 * T34_0)
            - T47;
        output[count].real = STB_KP623489801 * T12_0
            + T1_0
            - (STB_KP222520933 * T23_0
                  + STB_KP900968867 * T34_0)
            + T47;
        stb_real_t T58 = STB_KP623489801 * T50 + T53 - (STB_KP222520933 * T51 + STB_KP900968867 * T52);
        output[count].imag = STB_KP433883739 * T56 + STB_KP781831482 * T54 + STB_KP974927912 * T55 + T58;
        output[6 * count].imag = T58 - (STB_KP433883739 * T56 + STB_KP781831482 * T54 + STB_KP974927912 * T55);
        stb_real_t T49 = STB_KP974927912 * T41 + STB_KP433883739 * T38_0 - STB_KP781831482 * T44;
        output[4 * count].real = STB_KP623489801 * T23_0
            + T1_0
            - (STB_KP900968867 * T12_0
                  + STB_KP222520933 * T34_0)
            - T49;
        output[3 * count].real = STB_KP623489801 * T23_0
            + T1_0
            - (STB_KP900968867 * T12_0
                  + STB_KP222520933 * T34_0)
            + T49;
        stb_real_t T62 = STB_KP623489801 * T51 + T53 - (STB_KP900968867 * T50 + STB_KP222520933 * T52);
        output[3 * count].imag = STB_KP974927912 * T56 + STB_KP433883739 * T54 - STB_KP781831482 * T55 + T62;
        output[4 * count].imag = T62
            - (STB_KP974927912 * T56
                  + STB_KP433883739 * T54
                  - STB_KP781831482 * T55);
        ++m;
        ++output;
        tw += 6;
    }
}

void stb_radix_8_dit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T2 = out[4 * count].real;
    stb_real_t T3 = out[0].real + T2;
    stb_real_t T23 = out[0].real - T2;
    stb_real_t T16 = out[0].imag;
    stb_real_t T17 = out[4 * count].imag;
    stb_real_t T4 = out[2 * count].real;
    stb_real_t T5 = out[6 * count].real;
    stb_real_t T19 = out[2 * count].imag;
    stb_real_t T20 = out[6 * count].imag;
    stb_real_t T11 = out[7 * count].real;
    stb_real_t T12 = out[3 * count].real;
    stb_real_t T32 = out[7 * count].imag;
    stb_real_t T33 = out[3 * count].imag;
    stb_real_t T35 = T11 - T12 - (T32 - T33);
    stb_real_t T43 = T11 - T12 + T32 - T33;
    stb_real_t T8 = out[count].real;
    stb_real_t T9 = out[5 * count].real;
    stb_real_t T27 = out[count].imag;
    stb_real_t T28 = out[5 * count].imag;
    stb_real_t T30 = T8 - T9 + T27 - T28;
    stb_real_t T42 = T27 - T28 - (T8 - T9);
    stb_real_t T7 = T3 + T4 + T5;
    stb_real_t T14 = T8 + T9 + T11 + T12;
    out[4 * count].real = T7 - T14;
    out[0].real = T7 + T14;
    out[4 * count].imag = T16 + T17 + T19 + T20 - (T27 + T28 + T32 + T33);
    out[0].imag = T16 + T17 + T19 + T20 + T27 + T28 + T32 + T33;
    stb_real_t T15 = T11 + T12 - (T8 + T9);
    stb_real_t T22 = T16 + T17 - (T19 + T20);
    out[2 * count].imag = T15 + T22;
    out[6 * count].imag = T22 - T15;
    stb_real_t T47 = T3 - (T4 + T5);
    stb_real_t T50 = T27 + T28 - (T32 + T33);
    out[6 * count].real = T47 - T50;
    out[2 * count].real = T47 + T50;
    stb_real_t T36 = (T30 + T35) * STB_KP707106781;
    out[5 * count].real = T23 + T19 - T20 - T36;
    out[count].real = T23 + T19 - T20 + T36;
    stb_real_t T46 = (T42 + T43) * STB_KP707106781;
    out[5 * count].imag = T16 - T17 - (T4 - T5) - T46;
    out[count].imag = T16 - T17 - (T4 - T5) + T46;
    stb_real_t T39 = T4 - T5 + T16 - T17;
    stb_real_t T40 = (T35 - T30) * STB_KP707106781;
    out[7 * count].imag = T39 - T40;
    out[3 * count].imag = T39 + T40;
    stb_real_t T41 = T23 - (T19 - T20);
    stb_real_t T44 = (T42 - T43) * STB_KP707106781;
    out[7 * count].real = T41 - T44;
    out[3 * count].real = T41 + T44;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T70 = output[0].imag;
        stb_real_t T3_0 = output[4 * count].real;
        stb_real_t T5_0 = output[4 * count].imag;
        stb_real_t T2_0 = tw[3].real;
        stb_real_t T4_0 = tw[3].imag;
        stb_real_t T6_0 = T4_0 * T5_0 + T2_0 * T3_0;
        stb_real_t T7_0 = output[0].real + T6_0;
        stb_real_t T76 = T70 - (T2_0 * T5_0 - T4_0 * T3_0);
        stb_real_t T43_0 = output[0].real - T6_0;
        stb_real_t T71 = T2_0 * T5_0 - T4_0 * T3_0 + T70;
        stb_real_t T32_0 = output[7 * count].real;
        stb_real_t T34_0 = output[7 * count].imag;
        stb_real_t T31_0 = tw[6].real;
        stb_real_t T33_0 = tw[6].imag;
        stb_real_t T37_0 = output[3 * count].real;
        stb_real_t T39_0 = output[3 * count].imag;
        stb_real_t T36_0 = tw[2].real;
        stb_real_t T38_0 = tw[2].imag;
        stb_real_t T40_0 = T38_0 * T39_0 + T36_0 * T37_0;
        stb_real_t T55 = T36_0 * T39_0 - T38_0 * T37_0;
        stb_real_t T41_0 = T33_0 * T34_0 + T31_0 * T32_0 + T40_0;
        stb_real_t T65 = T31_0 * T34_0 - T33_0 * T32_0 + T55;
        stb_real_t T53 = T33_0 * T34_0 + T31_0 * T32_0 - T40_0;
        stb_real_t T56 = T31_0 * T34_0 - T33_0 * T32_0 - T55;
        stb_real_t T9_0 = output[2 * count].real;
        stb_real_t T11_0 = output[2 * count].imag;
        stb_real_t T8_0 = tw[1].real;
        stb_real_t T10_0 = tw[1].imag;
        stb_real_t T14_0 = output[6 * count].real;
        stb_real_t T16_0 = output[6 * count].imag;
        stb_real_t T13_0 = tw[5].real;
        stb_real_t T15_0 = tw[5].imag;
        stb_real_t T17_0 = T15_0 * T16_0 + T13_0 * T14_0;
        stb_real_t T45_0 = T13_0 * T16_0 - T15_0 * T14_0;
        stb_real_t T18_0 = T10_0 * T11_0 + T8_0 * T9_0 + T17_0;
        stb_real_t T77 = T10_0 * T11_0 + T8_0 * T9_0 - T17_0;
        stb_real_t T46_0 = T8_0 * T11_0 - T10_0 * T9_0 - T45_0;
        stb_real_t T68 = T8_0 * T11_0 - T10_0 * T9_0 + T45_0;
        stb_real_t T21_0 = output[count].real;
        stb_real_t T23_0 = output[count].imag;
        stb_real_t T22_0 = tw[0].imag;
        stb_real_t T26_0 = output[5 * count].real;
        stb_real_t T28_0 = output[5 * count].imag;
        stb_real_t T25_0 = tw[4].real;
        stb_real_t T27_0 = tw[4].imag;
        stb_real_t T29_0 = T27_0 * T28_0 + T25_0 * T26_0;
        stb_real_t T50_0 = T25_0 * T28_0 - T27_0 * T26_0;
        stb_real_t T30_0 = T22_0 * T23_0 + tw[0].real * T21_0 + T29_0;
        stb_real_t T64 = tw[0].real * T23_0 - T22_0 * T21_0 + T50_0;
        stb_real_t T48_0 = T22_0 * T23_0 + tw[0].real * T21_0 - T29_0;
        stb_real_t T51_0 = tw[0].real * T23_0 - T22_0 * T21_0 - T50_0;
        stb_real_t T19_0 = T7_0 + T18_0;
        stb_real_t T42_0 = T30_0 + T41_0;
        output[4 * count].real = T19_0 - T42_0;
        output[0].real = T19_0 + T42_0;
        output[0].imag = T64 + T65 + T68 + T71;
        output[4 * count].imag = T68 + T71 - (T64 + T65);
        output[6 * count].real = T7_0 - T18_0 - (T64 - T65);
        output[2 * count].real = T7_0 - T18_0 + T64 - T65;
        output[2 * count].imag = T41_0 - T30_0 + T71 - T68;
        output[6 * count].imag = T71 - T68 - (T41_0 - T30_0);
        stb_real_t T62 = (T51_0 - T48_0 - (T53 + T56)) * STB_KP707106781;
        stb_real_t T75 = (T51_0 - T48_0 + T53 + T56) * STB_KP707106781;
        output[7 * count].real = T43_0 - T46_0 - T62;
        output[5 * count].imag = T76 - T77 - T75;
        output[3 * count].real = T43_0 - T46_0 + T62;
        output[count].imag = T75 + T76 - T77;
        stb_real_t T58 = (T48_0 + T51_0 + T53 - T56) * STB_KP707106781;
        stb_real_t T79 = (T53 - T56 - (T48_0 + T51_0)) * STB_KP707106781;
        output[5 * count].real = T43_0 + T46_0 - T58;
        output[7 * count].imag = T77 + T76 - T79;
        output[count].real = T43_0 + T46_0 + T58;
        output[3 * count].imag = T79 + T77 + T76;
        ++m;
        ++output;
        tw += 7;
    }
}

void stb_radix_2_idit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T1 = out[0].real;
    stb_real_t T2 = out[count].real;
    out[count].real = out[0].real - T2;
    out[0].real = T1 + T2;
    stb_real_t T3 = out[0].imag;
    stb_real_t T4 = out[count].imag;
    out[count].imag = T3 - T4;
    out[0].imag = T3 + T4;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T8 = output[0].imag;
        stb_real_t T3_0 = output[count].real;
        stb_real_t T5 = output[count].imag;
        stb_real_t T4_0 = tw[0].imag;
        stb_real_t T6 = tw[0].real * T3_0 - T4_0 * T5;
        stb_real_t T7 = tw[0].real * T5 + T4_0 * T3_0;
        output[count].real = output[0].real - T6;
        output[count].imag = T8 - T7;
        output[0].real = T1_0 + T6;
        output[0].imag = T7 + T8;
        ++m;
        ++output;
        ++tw;
    }
}

void stb_radix_3_idit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T1 = out[0].real;
    stb_real_t T6 = out[0].imag;
    stb_real_t T2 = out[count].real;
    stb_real_t T3 = out[2 * count].real;
    stb_real_t T5 = (T2 - T3) * STB_KP866025403;
    stb_real_t T7 = out[count].imag;
    stb_real_t T8 = out[2 * count].imag;
    stb_real_t T12 = (T8 - T7) * STB_KP866025403;
    out[0].real = out[0].real + T2 + T3;
    out[0].imag = T6 + T7 + T8;
    stb_real_t T10 = T6 - STB_KP5 * (T7 + T8);
    out[count].imag = T5 + T10;
    out[2 * count].imag = T10 - T5;
    stb_real_t T11 = T1 - STB_KP5 * (T2 + T3);
    out[2 * count].real = T11 - T12;
    out[count].real = T11 + T12;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T18 = output[0].imag;
        stb_real_t T3_0 = output[count].real;
        stb_real_t T5_0 = output[count].imag;
        stb_real_t T4_0 = tw[0].imag;
        stb_real_t T6_0 = tw[0].real * T3_0 - T4_0 * T5_0;
        stb_real_t T15 = tw[0].real * T5_0 + T4_0 * T3_0;
        stb_real_t T8_0 = output[2 * count].real;
        stb_real_t T10_0 = output[2 * count].imag;
        stb_real_t T7_0 = tw[1].real;
        stb_real_t T9_0 = tw[1].imag;
        stb_real_t T11_0 = T7_0 * T8_0 - T9_0 * T10_0;
        stb_real_t T14 = T7_0 * T10_0 + T9_0 * T8_0;
        stb_real_t T12_0 = T6_0 + T11_0;
        stb_real_t T17 = T15 + T14;
        output[0].real = output[0].real + T12_0;
        output[0].imag = T17 + T18;
        stb_real_t T16 = (T14 - T15) * STB_KP866025403;
        output[2 * count].real = T1_0 - STB_KP5 * T12_0 - T16;
        output[count].real = T1_0 - STB_KP5 * T12_0 + T16;
        stb_real_t T20 = T18 - STB_KP5 * T17;
        output[count].imag = (T6_0 - T11_0) * STB_KP866025403 + T20;
        output[2 * count].imag = T20 - (T6_0 - T11_0) * STB_KP866025403;
        ++m;
        ++output;
        tw += 2;
    }
}

void stb_radix_4_idit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T2 = out[2 * count].real;
    stb_real_t T3 = out[0].real + T2;
    stb_real_t T11 = out[0].real - T2;
    stb_real_t T8 = out[0].imag;
    stb_real_t T9 = out[2 * count].imag;
    stb_real_t T4 = out[count].real;
    stb_real_t T5 = out[3 * count].real;
    stb_real_t T12 = out[count].imag;
    stb_real_t T13 = out[3 * count].imag;
    out[2 * count].real = T3 - (T4 + T5);
    out[2 * count].imag = T8 + T9 - (T12 + T13);
    out[0].real = T3 + T4 + T5;
    out[0].imag = T8 + T9 + T12 + T13;
    out[count].imag = T4 - T5 + T8 - T9;
    out[count].real = T11 - (T12 - T13);
    out[3 * count].imag = T8 - T9 - (T4 - T5);
    out[3 * count].real = T11 + T12 - T13;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T25 = output[0].imag;
        stb_real_t T3_0 = output[2 * count].real;
        stb_real_t T5_0 = output[2 * count].imag;
        stb_real_t T2_0 = tw[1].real;
        stb_real_t T4_0 = tw[1].imag;
        stb_real_t T6_0 = T2_0 * T3_0 - T4_0 * T5_0;
        stb_real_t T24 = T2_0 * T5_0 + T4_0 * T3_0;
        stb_real_t T9_0 = output[count].real;
        stb_real_t T11_0 = output[count].imag;
        stb_real_t T10_0 = tw[0].imag;
        stb_real_t T12_0 = tw[0].real * T9_0 - T10_0 * T11_0;
        stb_real_t T20 = tw[0].real * T11_0 + T10_0 * T9_0;
        stb_real_t T14_0 = output[3 * count].real;
        stb_real_t T16_0 = output[3 * count].imag;
        stb_real_t T13_0 = tw[2].real;
        stb_real_t T15_0 = tw[2].imag;
        stb_real_t T17 = T13_0 * T14_0 - T15_0 * T16_0;
        stb_real_t T21 = T13_0 * T16_0 + T15_0 * T14_0;
        stb_real_t T18 = T12_0 + T17;
        output[2 * count].real = output[0].real + T6_0 - T18;
        output[0].real = T1_0 + T6_0 + T18;
        output[0].imag = T20 + T21 + T24 + T25;
        output[2 * count].imag = T24 + T25 - (T20 + T21);
        output[count].real = T1_0 - T6_0 - (T20 - T21);
        output[3 * count].real = T1_0 - T6_0 + T20 - T21;
        output[count].imag = T12_0 - T17 + T25 - T24;
        output[3 * count].imag = T25 - T24 - (T12_0 - T17);
        ++m;
        ++output;
        tw += 3;
    }
}

void stb_radix_5_idit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T1 = out[0].real;
    stb_real_t T27 = out[0].imag;
    stb_real_t T2 = out[count].real;
    stb_real_t T3 = out[4 * count].real;
    stb_real_t T5 = out[2 * count].real;
    stb_real_t T6 = out[3 * count].real;
    stb_real_t T8 = T2 + T3 + T5 + T6;
    stb_real_t T10 = (T2 + T3 - (T5 + T6)) * STB_KP559016994;
    stb_real_t T12 = out[count].imag;
    stb_real_t T13 = out[4 * count].imag;
    stb_real_t T15 = out[2 * count].imag;
    stb_real_t T16 = out[3 * count].imag;
    stb_real_t T28 = T12 + T13 + T15 + T16;
    stb_real_t T26 = (T12 + T13 - (T15 + T16)) * STB_KP559016994;
    out[0].real = out[0].real + T8;
    out[0].imag = T27 + T28;
    stb_real_t T18 = STB_KP587785252 * (T12 - T13) - STB_KP951056516 * (T15 - T16);
    stb_real_t T20 = STB_KP587785252 * (T15 - T16) + STB_KP951056516 * (T12 - T13);
    stb_real_t T19 = T10 + T1 - STB_KP25 * T8;
    out[2 * count].real = T1 - STB_KP25 * T8 - T10 - T18;
    out[4 * count].real = T19 + T20;
    out[3 * count].real = T1 - STB_KP25 * T8 - T10 + T18;
    out[count].real = T19 - T20;
    stb_real_t T31 = STB_KP587785252 * (T2 - T3) - STB_KP951056516 * (T5 - T6);
    stb_real_t T30 = T26 + T27 - STB_KP25 * T28;
    stb_real_t T32 = T27 - STB_KP25 * T28 - T26;
    out[count].imag = STB_KP587785252 * (T5 - T6) + STB_KP951056516 * (T2 - T3) + T30;
    out[3 * count].imag = T32 - T31;
    out[4 * count].imag = T30 - (STB_KP587785252 * (T5 - T6) + STB_KP951056516 * (T2 - T3));
    out[2 * count].imag = T31 + T32;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T40 = output[0].imag;
        stb_real_t T3_0 = output[count].real;
        stb_real_t T5_0 = output[count].imag;
        stb_real_t T4_0 = tw[0].imag;
        stb_real_t T19_0 = output[3 * count].real;
        stb_real_t T21_0 = output[3 * count].imag;
        stb_real_t T18_0 = tw[2].real;
        stb_real_t T20_0 = tw[2].imag;
        stb_real_t T22_0 = T18_0 * T19_0 - T20_0 * T21_0;
        stb_real_t T32_0 = T18_0 * T21_0 + T20_0 * T19_0;
        stb_real_t T8_0 = output[4 * count].real;
        stb_real_t T10_0 = output[4 * count].imag;
        stb_real_t T7_0 = tw[3].real;
        stb_real_t T9_0 = tw[3].imag;
        stb_real_t T11_0 = T7_0 * T8_0 - T9_0 * T10_0;
        stb_real_t T29_0 = T7_0 * T10_0 + T9_0 * T8_0;
        stb_real_t T14_0 = output[2 * count].real;
        stb_real_t T16_0 = output[2 * count].imag;
        stb_real_t T13_0 = tw[1].real;
        stb_real_t T15_0 = tw[1].imag;
        stb_real_t T30_0 = tw[0].real * T5_0 + T4_0 * T3_0 - T29_0;
        stb_real_t T33 = T13_0 * T16_0 + T15_0 * T14_0 - T32_0;
        stb_real_t T42 = T13_0 * T14_0 - T15_0 * T16_0 - T22_0;
        stb_real_t T41 = tw[0].real * T3_0 - T4_0 * T5_0 - T11_0;
        stb_real_t T37 = tw[0].real * T5_0 + T4_0 * T3_0 + T29_0;
        stb_real_t T38 = T13_0 * T16_0 + T15_0 * T14_0 + T32_0;
        stb_real_t T39 = T37 + T38;
        stb_real_t T12_0 = tw[0].real * T3_0 - T4_0 * T5_0 + T11_0;
        stb_real_t T23_0 = T13_0 * T14_0 - T15_0 * T16_0 + T22_0;
        stb_real_t T24_0 = T12_0 + T23_0;
        output[0].real = output[0].real + T24_0;
        output[0].imag = T39 + T40;
        stb_real_t T34 = STB_KP587785252 * T30_0 - STB_KP951056516 * T33;
        stb_real_t T36 = STB_KP587785252 * T33 + STB_KP951056516 * T30_0;
        stb_real_t T27_0 = T1_0 - STB_KP25 * T24_0 - (T12_0 - T23_0) * STB_KP559016994;
        stb_real_t T35 = (T12_0 - T23_0) * STB_KP559016994 + T1_0 - STB_KP25 * T24_0;
        output[2 * count].real = T27_0 - T34;
        output[4 * count].real = T35 + T36;
        output[3 * count].real = T27_0 + T34;
        output[count].real = T35 - T36;
        stb_real_t T47 = STB_KP587785252 * T41 - STB_KP951056516 * T42;
        stb_real_t T46 = (T37 - T38) * STB_KP559016994 + T40 - STB_KP25 * T39;
        stb_real_t T48 = T40 - STB_KP25 * T39 - (T37 - T38) * STB_KP559016994;
        output[count].imag = STB_KP587785252 * T42 + STB_KP951056516 * T41 + T46;
        output[3 * count].imag = T48 - T47;
        output[4 * count].imag = T46 - (STB_KP587785252 * T42 + STB_KP951056516 * T41);
        output[2 * count].imag = T47 + T48;
        ++m;
        ++output;
        tw += 4;
    }
}

void stb_radix_6_idit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T2 = out[3 * count].real;
    stb_real_t T3 = out[0].real - T2;
    stb_real_t T11 = out[0].real + T2;
    stb_real_t T24 = out[0].imag;
    stb_real_t T25 = out[3 * count].imag;
    stb_real_t T4 = out[2 * count].real;
    stb_real_t T5 = out[5 * count].real;
    stb_real_t T7 = out[4 * count].real;
    stb_real_t T8 = out[count].real;
    stb_real_t T10 = T4 - T5 + T7 - T8;
    stb_real_t T14 = T4 + T5 + T7 + T8;
    stb_real_t T16 = out[4 * count].imag;
    stb_real_t T17 = out[count].imag;
    stb_real_t T19 = out[2 * count].imag;
    stb_real_t T20 = out[5 * count].imag;
    stb_real_t T27 = T19 - T20 + T16 - T17;
    stb_real_t T32 = T19 + T20 + T16 + T17;
    out[3 * count].real = T3 + T10;
    out[3 * count].imag = T24 - T25 + T27;
    out[0].real = T11 + T14;
    out[0].imag = T24 + T25 + T32;
    stb_real_t T22 = (T16 - T17 - (T19 - T20)) * STB_KP866025403;
    out[5 * count].real = T3 - STB_KP5 * T10 - T22;
    out[count].real = T3 - STB_KP5 * T10 + T22;
    stb_real_t T23 = (T4 - T5 - (T7 - T8)) * STB_KP866025403;
    stb_real_t T28 = T24 - T25 - STB_KP5 * T27;
    out[count].imag = T23 + T28;
    out[5 * count].imag = T28 - T23;
    stb_real_t T33 = T24 + T25 - STB_KP5 * T32;
    stb_real_t T34 = (T4 + T5 - (T7 + T8)) * STB_KP866025403;
    out[2 * count].imag = T33 - T34;
    out[4 * count].imag = T34 + T33;
    stb_real_t T36 = (T16 + T17 - (T19 + T20)) * STB_KP866025403;
    out[2 * count].real = T11 - STB_KP5 * T14 - T36;
    out[4 * count].real = T11 - STB_KP5 * T14 + T36;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T49 = output[0].imag;
        stb_real_t T3_0 = output[3 * count].real;
        stb_real_t T5_0 = output[3 * count].imag;
        stb_real_t T2_0 = tw[2].real;
        stb_real_t T4_0 = tw[2].imag;
        stb_real_t T6_0 = T2_0 * T3_0 - T4_0 * T5_0;
        stb_real_t T54 = T49 - (T2_0 * T5_0 + T4_0 * T3_0);
        stb_real_t T31_0 = output[0].real + T6_0;
        stb_real_t T50 = T2_0 * T5_0 + T4_0 * T3_0 + T49;
        stb_real_t T20_0 = output[4 * count].real;
        stb_real_t T22_0 = output[4 * count].imag;
        stb_real_t T19_0 = tw[3].real;
        stb_real_t T21_0 = tw[3].imag;
        stb_real_t T25_0 = output[count].real;
        stb_real_t T27_0 = output[count].imag;
        stb_real_t T26_0 = tw[0].imag;
        stb_real_t T28_0 = tw[0].real * T25_0 - T26_0 * T27_0;
        stb_real_t T37 = tw[0].real * T27_0 + T26_0 * T25_0;
        stb_real_t T29_0 = T19_0 * T20_0 - T21_0 * T22_0 - T28_0;
        stb_real_t T44 = T19_0 * T22_0 + T21_0 * T20_0 + T37;
        stb_real_t T33_0 = T19_0 * T20_0 - T21_0 * T22_0 + T28_0;
        stb_real_t T9_0 = output[2 * count].real;
        stb_real_t T11_0 = output[2 * count].imag;
        stb_real_t T8_0 = tw[1].real;
        stb_real_t T10_0 = tw[1].imag;
        stb_real_t T14_0 = output[5 * count].real;
        stb_real_t T16_0 = output[5 * count].imag;
        stb_real_t T13_0 = tw[4].real;
        stb_real_t T15_0 = tw[4].imag;
        stb_real_t T17_0 = T13_0 * T14_0 - T15_0 * T16_0;
        stb_real_t T40 = T13_0 * T16_0 + T15_0 * T14_0;
        stb_real_t T45 = T8_0 * T11_0 + T10_0 * T9_0 + T40;
        stb_real_t T32_0 = T8_0 * T9_0 - T10_0 * T11_0 + T17_0;
        stb_real_t T41 = T8_0 * T11_0 + T10_0 * T9_0 - T40;
        stb_real_t T42 = (T19_0 * T22_0 + T21_0 * T20_0 - T37 - T41) * STB_KP866025403;
        stb_real_t T30_0 = T8_0 * T9_0 - T10_0 * T11_0 - T17_0 + T29_0;
        stb_real_t T35_0 = output[0].real - T6_0 - STB_KP5 * T30_0;
        output[3 * count].real = output[0].real - T6_0 + T30_0;
        output[count].real = T35_0 + T42;
        output[5 * count].real = T35_0 - T42;
        stb_real_t T53 = (T8_0 * T9_0 - T10_0 * T11_0 - T17_0 - T29_0) * STB_KP866025403;
        stb_real_t T55 = T41 + T19_0 * T22_0 + T21_0 * T20_0 - T37;
        stb_real_t T56 = T54 - STB_KP5 * T55;
        output[count].imag = T53 + T56;
        output[3 * count].imag = T55 + T54;
        output[5 * count].imag = T56 - T53;
        stb_real_t T46 = (T44 - T45) * STB_KP866025403;
        stb_real_t T43 = T31_0 - STB_KP5 * (T32_0 + T33_0);
        output[0].real = T31_0 + T32_0 + T33_0;
        output[4 * count].real = T43 + T46;
        output[2 * count].real = T43 - T46;
        stb_real_t T52 = (T32_0 - T33_0) * STB_KP866025403;
        stb_real_t T51 = T50 - STB_KP5 * (T45 + T44);
        output[0].imag = T45 + T44 + T50;
        output[4 * count].imag = T52 + T51;
        output[2 * count].imag = T51 - T52;
        ++m;
        ++output;
        tw += 5;
    }
}

void stb_radix_7_idit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T1 = out[0].real;
    stb_real_t T15 = out[0].imag;
    stb_real_t T2 = out[count].real;
    stb_real_t T3 = out[6 * count].real;
    stb_real_t T22 = out[count].imag;
    stb_real_t T23 = out[6 * count].imag;
    stb_real_t T5 = out[2 * count].real;
    stb_real_t T6 = out[5 * count].real;
    stb_real_t T16 = out[2 * count].imag;
    stb_real_t T17 = out[5 * count].imag;
    stb_real_t T8 = out[3 * count].real;
    stb_real_t T9 = out[4 * count].real;
    stb_real_t T19 = out[3 * count].imag;
    stb_real_t T20 = out[4 * count].imag;
    out[0].real = out[0].real + T2 + T3 + T5 + T6 + T8 + T9;
    out[0].imag = T15 + T22 + T23 + T16 + T17 + T19 + T20;
    stb_real_t T25 = STB_KP623489801 * (T16 + T17) + T15 - (STB_KP900968867 * (T22 + T23) + STB_KP222520933 * (T19 + T20));
    out[3 * count].imag = STB_KP974927912 * (T8 - T9)
        + STB_KP433883739 * (T2 - T3)
        - STB_KP781831482 * (T5 - T6)
        + T25;
    out[4 * count].imag = T25
        - (STB_KP974927912 * (T8 - T9)
              + STB_KP433883739 * (T2 - T3)
              - STB_KP781831482 * (T5 - T6));
    stb_real_t T38 = STB_KP974927912 * (T20 - T19) + STB_KP433883739 * (T23 - T22) - STB_KP781831482 * (T17 - T16);
    out[4 * count].real = STB_KP623489801 * (T5 + T6)
        + T1
        - (STB_KP900968867 * (T2 + T3)
              + STB_KP222520933 * (T8 + T9))
        - T38;
    out[3 * count].real = STB_KP623489801 * (T5 + T6)
        + T1
        - (STB_KP900968867 * (T2 + T3)
              + STB_KP222520933 * (T8 + T9))
        + T38;
    stb_real_t T27 = STB_KP623489801 * (T22 + T23) + T15 - (STB_KP222520933 * (T16 + T17) + STB_KP900968867 * (T19 + T20));
    out[count].imag = STB_KP433883739 * (T8 - T9)
        + STB_KP781831482 * (T2 - T3)
        + STB_KP974927912 * (T5 - T6)
        + T27;
    out[6 * count].imag = T27
        - (STB_KP433883739 * (T8 - T9)
              + STB_KP781831482 * (T2 - T3)
              + STB_KP974927912 * (T5 - T6));
    stb_real_t T36 = STB_KP433883739 * (T20 - T19) + STB_KP781831482 * (T23 - T22) + STB_KP974927912 * (T17 - T16);
    out[6 * count].real = STB_KP623489801 * (T2 + T3)
        + T1
        - (STB_KP222520933 * (T5 + T6)
              + STB_KP900968867 * (T8 + T9))
        - T36;
    out[count].real = STB_KP623489801 * (T2 + T3)
        + T1
        - (STB_KP222520933 * (T5 + T6)
              + STB_KP900968867 * (T8 + T9))
        + T36;
    stb_real_t T28 = STB_KP974927912 * (T2 - T3) - STB_KP781831482 * (T8 - T9) - STB_KP433883739 * (T5 - T6);
    stb_real_t T29 = STB_KP623489801 * (T19 + T20) + T15 - (STB_KP222520933 * (T22 + T23) + STB_KP900968867 * (T16 + T17));
    out[2 * count].imag = T28 + T29;
    out[5 * count].imag = T29 - T28;
    stb_real_t T34 = STB_KP974927912 * (T23 - T22) - STB_KP781831482 * (T20 - T19) - STB_KP433883739 * (T17 - T16);
    stb_real_t T30 = STB_KP623489801 * (T8 + T9) + T1 - (STB_KP222520933 * (T2 + T3) + STB_KP900968867 * (T5 + T6));
    out[5 * count].real = T30 - T34;
    out[2 * count].real = T30 + T34;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T1_0 = output[0].real;
        stb_real_t T53 = output[0].imag;
        stb_real_t T3_0 = output[count].real;
        stb_real_t T5_0 = output[count].imag;
        stb_real_t T4_0 = tw[0].imag;
        stb_real_t T8_0 = output[6 * count].real;
        stb_real_t T10_0 = output[6 * count].imag;
        stb_real_t T7_0 = tw[5].real;
        stb_real_t T9_0 = tw[5].imag;
        stb_real_t T11_0 = T7_0 * T8_0 - T9_0 * T10_0;
        stb_real_t T12_0 = tw[0].real * T3_0 - T4_0 * T5_0 + T11_0;
        stb_real_t T54 = tw[0].real * T3_0 - T4_0 * T5_0 - T11_0;
        stb_real_t T38_0 = T7_0 * T10_0 + T9_0 * T8_0 - (tw[0].real * T5_0 + T4_0 * T3_0);
        stb_real_t T50 = tw[0].real * T5_0 + T4_0 * T3_0 + T7_0 * T10_0 + T9_0 * T8_0;
        stb_real_t T14_0 = output[2 * count].real;
        stb_real_t T16_0 = output[2 * count].imag;
        stb_real_t T13_0 = tw[1].real;
        stb_real_t T15_0 = tw[1].imag;
        stb_real_t T19_0 = output[5 * count].real;
        stb_real_t T21_0 = output[5 * count].imag;
        stb_real_t T18_0 = tw[4].real;
        stb_real_t T20_0 = tw[4].imag;
        stb_real_t T22_0 = T18_0 * T19_0 - T20_0 * T21_0;
        stb_real_t T23_0 = T13_0 * T14_0 - T15_0 * T16_0 + T22_0;
        stb_real_t T55 = T13_0 * T14_0 - T15_0 * T16_0 - T22_0;
        stb_real_t T44 = T18_0 * T21_0 + T20_0 * T19_0 - (T13_0 * T16_0 + T15_0 * T14_0);
        stb_real_t T51 = T13_0 * T16_0 + T15_0 * T14_0 + T18_0 * T21_0 + T20_0 * T19_0;
        stb_real_t T25_0 = output[3 * count].real;
        stb_real_t T27_0 = output[3 * count].imag;
        stb_real_t T24_0 = tw[2].real;
        stb_real_t T26_0 = tw[2].imag;
        stb_real_t T30_0 = output[4 * count].real;
        stb_real_t T32_0 = output[4 * count].imag;
        stb_real_t T29_0 = tw[3].real;
        stb_real_t T31_0 = tw[3].imag;
        stb_real_t T33_0 = T29_0 * T30_0 - T31_0 * T32_0;
        stb_real_t T34_0 = T24_0 * T25_0 - T26_0 * T27_0 + T33_0;
        stb_real_t T56 = T24_0 * T25_0 - T26_0 * T27_0 - T33_0;
        stb_real_t T41 = T29_0 * T32_0 + T31_0 * T30_0 - (T24_0 * T27_0 + T26_0 * T25_0);
        stb_real_t T52 = T24_0 * T27_0 + T26_0 * T25_0 + T29_0 * T32_0 + T31_0 * T30_0;
        output[0].real = output[0].real + T12_0 + T23_0 + T34_0;
        output[0].imag = T50 + T51 + T52 + T53;
        stb_real_t T45 = STB_KP974927912 * T38_0 - STB_KP781831482 * T41 - STB_KP433883739 * T44;
        output[5 * count].real = STB_KP623489801 * T34_0
            + T1_0
            - (STB_KP222520933 * T12_0
                  + STB_KP900968867 * T23_0)
            - T45;
        output[2 * count].real = STB_KP623489801 * T34_0
            + T1_0
            - (STB_KP222520933 * T12_0
                  + STB_KP900968867 * T23_0)
            + T45;
        stb_real_t T60 = STB_KP623489801 * T52 + T53 - (STB_KP222520933 * T50 + STB_KP900968867 * T51);
        output[2 * count].imag = STB_KP974927912 * T54 - STB_KP781831482 * T56 - STB_KP433883739 * T55 + T60;
        output[5 * count].imag = T60 - (STB_KP974927912 * T54 - STB_KP781831482 * T56 - STB_KP433883739 * T55);
        stb_real_t T47 = STB_KP433883739 * T41 + STB_KP781831482 * T38_0 + STB_KP974927912 * T44;
        output[6 * count].real = STB_KP623489801 * T12_0
            + T1_0
            - (STB_KP222520933 * T23_0
                  + STB_KP900968867 * T34_0)
            - T47;
        output[count].real = STB_KP623489801 * T12_0
            + T1_0
            - (STB_KP222520933 * T23_0
                  + STB_KP900968867 * T34_0)
            + T47;
        stb_real_t T58 = STB_KP623489801 * T50 + T53 - (STB_KP222520933 * T51 + STB_KP900968867 * T52);
        output[count].imag = STB_KP433883739 * T56 + STB_KP781831482 * T54 + STB_KP974927912 * T55 + T58;
        output[6 * count].imag = T58 - (STB_KP433883739 * T56 + STB_KP781831482 * T54 + STB_KP974927912 * T55);
        stb_real_t T49 = STB_KP974927912 * T41 + STB_KP433883739 * T38_0 - STB_KP781831482 * T44;
        output[4 * count].real = STB_KP623489801 * T23_0
            + T1_0
            - (STB_KP900968867 * T12_0
                  + STB_KP222520933 * T34_0)
            - T49;
        output[3 * count].real = STB_KP623489801 * T23_0
            + T1_0
            - (STB_KP900968867 * T12_0
                  + STB_KP222520933 * T34_0)
            + T49;
        stb_real_t T62 = STB_KP623489801 * T51 + T53 - (STB_KP900968867 * T50 + STB_KP222520933 * T52);
        output[3 * count].imag = STB_KP974927912 * T56 + STB_KP433883739 * T54 - STB_KP781831482 * T55 + T62;
        output[4 * count].imag = T62
            - (STB_KP974927912 * T56
                  + STB_KP433883739 * T54
                  - STB_KP781831482 * T55);
        ++m;
        ++output;
        tw += 6;
    }
}

void stb_radix_8_idit(const cmplx* tw, cmplx* out, int count)
{
    int m;
    cmplx* output;
    stb_real_t T2 = out[4 * count].real;
    stb_real_t T3 = out[0].real + T2;
    stb_real_t T37 = out[0].real - T2;
    stb_real_t T16 = out[0].imag;
    stb_real_t T17 = out[4 * count].imag;
    stb_real_t T4 = out[2 * count].real;
    stb_real_t T5 = out[6 * count].real;
    stb_real_t T19 = out[2 * count].imag;
    stb_real_t T20 = out[6 * count].imag;
    stb_real_t T11 = out[7 * count].real;
    stb_real_t T12 = out[3 * count].real;
    stb_real_t T32 = out[7 * count].imag;
    stb_real_t T33 = out[3 * count].imag;
    stb_real_t T35 = T11 - T12 + T32 - T33;
    stb_real_t T43 = T32 - T33 - (T11 - T12);
    stb_real_t T8 = out[count].real;
    stb_real_t T9 = out[5 * count].real;
    stb_real_t T27 = out[count].imag;
    stb_real_t T28 = out[5 * count].imag;
    stb_real_t T30 = T8 - T9 - (T27 - T28);
    stb_real_t T42 = T8 - T9 + T27 - T28;
    stb_real_t T7 = T3 + T4 + T5;
    stb_real_t T14 = T8 + T9 + T11 + T12;
    out[4 * count].real = T7 - T14;
    out[0].real = T7 + T14;
    out[4 * count].imag = T16 + T17 + T19 + T20 - (T27 + T28 + T32 + T33);
    out[0].imag = T16 + T17 + T19 + T20 + T27 + T28 + T32 + T33;
    stb_real_t T15 = T8 + T9 - (T11 + T12);
    stb_real_t T22 = T16 + T17 - (T19 + T20);
    out[2 * count].imag = T15 + T22;
    out[6 * count].imag = T22 - T15;
    stb_real_t T51 = T3 - (T4 + T5);
    stb_real_t T52 = T32 + T33 - (T27 + T28);
    out[6 * count].real = T51 - T52;
    out[2 * count].real = T51 + T52;
    stb_real_t T36 = (T30 - T35) * STB_KP707106781;
    out[7 * count].imag = T16 - T17 - (T4 - T5) - T36;
    out[3 * count].imag = T16 - T17 - (T4 - T5) + T36;
    stb_real_t T46 = (T43 - T42) * STB_KP707106781;
    out[7 * count].real = T37 + T19 - T20 - T46;
    out[3 * count].real = T37 + T19 - T20 + T46;
    stb_real_t T39 = T37 - (T19 - T20);
    stb_real_t T40 = (T30 + T35) * STB_KP707106781;
    out[5 * count].real = T39 - T40;
    out[count].real = T39 + T40;
    stb_real_t T41 = T4 - T5 + T16 - T17;
    stb_real_t T44 = (T42 + T43) * STB_KP707106781;
    out[5 * count].imag = T41 - T44;
    out[count].imag = T41 + T44;
    output = out + 1;
    m = 1;
    while (m < count) {
        stb_real_t T70 = output[0].imag;
        stb_real_t T3_0 = output[4 * count].real;
        stb_real_t T5_0 = output[4 * count].imag;
        stb_real_t T2_0 = tw[3].real;
        stb_real_t T4_0 = tw[3].imag;
        stb_real_t T6_0 = T2_0 * T3_0 - T4_0 * T5_0;
        stb_real_t T7_0 = output[0].real + T6_0;
        stb_real_t T77 = T70 - (T2_0 * T5_0 + T4_0 * T3_0);
        stb_real_t T43_0 = output[0].real - T6_0;
        stb_real_t T71 = T2_0 * T5_0 + T4_0 * T3_0 + T70;
        stb_real_t T32_0 = output[7 * count].real;
        stb_real_t T34_0 = output[7 * count].imag;
        stb_real_t T31_0 = tw[6].real;
        stb_real_t T33_0 = tw[6].imag;
        stb_real_t T37_0 = output[3 * count].real;
        stb_real_t T39_0 = output[3 * count].imag;
        stb_real_t T36_0 = tw[2].real;
        stb_real_t T38_0 = tw[2].imag;
        stb_real_t T40_0 = T36_0 * T37_0 - T38_0 * T39_0;
        stb_real_t T55 = T36_0 * T39_0 + T38_0 * T37_0;
        stb_real_t T41_0 = T31_0 * T32_0 - T33_0 * T34_0 + T40_0;
        stb_real_t T64 = T31_0 * T34_0 + T33_0 * T32_0 + T55;
        stb_real_t T53 = T31_0 * T32_0 - T33_0 * T34_0 - T40_0;
        stb_real_t T56 = T31_0 * T34_0 + T33_0 * T32_0 - T55;
        stb_real_t T9_0 = output[2 * count].real;
        stb_real_t T11_0 = output[2 * count].imag;
        stb_real_t T8_0 = tw[1].real;
        stb_real_t T10_0 = tw[1].imag;
        stb_real_t T14_0 = output[6 * count].real;
        stb_real_t T16_0 = output[6 * count].imag;
        stb_real_t T13_0 = tw[5].real;
        stb_real_t T15_0 = tw[5].imag;
        stb_real_t T17_0 = T13_0 * T14_0 - T15_0 * T16_0;
        stb_real_t T45_0 = T13_0 * T16_0 + T15_0 * T14_0;
        stb_real_t T18_0 = T8_0 * T9_0 - T10_0 * T11_0 + T17_0;
        stb_real_t T76 = T8_0 * T9_0 - T10_0 * T11_0 - T17_0;
        stb_real_t T46_0 = T8_0 * T11_0 + T10_0 * T9_0 - T45_0;
        stb_real_t T68 = T8_0 * T11_0 + T10_0 * T9_0 + T45_0;
        stb_real_t T21_0 = output[count].real;
        stb_real_t T23_0 = output[count].imag;
        stb_real_t T22_0 = tw[0].imag;
        stb_real_t T26_0 = output[5 * count].real;
        stb_real_t T28_0 = output[5 * count].imag;
        stb_real_t T25_0 = tw[4].real;
        stb_real_t T27_0 = tw[4].imag;
        stb_real_t T29_0 = T25_0 * T26_0 - T27_0 * T28_0;
        stb_real_t T50_0 = T25_0 * T28_0 + T27_0 * T26_0;
        stb_real_t T30_0 = tw[0].real * T21_0 - T22_0 * T23_0 + T29_0;
        stb_real_t T65 = tw[0].real * T23_0 + T22_0 * T21_0 + T50_0;
        stb_real_t T48_0 = tw[0].real * T21_0 - T22_0 * T23_0 - T29_0;
        stb_real_t T51_0 = tw[0].real * T23_0 + T22_0 * T21_0 - T50_0;
        stb_real_t T19_0 = T7_0 + T18_0;
        stb_real_t T42_0 = T30_0 + T41_0;
        output[4 * count].real = T19_0 - T42_0;
        output[0].real = T19_0 + T42_0;
        output[0].imag = T65 + T64 + T68 + T71;
        output[4 * count].imag = T68 + T71 - (T65 + T64);
        output[6 * count].real = T7_0 - T18_0 - (T64 - T65);
        output[2 * count].real = T7_0 - T18_0 + T64 - T65;
        output[2 * count].imag = T30_0 - T41_0 + T71 - T68;
        output[6 * count].imag = T71 - T68 - (T30_0 - T41_0);
        stb_real_t T62 = (T56 - T53 - (T48_0 + T51_0)) * STB_KP707106781;
        stb_real_t T75 = (T48_0 + T51_0 + T56 - T53) * STB_KP707106781;
        output[7 * count].real = T43_0 + T46_0 - T62;
        output[5 * count].imag = T76 + T77 - T75;
        output[3 * count].real = T43_0 + T46_0 + T62;
        output[count].imag = T75 + T76 + T77;
        stb_real_t T58 = (T48_0 - T51_0 + T53 + T56) * STB_KP707106781;
        stb_real_t T79 = (T48_0 - T51_0 - (T53 + T56)) * STB_KP707106781;
        output[5 * count].real = T43_0 - T46_0 - T58;
        output[7 * count].imag = T77 - T76 - T79;
        output[count].real = T43_0 - T46_0 + T58;
        output[3 * count].imag = T79 + T77 - T76;
        ++m;
        ++output;
        tw += 7;
    }
}

void stb_recursive_mixed_radix_dit(const stb_fft_plan* plan, int stage, cmplx* in, cmplx* out, int sign)
{
    const int radix = plan->stages.radix[stage];
    const int count = plan->stages.remainder[stage];
    cmplx* twiddles = plan->twiddles;
    const int tw_offset = plan->stages.offsets[stage];
    cmplx* tw_sequential = &plan->twiddles_ordered[tw_offset];
    if (count == 1) {
        for (int i = 0; i < radix; i++) {
            out[i] = in[i * sign];
        }
    } else {
        const int cur_sign = sign * radix;
        for (int i = 0; i < radix; ++i) {
            stb_recursive_mixed_radix_dit(plan, stage + 1, &in[sign * i], &out[count * i], cur_sign);
        }
    }
    switch (radix) {
    case 2:
        stb_radix_2_dit(tw_sequential, out, count);
        break;
    case 3:
        stb_radix_3_dit(tw_sequential, out, count);
        break;
    case 4:
        stb_radix_4_dit(tw_sequential, out, count);
        break;
    case 5:
        stb_radix_5_dit(tw_sequential, out, count);
        break;
    case 6:
        stb_radix_6_dit(tw_sequential, out, count);
        break;
    case 7:
        stb_radix_7_dit(tw_sequential, out, count);
        break;
    case 8:
        stb_radix_8_dit(tw_sequential, out, count);
        break;
    default:
        stb_general_dit(twiddles, out, count, sign, radix, plan->N, 0);
        break;
    }
}

void stb_recursive_mixed_radix_idit(const stb_fft_plan* plan, int stage, cmplx* in, cmplx* out, int sign)
{
    const int radix = plan->stages.radix[stage];
    const int count = plan->stages.remainder[stage];
    cmplx* twiddles = plan->twiddles;
    const int tw_offset = plan->stages.offsets[stage];
    cmplx* tw_sequential = &plan->twiddles_ordered[tw_offset];
    if (count == 1) {
        for (int i = 0; i < radix; i++) {
            out[i] = in[i * sign];
        }
    } else {
        const int cur_sign = sign * radix;
        for (int i = 0; i < radix; ++i) {
            stb_recursive_mixed_radix_idit(plan, stage + 1, &in[sign * i], &out[count * i],
                cur_sign);
        }
    }
    switch (radix) {
    case 2:
        stb_radix_2_idit(tw_sequential, out, count);
        break;
    case 3:
        stb_radix_3_idit(tw_sequential, out, count);
        break;
    case 4:
        stb_radix_4_idit(tw_sequential, out, count);
        break;
    case 5:
        stb_radix_5_idit(tw_sequential, out, count);
        break;
    case 6:
        stb_radix_6_idit(tw_sequential, out, count);
        break;
    case 7:
        stb_radix_7_idit(tw_sequential, out, count);
        break;
    case 8:
        stb_radix_8_idit(tw_sequential, out, count);
        break;
    default:
        stb_general_dit(twiddles, out, count, sign, radix, plan->N, 1);
        break;
    }
}

void stb_fft_exec(const stb_fft_plan* plan, cmplx* in, cmplx* out)
{
    if (in == out) {
        cmplx* buf = (cmplx*)calloc(plan->N, sizeof(cmplx));
        if (buf == NULL)
            return;
        stb_recursive_mixed_radix_dit(plan, 0, in, buf, 1);
        memcpy(out, buf, sizeof(cmplx) * plan->N);
        free(buf);
    } else {
        stb_recursive_mixed_radix_dit(plan, 0, in, out, 1);
    }
}

void stb_ifft_exec(const stb_fft_plan* plan, cmplx* in, cmplx* out)
{
    if (in == out) {
        cmplx* buf = (cmplx*)calloc(plan->N, sizeof(cmplx));
        if (buf == NULL)
            return;
        stb_recursive_mixed_radix_idit(plan, 0, in, buf, 1);
        memcpy(out, buf, sizeof(cmplx) * plan->N);
        free(buf);
    } else {
        stb_recursive_mixed_radix_idit(plan, 0, in, out, 1);
    }
}

void stb_fft_r2c_exec(stb_fft_real_plan* plan, const stb_real_t* input, cmplx* output)
{
    int n = plan->half_plan->N;
    stb_fft_exec(plan->half_plan, (cmplx*)input, plan->buffer);
    cmplx tdc = plan->buffer[0];
    output[0].real = tdc.imag + tdc.real;
    output[n].real = tdc.real - tdc.imag;
    output[0].imag = 0.0;
    output[n].imag = output[0].imag;
    for (int c = 1; c <= n / 2; ++c) {
        cmplx t = plan->buffer[n - c];
        t.imag = -t.imag;
        stb_real_t r = (plan->twiddles[c - 1].real * (plan->buffer[c].real - t.real)) - (plan->twiddles[c - 1].imag * (plan->buffer[c].imag - t.imag));
        stb_real_t i = (plan->twiddles[c - 1].real * (plan->buffer[c].imag - t.imag)) + ((plan->buffer[c].real - t.real) * plan->twiddles[c - 1].imag);
        output[c].real = (r + (t.real + plan->buffer[c].real)) * STB_KP5;
        output[c].imag = (i + (t.imag + plan->buffer[c].imag)) * STB_KP5;
        output[n - c].real = ((t.real + plan->buffer[c].real) - r) * STB_KP5;
        output[n - c].imag = (i - (t.imag + plan->buffer[c].imag)) * STB_KP5;
    }
}

void stb_fft_c2r_exec(stb_fft_real_plan* plan, const cmplx* input, stb_real_t* output)
{
    int n = plan->half_plan->N;
    plan->buffer[0].real = input[n].real + input[0].real;
    plan->buffer[0].imag = input[0].real - input[n].real;
    for (int c = 1; c <= n / 2; ++c) {
        cmplx t = input[n - c];
        t.imag = -t.imag;
        stb_real_t r = (plan->twiddles[c - 1].real * (input[c].real - t.real)) + (plan->twiddles[c - 1].imag * (input[c].imag - t.imag));
        stb_real_t i = (plan->twiddles[c - 1].real * (input[c].imag - t.imag)) - ((input[c].real - t.real) * plan->twiddles[c - 1].imag);
        plan->buffer[c].real = (r + t.real + input[c].real);
        plan->buffer[c].imag = (i + t.imag + input[c].imag);
        plan->buffer[n - c].real = (t.real + input[c].real - r);
        plan->buffer[n - c].imag = (i - t.imag - input[c].imag);
    }
    stb_ifft_exec(plan->half_plan, plan->buffer, (cmplx*)output);
}

void STB_FFT(cmplx* input, cmplx* output, int n)
{
    if (n < 2) {
        output[0] = input[0];
        return;
    }
    stb_fft_plan* plan = stb_fft_plan_dft_1d(n);
    if (plan != NULL) {
        stb_fft_exec(plan, input, output);
        free(plan);
    }
}

void STB_IFFT(cmplx* input, cmplx* output, int n)
{
    if (n < 2) {
        output[0] = input[0];
        return;
    }
    stb_fft_plan* plan = stb_fft_plan_dft_1d(n);
    if (plan != NULL) {
        stb_ifft_exec(plan, input, output);
        free(plan);
    }
}

void STB_FFT_R2C(stb_real_t* input, cmplx* output, int n)
{
    if (n < 2) {
        output[0].real = input[0];
        return;
    }
    stb_fft_real_plan* stb_fft_plan = stb_fft_real_plan_dft_1d(n);
    if (stb_fft_plan != NULL) {
        stb_fft_r2c_exec(stb_fft_plan, input, output);
        free(stb_fft_plan);
    }
}

void STB_IFFT_C2R(cmplx* input, stb_real_t* output, int n)
{
    if (n < 2) {
        output[0] = input[0].real;
        return;
    }
    stb_fft_real_plan* stb_fft_plan = stb_fft_real_plan_dft_1d(n);
    if (stb_fft_plan != NULL) {
        stb_fft_c2r_exec(stb_fft_plan, input, output);
        free(stb_fft_plan);
    }
}

#endif

// REVISION HISTORY
// v0.11 - 2018-12-28
//   - Initial versioned release.

/*
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.
In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
For more information, please refer to <http://unlicense.org/>
*/