
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define DR_WAV_IMPLEMENTATION

#include "dr_wav.h"

#define DR_MP3_IMPLEMENTATION


#include "dr_mp3.h"

#include "timing.h"


#define STB_FFT_IMPLEMENTAION

#include "stb_fft.h"

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif


void wavWrite_f32(char *filename, float *buffer, int sampleRate, uint32_t totalSampleCount, uint32_t channels)
{
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = channels;
    format.sampleRate = (drwav_uint32) sampleRate;
    format.bitsPerSample = 32;
    drwav *pWav = drwav_open_file_write(filename, &format);
    if (pWav) {
        drwav_uint64 samplesWritten = drwav_write(pWav, totalSampleCount, buffer);
        drwav_uninit(pWav);
        if (samplesWritten != totalSampleCount) {
            fprintf(stderr, "write file [%s] error.\n", filename);
            exit(1);
        }
    }
}

float *wavRead_f32(const char *filename, uint32_t *sampleRate, uint64_t *sampleCount, uint32_t *channels)
{
    drwav_uint64 totalSampleCount = 0;
    float *input = drwav_open_file_and_read_pcm_frames_f32(filename, channels, sampleRate, &totalSampleCount);
    if (input == NULL) {
        drmp3_config pConfig;
        input = drmp3_open_file_and_read_f32(filename, &pConfig, &totalSampleCount);
        if (input != NULL) {
            *channels = pConfig.outputChannels;
            *sampleRate = pConfig.outputSampleRate;
        }
    }
    if (input == NULL) {
        fprintf(stderr, "read file [%s] error.\n", filename);
        exit(1);
    }
    *sampleCount = totalSampleCount * (*channels);
    return input;
}

void splitpath(const char *path, char *drv, char *dir, char *name, char *ext)
{
    const char *end;
    const char *p;
    const char *s;
    if (path[0] && path[1] == ':') {
        if (drv) {
            *drv++ = *path++;
            *drv++ = *path++;
            *drv = '\0';
        }
    }
    else if (drv)
        *drv = '\0';
    for (end = path; *end && *end != ':';)
        end++;
    for (p = end; p > path && *--p != '\\' && *p != '/';)
        if (*p == '.') {
            end = p;
            break;
        }
    if (ext)
        for (s = end; (*ext = *s++);)
            ext++;
    for (p = end; p > path;)
        if (*--p == '\\' || *p == '/') {
            p++;
            break;
        }
    if (name) {
        for (s = p; s < end;)
            *name++ = *s++;
        *name = '\0';
    }
    if (dir) {
        for (s = path; s < p;)
            *dir++ = *s++;
        *dir = '\0';
    }
}

typedef struct
{
    int frameSize;
    int windowSize;
    int freq_size;
    int sampleRate;
    float *windowing;
    stb_fft_real_plan *realPlan;
    float *fifo;
    float *synthesis_mem;
    float *smooth_mem;
    cmplx *samples;
    float *noise_mem;
    int noise_count;
} SimpleDenoiseHandle;

void SimpleDenoise_Free(SimpleDenoiseHandle *handle)
{
    if (handle) {
        if (handle->fifo) {
            free(handle->fifo);
            handle->fifo = NULL;
        }
        if (handle->smooth_mem) {
            free(handle->smooth_mem);
            handle->smooth_mem = NULL;
        }

        if (handle->realPlan) {
            free(handle->realPlan);
            handle->realPlan = NULL;
        }

        if (handle->samples) {
            free(handle->samples);
            handle->samples = NULL;
        }

        if (handle->synthesis_mem) {
            free(handle->synthesis_mem);
            handle->synthesis_mem = NULL;
        }

        if (handle->windowing) {
            free(handle->windowing);
            handle->windowing = NULL;
        }
        if (handle->noise_mem) {
            free(handle->noise_mem);
            handle->noise_mem = NULL;
        }
    }
}

int SimpleDenoise_Init(SimpleDenoiseHandle *handle, size_t sampleRate, size_t ms)
{
    if (handle) {
        size_t nfft = (ms * sampleRate / 1000);
        nfft += nfft % 2;
        handle->frameSize = nfft;
        handle->freq_size = handle->frameSize + 1;
        handle->windowSize = nfft * 2;
        handle->sampleRate = sampleRate;
        handle->windowing = (float *) calloc(handle->frameSize, sizeof(float));
        handle->realPlan = stb_fft_real_plan_dft_1d(handle->windowSize);
        handle->fifo = (float *) calloc(handle->windowSize, sizeof(float));
        handle->samples = (cmplx *) calloc(handle->windowSize, sizeof(cmplx));
        handle->synthesis_mem = (float *) calloc(handle->frameSize, sizeof(float));
        handle->noise_mem = (float *) calloc(handle->freq_size, sizeof(float));
        handle->smooth_mem = (float *) calloc(handle->freq_size, sizeof(float));
        handle->noise_count = 0;
        if ((handle->fifo == NULL) || (handle->realPlan == NULL)
            || (handle->noise_mem == NULL)
            || (handle->smooth_mem == NULL)
            || (handle->samples == NULL)
            || (handle->synthesis_mem == NULL) || (handle->windowing == NULL)
                ) {
            SimpleDenoise_Free(handle);
            return 0;
        }
        for (size_t i = 0; i < handle->frameSize; i++) {
            double t = sin(.5 * M_PI * (i + .5) / handle->frameSize);
            handle->windowing[i] = (float) sin(.5 * M_PI * t * t);
        }
        return 1;
    }
    return 0;
}

int Simple_NoiseEstimator(SimpleDenoiseHandle *handle, const float *input, int Sampling)
{
    if (handle == NULL) {
        return -1;
    }
    float *fifo = handle->fifo;
    float *noise_mem = handle->noise_mem;
    if (Sampling == 0) {
        if (input == NULL) {
            return -1;
        }
        float norm = 1.0f / handle->windowSize;
        for (size_t i = 0; i < handle->frameSize; i++) {
            fifo[i] *= handle->windowing[i] * norm;
            fifo[handle->windowSize - 1 - i] = input[handle->frameSize - 1 - i] * handle->windowing[i] * norm;
        }
        stb_fft_r2c_exec(handle->realPlan, fifo, handle->samples);
        float mag = 0;
        for (size_t i = 0; i < handle->freq_size; i++) {
            mag =
                (handle->samples[i].real * handle->samples[i].real + handle->samples[i].imag * handle->samples[i].imag);
            noise_mem[i] += mag;
        }
        handle->noise_count++;
    }
    else {
        float noise_norm = 1.0 / max(1, handle->noise_count);
        for (size_t i = 0; i < handle->freq_size; i++) {
            noise_mem[i] *= noise_norm;
        }
    }
    return 1;
}

int SimpleDenoise_Proc(SimpleDenoiseHandle *handle, const float *input, float *output)
{
    if ((input == NULL) || (handle == NULL) || (output == NULL)) {
        return -1;
    }
    float *fifo = handle->fifo;
    float *synthesis_mem = handle->synthesis_mem;
    float norm = 1.0f / handle->windowSize;
    for (size_t i = 0; i < handle->frameSize; i++) {
        fifo[i] *= handle->windowing[i] * norm;
        fifo[handle->windowSize - 1 - i] = input[handle->frameSize - 1 - i] * handle->windowing[i] * norm;
    }
    stb_fft_r2c_exec(handle->realPlan, fifo, handle->samples);
    float smooth = 0.98;
    float inv_smooth = 1.0f - smooth;
    for (size_t i = 0; i < handle->freq_size; i++) {
        const float mag =
            (handle->samples[i].real * handle->samples[i].real + handle->samples[i].imag * handle->samples[i].imag);
        const float smooth_mag = smooth * handle->smooth_mem[i] + inv_smooth * mag;
        const float gain = smooth_mag / max(handle->noise_mem[i], smooth_mag);
        handle->samples[i].real *= gain;
        handle->samples[i].imag *= gain;
        handle->smooth_mem[i] = mag * gain * gain;
    }
    stb_fft_c2r_exec(handle->realPlan, handle->samples, fifo);
    for (size_t i = 0; i < handle->frameSize; i++) {
        output[i] = fifo[i] * handle->windowing[i] + synthesis_mem[i];
        fifo[handle->windowSize - 1 - i] *= handle->windowing[i];
    }
    memcpy(synthesis_mem, fifo + handle->frameSize, handle->frameSize * sizeof(float));
    memcpy(fifo, input, handle->frameSize * sizeof(float));
    return 1;
}

void printUsage()
{
    printf("usage:\n");
    printf("./SimpleDenoise input.wav\n");
    printf("./SimpleDenoise input.mp3\n");
    printf("or\n");
    printf("./SimpleDenoise input.wav output.wav\n");
    printf("./SimpleDenoise input.mp3 output.wav\n");
    printf("press any key to exit.\n");
    getchar();
}

void simpleDenoise(char *in_file, char *out_file)
{
    if (in_file == NULL || out_file == NULL) {
        printUsage();
        return;
    }
    uint32_t sampleRate = 0;
    uint64_t sampleCount = 0;
    uint32_t channels = 0;
    float *input = wavRead_f32(in_file, &sampleRate, &sampleCount, &channels);
    if (input) {
        size_t ms = 20;
        SimpleDenoiseHandle *handle = (SimpleDenoiseHandle *) malloc(sizeof(SimpleDenoiseHandle));
        if (handle) {
            double startTime = now();
            if (SimpleDenoise_Init(handle, sampleRate, ms) == 1) {
                uint64_t frames = (sampleCount / handle->frameSize);
                int remainingSample = (sampleCount % handle->frameSize);
                float *output = (float *) calloc(sampleCount + handle->frameSize, sizeof(float));
                if (output) {
                    float *inBuffer = input;
                    int sampling = 0;
                    for (int n = 0; n < frames; ++n) {
                        Simple_NoiseEstimator(handle, inBuffer, sampling);
                        inBuffer += handle->frameSize;
                    }
                    sampling = 1;
                    Simple_NoiseEstimator(handle, inBuffer, sampling);
                    inBuffer = input;
                    float *outBuffer = output;
                    for (int n = 0; n < frames; ++n) {
                        if (SimpleDenoise_Proc(handle, inBuffer, outBuffer) == 1)
                            outBuffer += handle->frameSize;
                        inBuffer += handle->frameSize;
                    }
                    if (remainingSample != 0) {
                        memcpy(outBuffer, handle->synthesis_mem, sizeof(float) * handle->frameSize);
                    } else {
                        float *buffer = (float *) calloc(handle->frameSize * 2, sizeof(float));
                        if (buffer) {
                            memcpy(buffer, inBuffer, sizeof(float) * remainingSample);
                            SimpleDenoise_Proc(handle, buffer, outBuffer);
                            outBuffer += handle->frameSize;
                            memcpy(outBuffer, handle->synthesis_mem, sizeof(float) * remainingSample);
                            free(buffer);
                        }
                    }
                    double time_interval = calcElapsed(startTime, now());
                    printf("time interval: %f ms\n ", (time_interval * 1000));
                    wavWrite_f32(out_file, output, sampleRate, (uint32_t) sampleCount, channels);
                    free(output);
                }
            }
            SimpleDenoise_Free(handle);
            free(handle);
        }
        free(input);
    }
}

int main(int argc, char *argv[])
{
    printf("Audio Processing\n");
    printf("blog:http://cpuimage.cnblogs.com/\n");
    printf("Audio Simple Denoise\n");
    if (argc < 2) {
        printUsage();
        return -1;
    }
    char *in_file = argv[1];
    if (argc > 2) {
        char *out_file = argv[2];
        simpleDenoise(in_file, out_file);
    }
    else {
        char drive[3];
        char dir[256];
        char fname[256];
        char ext[256];
        char out_file[1024];
        splitpath(in_file, drive, dir, fname, ext);
        sprintf(out_file, "%s%s%s_out.wav", drive, dir, fname);
        simpleDenoise(in_file, out_file);
    }
    printf("done.\n");
    return 0;
}

#ifdef __cplusplus
}
#endif
