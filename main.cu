#define _CRT_SECURE_NO_WARNINGS
#define D_SCL_SECURE_NO_WARNINGS

#pragma once
#include "online_stft_arrayfire.hpp"
#include "cpp_buffer_to_wave/WaveManager_revision\WaveManager.h"
#include "complex_double_binary_to_matlab\binary2matlab_converter.hpp"


// var. for test 

typedef af::array arr;
typedef af::dim4 dim4;
typedef af::seq seq;


int main(void)
{
    const int sfft = 512;
    const int nfreq = sfft / 2+1;
    const int soverlap = sfft / 2;
    const int sshift = sfft - soverlap;
    const int nCh = 7;

    stft_AF* pstft = new stft_AF(sfft, nCh, sshift);

    float** buffer;
    int nBlocks = wm::wread("wav_7ch.wav", buffer);

    int nbuffers = nBlocks / sshift;

    arr samples = af::constant(0.0f, nBlocks, nCh, c32);

    arr frames = af::constant(0.0f, nfreq, nCh, nbuffers, c32);

    arr processed = af::constant(0.0f, nBlocks, nCh, f32);

    
    for (int i = 0; i < nCh; i++)
    {
        samples(af::span, i) = arr(nBlocks, (&buffer[i][0]), afHost).as(c32);
    }
    
    arr frame = af::constant(0.0f, nfreq, c32);

    af::timer::start();
    for (int i = 0; i < nbuffers; i++)
    {
        int start = i * sshift;
        
        arr arrived = samples(seq(start, start + sshift - 1), af::span);
        
        frame = pstft->stft(arrived);
        frames(af::span, af::span, i) = frame;
        processed(seq(start, start + sshift - 1), af::span) = pstft->istft(frame, false);
    }
    printf("elapsed seconds per 1 frame: %g\n", af::timer::stop() / (double)nbuffers);

    float** ptr7ch = new float*[nCh];

    for (int i = 0; i < nCh; i++)
        ptr7ch[i] = (processed(af::span, i)).host<float>();

    wm::wwrite("ch7_processed.wav", ptr7ch, nBlocks, nCh);

    processed.unlock();

    delete pstft;

    delete buffer;

    for (int i = 0; i < nCh; i++)
        delete ptr7ch[i];

    delete ptr7ch;

    return 0;
}