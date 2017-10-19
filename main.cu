#define _CRT_SECURE_NO_WARNINGS
#define D_SCL_SECURE_NO_WARNINGS

#include <arrayfire.h>
#include "cpp_buffer_to_wave/WaveManager_revision\WaveManager.h"
#include "complex_double_binary_to_matlab\binary2matlab_converter.hpp"

// var. for test 

typedef af::array arr;
typedef af::dim4 dim4;
typedef af::seq seq;

class stft_AF
{
    typedef af::array arr;
    typedef af::dim4 dim4;
    typedef af::seq seq;

    int sfft;
    int nch;
    int soverlap;
    int sfreq;
    

    arr window;
    
    std::vector<double> vHanWinCoeff;
    arr hanWinCoeff;

public:

    stft_AF(const int sfft, const int nch, const int soverlap)
    {
        // make hanning window
        initialize(sfft, nch, soverlap);
    }

    void initialize(const int _sfft, const int _nch, const int _soverlap)
    {
        sfft = _sfft;
        soverlap = _soverlap;
        sfreq = _sfft / 2 + 1;
        nch = _nch;

        window = af::constant(0.0f, sfft, nch, c32);

        make_hanning_window();
    }

    void make_hanning_window()
    {
        std::cout << "make_hanning_window\n";

        vHanWinCoeff.resize(sfft);
        
        std::cout << "make_hanning_window\n";

        for (int sample = 0; sample < sfft; sample++)
            vHanWinCoeff[sample] = 0.5 * (1.0 - cos(2.0 * 3.14159265358979323846*(double)(sample) / ((double)sfft)));

        std::cout << "make_hanning_window\n";
        hanWinCoeff = arr(sfft, 1, &vHanWinCoeff[0]).as(c32);
    }

    ~stft_AF()
    {

    }

    // in: time domain windowed. dim: (soverlap x nch)
    // out: freq domain frame. dim: (sfreq x nch)
    arr stft(arr& in)
    {
        //std::cout << in.dims() << std::endl;
        //std::cout << window.dims() << std::endl;

        //std::cout << "sfft: " << sfft << std::endl;
        //std::cout << "soverlap: " << soverlap << std::endl;
        //

        //std::cout << "here2\n";
        // shifting
        window(seq(0, (sfft-soverlap-1)), af::span) = window(seq(soverlap, af::end), af::span);
        //std::cout << "here2\n";
        // copying new samples
        window(seq((sfft - soverlap), af::end), af::span) = in;
        //std::cout << "here2\n";        // batch windowing
        arr windowed = af::batchFunc(window, hanWinCoeff, af::operator*);
        //std::cout << "here2\n";
        // batch fft
        arr fftd = af::fft(windowed, sfft);
        //std::cout << "here2\n";
        // return cropped last (1 + half) -> size: (sfft / 2 + 1)
        //arr cropped = fftd(seq(0, sfreq), af::span);

        return fftd(seq(0, sfreq-1), af::span);
    }
};



int main(void)
{
    const int sfft = 512;
    const int sfreq = sfft / 2+1;
    const int soverlap = sfft / 2;

    const int nCh = 7;

    stft_AF* pstft = new stft_AF(sfft, nCh, soverlap);

    float** buffer;
    int nSamples = wm::wread("wav_7ch.wav", buffer);

    int nbuffers = nSamples / soverlap;

    //std::cout << "nbuffers: " << nbuffers << std::endl;

    arr frames = arr(sfreq, nCh, nbuffers, c32);

    arr samples = af::constant(0.0f, nSamples, nCh, c32);
    //std::cout << "herehere\n";
    for (int i = 0; i < nCh; i++)
    {
        samples(af::span, i) = arr(nSamples, (&buffer[i][0]), afHost).as(c32);
    }
    //std::cout << "herehere\n";

    af::timer::start();
    for (int i = 0; i < nbuffers; i++)
    {
        int start = i * soverlap;
        //std::cout << "here\n";
        arr arrived = samples(seq(start, start + soverlap - 1), af::span);
        //std::cout << "here\n";
        frames(af::span, af::span, i) = pstft->stft(arrived);
    }
    printf("elapsed seconds per 1 frame: %g\n", af::timer::stop() / (double)nbuffers);


    //std::cout << "herehere\n";
    af::write_complex64_binary(frames, "stfted");
    

    delete pstft;

    delete buffer;

    return 0;
}