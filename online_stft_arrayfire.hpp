#pragma once
#include <arrayfire.h>

class stft_AF
{
    typedef af::array arr;
    typedef af::dim4 dim4;
    typedef af::seq seq;

    int sfft;
    int nch;
    int sshift;
    int nfreq;

    arr window;
    arr outWindow;

    std::vector<double> vHanHalfWinCoeff;
    arr hanHalfWinCoeff;

    arr hanWinCoeff;

public:

    stft_AF(const int sfft, const int nch, const int sshift)
    {
        // make hanning window
        initialize(sfft, nch, sshift);
    }

    void initialize(const int _sfft, const int _nch, const int _sshift)
    {
        sfft = _sfft;
        sshift = _sshift;
        nfreq = _sfft / 2 + 1;
        nch = _nch;

        window = af::constant(0.0f, sfft, nch, c32);
        outWindow = af::constant(0.0f, sfft, nch, f32);

        make_hanning_window();
    }

    void make_hanning_window()
    {
        std::cout << "make_hanning_window\n";

        vHanHalfWinCoeff.resize(sfft / 2 + 1);

        for (int sample = 0; sample < sfft / 2 + 1; sample++)
            vHanHalfWinCoeff[sample] = 0.5 * (1.0 - cos(2.0 * 3.14159265358979323846*(double)(sample) / ((double)sfft)));

        hanWinCoeff = af::constant(0.0f, sfft, c32);
        hanHalfWinCoeff = arr(sfft / 2 + 1, &vHanHalfWinCoeff[0]).as(c32);

        hanWinCoeff(seq(0, sfft / 2)) = hanHalfWinCoeff;
        hanWinCoeff(seq(sfft / 2 + 1, sfft - 1, 1)) = hanHalfWinCoeff(seq(sfft / 2 - 1, 1, -1));

        hanWinCoeff = sshift*hanWinCoeff / af::tile(af::sum(hanWinCoeff), sfft, 1);
    }

    ~stft_AF()
    {

    }

    // in: time domain windowed. dim: (sshift x nch)
    // out: freq domain frame. dim: (nfreq x nch)
    arr stft(arr& in)
    {
        // shifting
        window(seq(0, (sfft - sshift - 1)), af::span) = window(seq(sshift, af::end), af::span);

        // copying new samples
        window(seq((sfft - sshift), af::end), af::span) = in;
        arr windowed = af::batchFunc(window, hanWinCoeff, af::operator*);

        // batch fft
        arr fftd = af::fft(windowed, sfft);

        // return cropped last (1 + half) -> size: (sfft / 2 + 1)
        return fftd(seq(0, nfreq - 1), af::span);
    }

    // in: freq domain frame. dim: (nfreq x nch)
    // out: time domain windowed. dim: (sshift x nch)
    af::array istft(af::array& frame, const bool isEnd)
    {

        //shifting
        outWindow(seq(0, af::end - sshift), af::span) = outWindow(seq(sshift, af::end), af::span);

        outWindow(seq(af::end - sshift, af::end), af::span) = 0.0f;

        //adding
        outWindow = outWindow + af::real(ifft(af::join(0, frame, af::conjg(frame(seq(af::end - 1, 1, -1), af::span)))));

        return outWindow(seq(0, sshift - 1), af::span);
    }

};