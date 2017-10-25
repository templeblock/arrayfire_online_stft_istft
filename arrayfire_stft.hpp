#pragma once
#include <arrayfire.h>

namespace af
{
    class stft
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

        stft(const int sfft, const int nch, const int sshift);

        void initialize(const int _sfft, const int _nch, const int _sshift);

        void make_hanning_window();

        ~stft();

        // in: time domain windowed. dim: (sshift x nch)
        // out: freq domain frame. dim: (nfreq x nch)
        arr run(arr in);

        // in: freq domain frame. dim: (nfreq x nch)
        // out: time domain windowed. dim: (sshift x nch)
        arr inverse(arr& frame);

		arr batch_inverse(arr batch);
    };
}
