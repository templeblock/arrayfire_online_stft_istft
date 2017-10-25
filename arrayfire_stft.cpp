#pragma once

#include "arrayfire_stft.hpp"

namespace af
{
    stft::stft(const int sfft, const int nch, const int sshift)
    {
        // make hanning window
        initialize(sfft, nch, sshift);
    }

    void stft::initialize(const int _sfft, const int _nch, const int _sshift)
    {
        sfft = _sfft;
        sshift = _sshift;
        nfreq = _sfft / 2 + 1;
        nch = _nch;

        window = af::constant(0.0f, sfft, nch, c32);
        outWindow = af::constant(0.0f, sfft, nch, f32);

        make_hanning_window();
    }

    void stft::make_hanning_window()
    {
        std::cout << "make_hanning_window\n";

        vHanHalfWinCoeff.resize(sfft / 2 + 1);

        for (int sample = 0; sample < sfft / 2 + 1; sample++)
            vHanHalfWinCoeff[sample] = 0.5 * (1.0 - std::cos(2.0 * 3.14159265358979323846*(double)(sample) / ((double)sfft)));

        hanWinCoeff = af::constant(0.0f, sfft, c32);
        hanHalfWinCoeff = arr(sfft / 2 + 1, &vHanHalfWinCoeff[0]).as(c32);

        hanWinCoeff(seq(0, sfft / 2)) = hanHalfWinCoeff;
        hanWinCoeff(seq(sfft / 2 + 1, sfft - 1, 1)) = hanHalfWinCoeff(seq(sfft / 2 - 1, 1, -1));

        hanWinCoeff = sshift*hanWinCoeff / af::tile(af::sum(hanWinCoeff), sfft, 1);
    }

    stft::~stft()
    {

    }

    // in: time domain windowed. dim: (sshift x nch)
    // out: freq domain frame. dim: (nfreq x nch)

    af::array stft::run(arr in)
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

    af::array stft::inverse(af::array& frame)
    {

        //shifting
        outWindow(seq(0, af::end - sshift), af::span) = outWindow(seq(sshift, af::end), af::span);

        outWindow(seq(af::end - sshift, af::end), af::span) = 0.0f;

        //adding
        outWindow = outWindow + af::real(ifft(af::join(0, frame, af::conjg(frame(seq(af::end - 1, 1, -1), af::span)))));

        return outWindow(seq(0, sshift - 1), af::span);
    }


	// in: freq domain frames. dim: (nfreq x nframes x nch)
	// out: time signal. dim: (nsamples x nch)
	af::array stft::batch_inverse(af::array batch)
	{
		const int nframes = batch.dims(1);

		//cout << batch.dims() << endl;
		af::array out = af::constant(0.0f, sshift * nframes, nch, f32);

		for (size_t i = 0; i < nframes; i++)
		{
			out(seq(i * sshift, (i + 1)*sshift - 1), af::span) = inverse(arr(batch(af::span, i, af::span),nfreq, nch));
		}

		return out;
	}
}

