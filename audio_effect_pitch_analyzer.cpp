/*************************************************************************/
/*  audio_effect_spectrum_analyzer.cpp                                   */
/*************************************************************************/
/*                       This file is part of:                           */
/*                           GODOT ENGINE                                */
/*                      https://godotengine.org                          */
/*************************************************************************/
/* Copyright (c) 2007-2022 Juan Linietsky, Ariel Manzur.                 */
/* Copyright (c) 2014-2022 Godot Engine contributors (cf. AUTHORS.md).   */
/*                                                                       */
/* Permission is hereby granted, free of charge, to any person obtaining */
/* a copy of this software and associated documentation files (the       */
/* "Software"), to deal in the Software without restriction, including   */
/* without limitation the rights to use, copy, modify, merge, publish,   */
/* distribute, sublicense, and/or sell copies of the Software, and to    */
/* permit persons to whom the Software is furnished to do so, subject to */
/* the following conditions:                                             */
/*                                                                       */
/* The above copyright notice and this permission notice shall be        */
/* included in all copies or substantial portions of the Software.       */
/*                                                                       */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.*/
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE     */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                */
/*************************************************************************/

#include "audio_effect_pitch_analyzer.h"
#include "servers/audio_server.h"
#include <algorithm>

template <typename T> using fft_arg = std::vector<std::complex<T>>;


int findMSB(int x)
{
    int p = 0;

    while (x > 1) {
        x>>= 1;
        ++p;
    }

    return p;
}

enum PeakCorrection {
    Quadratic,
    None,
};

/*
 *  Bit reverse an integer given a word of nb bits
 *  NOTE: Only works for 32-bit words max
 *  examples:
 *  10b      -> 01b
 *  101b     -> 101b
 *  1011b    -> 1101b
 *  0111001b -> 1001110b
 */
int bitr(uint32_t x, int nb)
{
    x = ( x               << 16) | ( x               >> 16);
    x = ((x & 0x00FF00FF) <<  8) | ((x & 0xFF00FF00) >>  8);
    x = ((x & 0x0F0F0F0F) <<  4) | ((x & 0xF0F0F0F0) >>  4);
    x = ((x & 0x33333333) <<  2) | ((x & 0xCCCCCCCC) >>  2);
    x = ((x & 0x55555555) <<  1) | ((x & 0xAAAAAAAA) >>  1);

    return ((x >> (32 - nb)) & (0xFFFFFFFF >> (32 - nb)));
}

/*
 * Computes a Fourier transform, i.e.,
 * xo[k] = 1/sqrt(N) sum(j=0 -> N-1) xi[j] exp(i 2pi j k / N)
 * with O(N log N) complexity using the butterfly technique
 *
 * NOTE: Only works for arrays whose size is a power-of-two.
 * 
 * License : * modified from "https://github.com/jdupuy/dj_fft" which was licenced under MIT.
 */
template <typename T> fft_arg<T> fft1d(const fft_arg<T> &xi, fft_arg<T> &xo, int dir)
{
	int cnt = xi.size();
    int msb = findMSB(cnt);
    T nrm = 1.0 / sqrt(T(cnt));
	xo.resize(cnt);

    // pre-process the input data
    for (int j = 0; j < cnt; ++j)
        xo[j] = nrm * xi[bitr(j, msb)];//dir == -1 ? nrm * xi[bitr(j, msb)] : xi[bitr(j, msb)];

    // fft passes
    for (int i = 0; i < msb; ++i) {
        int bm = 1 << i; // butterfly mask
        int bw = 2 << i; // butterfly width
        T ang = T(dir) * M_PI / T(bm); // precomputation

        // fft butterflies
        for (int j = 0; j < (cnt/2); ++j) {
            int i1 = ((j >> i) << (i + 1)) + j % bm; // left wing
            int i2 = i1 ^ bm;                        // right wing
            std::complex<T> z1 = std::polar(T(1), ang * T(i1 ^ bw)); // left wing rotation
            std::complex<T> z2 = std::polar(T(1), ang * T(i2 ^ bw)); // right wing rotation
            std::complex<T> tmp = xo[i1];

            xo[i1]+= z1 * xo[i2];
            xo[i2] = tmp + z2 * xo[i2];
        }
    }

    return xo;
}





void _autocorrelation(fft_arg<float> &frames, fft_arg<float> &scratch) {
	fft1d(frames, scratch, -1);
	for(auto it = scratch.begin(); it != scratch.end(); ++it)
	{
		float rl = it->real();
		float im = it->imag();
		*it = rl*rl + im*im;
	}
	fft1d(scratch, frames, 1);
}

// TODO : If needed the clarity threshold can be replaced by a Callable 
//   argument to give different thresholds for different frequencies.
//
// directly takes the array of complex data to avoid copying.
void AudioEffectPitchAnalyzerInstance::_find_peaks( const std::vector<std::complex<float>> &points,std::vector<peak_t> &out) {

	// find peaks.

	size_t idx = 0;
	auto count = points.size();
	out.resize(0);
	out.reserve(count);

	// skip over the initial positive values :
	while(idx != count && !(points[idx].real() < 0.f)){
		++idx;
	}

	if(idx == count) {
		return;
	}
	
	while(true) {
		// detect peaks	

		// skip over the negative parts afterwards :
		while(idx != count && !(points[idx].real() > 0.f)){
			++idx;
		}

		float max = 0;
		size_t midx = -1;

		// find the peak before going into negative territory again.
		while(idx != count && points[idx].real() > 0.f) {
			float v = points[idx].real();
			if (v > max) {
				max = v;
				midx = idx;
			}
			++idx;
		}

		if(idx == count) {
			return;
		}

		// add peak if we didn't reach the end.
		out.push_back(peak_t((float)midx, max));
	}
}

void AudioEffectPitchAnalyzerInstance::_limit_peaks(std::vector<peak_t> &peaks, float cl_thr, std::vector<peak_t> &scratch) {
	
	scratch.clear();
	scratch.swap(peaks);

	for(auto it = scratch.begin(); it != scratch.end(); ++it) {
		if(it->y >= cl_thr) {
			peaks.push_back(*it);
		}
	}
}

void AudioEffectPitchAnalyzerInstance::process(const AudioFrame *p_src_frames, AudioFrame *p_dst_frames, int p_frame_count) {
	uint64_t time = OS::get_singleton()->get_ticks_usec();


	//copy everything over first, since this only really does capture
	for (int i = 0; i < p_frame_count; i++) {
		p_dst_frames[i] = p_src_frames[i];
	}

	//capture spectrum
	while (p_frame_count) {
		int to_fill = fft_size - temporal_fft_pos;
		to_fill = MIN(to_fill, p_frame_count);
		const double to_fill_step = Math_TAU / (double)fft_size;


		for (int i = 0; i < to_fill; i++) { //left buffers
			// MATH : Window functions should not be used when doing autocorellation. (Hann smoothing was used in the original code)
			temporal_fft1[temporal_fft_pos] = p_src_frames->l;
			
			++p_src_frames;
			++temporal_fft_pos;
		}

		p_frame_count -= to_fill;

		if (temporal_fft_pos == fft_size) {
			int next = (fft_pos + 1) % fft_count;

			//time to do a FFT
			// FIXME : this is just experimental. i'm sure there is a more efficient way to do this.
			//float *dst_buff = (float *)peak_history[next].ptr(); //do not use write, avoid cow
			//AudioFrame *hw = (AudioFrame *)fft_history[next].ptr(); //do not use write, avoid cow
			_autocorrelation(temporal_fft1, temporal_fft2);
			std::vector<peak_t> &dst_peaks = peak_history[next];
			float cl = clarity_threshold * temporal_fft1[0].real();
			_find_peaks(temporal_fft1, dst_peaks);
			_limit_peaks(dst_peaks, cl, peak_scratch);

			// finalizing peak values :
			for(auto it = dst_peaks.begin(); it != dst_peaks.end(); ++it) {
				it->x = mix_rate/it->x;
			}

			fft_pos = next; //swap
			temporal_fft_pos = 0;
		}
	}

	//determine time of capture
	double remainer_sec = (temporal_fft_pos / mix_rate); //subtract remainder from mix time
	last_fft_time = time - uint64_t(remainer_sec * 1000000.0);
}

void AudioEffectPitchAnalyzerInstance::_bind_methods() {
	ClassDB::bind_method(D_METHOD("get_peaks"), &AudioEffectPitchAnalyzerInstance::get_peaks);
}

// Vector2 AudioEffectPitchAnalyzerInstance::get_magnitude_for_frequency_range(float p_begin, float p_end, MagnitudeMode p_mode) const {
PackedVector2Array AudioEffectPitchAnalyzerInstance::get_peaks() const {
	if (last_fft_time == 0) {
		return PackedVector2Array();
	}
	uint64_t time = OS::get_singleton()->get_ticks_usec();
	float diff = double(time - last_fft_time) / 1000000.0 + base->get_tap_back_pos();
	diff -= AudioServer::get_singleton()->get_output_latency();
	float fft_time_size = float(fft_size) / mix_rate;

	int fft_index = fft_pos;

	while (diff > fft_time_size) {
		diff -= fft_time_size;
		fft_index -= 1;
		if (fft_index < 0) {
			fft_index = fft_count - 1;
		}
	}

	const std::vector<peak_t> &peaks = peak_history[fft_index];
	PackedVector2Array ret;
	size_t cnt = peaks.size();
	if(cnt == 0) {
		return ret;
	}
	ret.resize(cnt);
	size_t s = ret.size();

	Vector2 *retw = ret.ptrw();
	memcpy(retw, peaks.data(), sizeof(Vector2)*cnt);
	return ret;
}

Ref<AudioEffectInstance> AudioEffectPitchAnalyzer::instantiate() {
	Ref<AudioEffectPitchAnalyzerInstance> ins;
	ins.instantiate();

	{
		ins->base = Ref<AudioEffectPitchAnalyzer>(this);
		static const int fft_sizes[FFT_SIZE_MAX] = { 256, 512, 1024, 2048, 4096 };
		ins->fft_size = fft_sizes[fft_size];
		ins->mix_rate = AudioServer::get_singleton()->get_mix_rate();
		ins->fft_count = (buffer_length / (float(ins->fft_size) / ins->mix_rate)) + 1;
		ins->fft_pos = 0;
		ins->last_fft_time = 0;
		ins->peak_history.resize(ins->fft_count);
		ins->temporal_fft1.resize(ins->fft_size);
		ins->temporal_fft2.resize(ins->fft_size);
		ins->scratchb1.resize(ins->fft_size);
		ins->scratchb2.resize(ins->fft_size);
		ins->peak_scratch.resize(ins->fft_size);
		ins->temporal_fft_pos = 0;
		ins->clarity_threshold = clarity_threshold;
		for (int i = 0; i < ins->fft_count; i++) {
			ins->peak_history[i].resize(ins->fft_size, Vector2(0,0));
		}
	}

	return ins;
}

void AudioEffectPitchAnalyzer::set_buffer_length(float p_seconds) {
	buffer_length = p_seconds;
}

float AudioEffectPitchAnalyzer::get_buffer_length() const {
	return buffer_length;
}

void AudioEffectPitchAnalyzer::set_tap_back_pos(float p_seconds) {
	tapback_pos = p_seconds;
}

float AudioEffectPitchAnalyzer::get_tap_back_pos() const {
	return tapback_pos;
}

void AudioEffectPitchAnalyzer::set_fft_size(FFTSize p_fft_size) {
	ERR_FAIL_INDEX(p_fft_size, FFT_SIZE_MAX);
	fft_size = p_fft_size;
}

AudioEffectPitchAnalyzer::FFTSize AudioEffectPitchAnalyzer::get_fft_size() const {
	auto ret = fft_size;
	
	return ret;
}

void AudioEffectPitchAnalyzer::set_clarity_threshold(float value) {
	clarity_threshold = value;
}

float AudioEffectPitchAnalyzer::get_clarity_threshold() const {
	return clarity_threshold;
}

void AudioEffectPitchAnalyzer::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_buffer_length", "seconds"), &AudioEffectPitchAnalyzer::set_buffer_length);
	ClassDB::bind_method(D_METHOD("get_buffer_length"), &AudioEffectPitchAnalyzer::get_buffer_length);

	ClassDB::bind_method(D_METHOD("set_tap_back_pos", "seconds"), &AudioEffectPitchAnalyzer::set_tap_back_pos);
	ClassDB::bind_method(D_METHOD("get_tap_back_pos"), &AudioEffectPitchAnalyzer::get_tap_back_pos);

	ClassDB::bind_method(D_METHOD("set_fft_size", "size"), &AudioEffectPitchAnalyzer::set_fft_size);
	ClassDB::bind_method(D_METHOD("get_fft_size"), &AudioEffectPitchAnalyzer::get_fft_size);

	ClassDB::bind_method(D_METHOD("set_clarity_threshold", "value"), &AudioEffectPitchAnalyzer::set_clarity_threshold);
	ClassDB::bind_method(D_METHOD("get_clarity_threshold"), &AudioEffectPitchAnalyzer::get_clarity_threshold);

	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "buffer_length", PROPERTY_HINT_RANGE, "0.1,4,0.1,suffix:s"), "set_buffer_length", "get_buffer_length");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "tap_back_pos", PROPERTY_HINT_RANGE, "0.1,4,0.1"), "set_tap_back_pos", "get_tap_back_pos");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "fft_size", PROPERTY_HINT_ENUM, "256,512,1024,2048,4096"), "set_fft_size", "get_fft_size");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "clarity_threshold", PROPERTY_HINT_RANGE, "0.001,200.0"), "set_clarity_threshold", "get_clarity_threshold");

	BIND_ENUM_CONSTANT(FFT_SIZE_256);
	BIND_ENUM_CONSTANT(FFT_SIZE_512);
	BIND_ENUM_CONSTANT(FFT_SIZE_1024);
	BIND_ENUM_CONSTANT(FFT_SIZE_2048);
	BIND_ENUM_CONSTANT(FFT_SIZE_4096);
	BIND_ENUM_CONSTANT(FFT_SIZE_MAX);
}

AudioEffectPitchAnalyzer::AudioEffectPitchAnalyzer() {
	buffer_length = 2;
	tapback_pos = 0.01;
	fft_size = FFT_SIZE_1024;
}
