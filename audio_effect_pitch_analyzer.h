/*************************************************************************/
/*  audio_effect_spectrum_analyzer.h                                     */
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

#pragma once

#include "servers/audio/audio_effect.h"
#include "core/templates/vector.h"
#include "core/math/vector2.h"

#include <utility>
#include <vector>
#include <complex>

class AudioEffectPitchAnalyzer;

class AudioEffectPitchAnalyzerInstance : public AudioEffectInstance {
	GDCLASS(AudioEffectPitchAnalyzerInstance, AudioEffectInstance);

public:

	// pair elements will contain frequency, value pairs. .x : frequency, .y : volume.
	using peak_t = Vector2;

private:

	friend class AudioEffectPitchAnalyzer;
	Ref<AudioEffectPitchAnalyzer> base;

	std::vector<std::vector<peak_t>> peak_history;

	// temporary buffers to process data.
	std::vector<std::complex<float>> temporal_fft1;
	std::vector<std::complex<float>> temporal_fft2;
	std::vector<peak_t> peak_scratch;
	std::vector<float> scratchb1;
	std::vector<float> scratchb2;
	int temporal_fft_pos;

	// properties of the audio stream and processing
	int fft_size;
	int fft_count;
	int fft_pos;
	float mix_rate;
	float clarity_threshold;
	uint64_t last_fft_time;

	static void _find_peaks(const std::vector<std::complex<float>> &points,std::vector<peak_t> &out);
	static void _limit_peaks(std::vector<peak_t> &peaks, float cl_trh, std::vector<peak_t> &scratch);

protected:
	static void _bind_methods();
public:
	virtual void process(const AudioFrame *p_src_frames, AudioFrame *p_dst_frames, int p_frame_count) override;
	PackedVector2Array get_peaks() const;
};

class AudioEffectPitchAnalyzer : public AudioEffect {
	GDCLASS(AudioEffectPitchAnalyzer, AudioEffect);

public:
	enum FFTSize {
		FFT_SIZE_256,
		FFT_SIZE_512,
		FFT_SIZE_1024,
		FFT_SIZE_2048,
		FFT_SIZE_4096,
		FFT_SIZE_MAX
	};

public:
	friend class AudioEffectPitchAnalyzerInstance;
	float buffer_length;
	float tapback_pos;
	float clarity_threshold;
	FFTSize fft_size;
protected:
	static void _bind_methods();

public:
	Ref<AudioEffectInstance> instantiate() override;
	void set_buffer_length(float p_seconds);
	float get_buffer_length() const;
	void set_tap_back_pos(float p_seconds);
	float get_tap_back_pos() const;
	void set_clarity_threshold(float value);
	float get_clarity_threshold() const;

	void set_fft_size(FFTSize);
	FFTSize get_fft_size() const;

	AudioEffectPitchAnalyzer();

};

VARIANT_ENUM_CAST(AudioEffectPitchAnalyzer::FFTSize);
