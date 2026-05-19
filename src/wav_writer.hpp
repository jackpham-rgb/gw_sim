// wav_writer.hpp -- minimal mono 16-bit PCM WAV file writer.
// No external dependencies; the WAV format itself is just a 44-byte
// header followed by raw little-endian PCM samples.

#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace gw {

// Writes `samples` (each in [-1, 1]) to `filename` as a mono 16-bit PCM WAV
// at the given sample_rate. Returns true on success.
inline bool write_wav(const std::string& filename,
                      const std::vector<double>& samples,
                      int sample_rate)
{
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    const uint32_t num_samples     = static_cast<uint32_t>(samples.size());
    const uint16_t num_channels    = 1;
    const uint16_t bits_per_sample = 16;
    const uint32_t byte_rate       = sample_rate * num_channels * bits_per_sample / 8;
    const uint16_t block_align     = num_channels * bits_per_sample / 8;
    const uint32_t data_size       = num_samples * num_channels * bits_per_sample / 8;
    const uint32_t chunk_size      = 36 + data_size;

    auto write_u32 = [&](uint32_t v) {
        f.write(reinterpret_cast<const char*>(&v), 4);
    };
    auto write_u16 = [&](uint16_t v) {
        f.write(reinterpret_cast<const char*>(&v), 2);
    };

    // RIFF header
    f.write("RIFF", 4);
    write_u32(chunk_size);
    f.write("WAVE", 4);

    // fmt sub-chunk
    f.write("fmt ", 4);
    write_u32(16);                // PCM fmt chunk size
    write_u16(1);                 // PCM format code
    write_u16(num_channels);
    write_u32(static_cast<uint32_t>(sample_rate));
    write_u32(byte_rate);
    write_u16(block_align);
    write_u16(bits_per_sample);

    // data sub-chunk
    f.write("data", 4);
    write_u32(data_size);

    // Samples
    for (double s : samples) {
        if (s >  1.0) s =  1.0;
        if (s < -1.0) s = -1.0;
        int16_t v = static_cast<int16_t>(s * 32767.0);
        f.write(reinterpret_cast<const char*>(&v), 2);
    }

    return f.good();
}

} // namespace gw
