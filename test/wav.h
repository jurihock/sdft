#pragma once

#include <stdbool.h>
#include <stdlib.h>

#include <dr_wav.h>

float* allocwav(const size_t size)
{
  #if defined(__cplusplus)
    return new float[size];
  #else
    return malloc(size * sizeof(float));
  #endif
}

void freewav(float* data)
{
  #if defined(__cplusplus)
    delete[] data;
  #else
    free(data);
  #endif
}

drwav_uint64 drwav_write_pcm_frames_f32_to_s32(drwav* wav, const drwav_uint64 samples, const float* data)
{
  drwav_uint64 size = samples * wav->channels;
  drwav_int32* data_int32 = (drwav_int32*)malloc(size * sizeof(drwav_int32));

  drwav_f32_to_s32(data_int32, data, size);
  size = drwav_write_pcm_frames(wav, samples, data_int32);

  free(data_int32);
  return size;
}

bool readwav(const char* path, float** data, size_t* size, size_t* samplerate)
{
  (*data) = NULL;
  (*size) = 0;
  (*samplerate) = 0;

  drwav wav;

  if (drwav_init_file(&wav, path, NULL) != DRWAV_TRUE)
  {
    return false;
  }

  const size_t samples = wav.totalPCMFrameCount;
  const size_t channels = wav.channels;
  const size_t bytes = samples * channels * sizeof(float);

  if (bytes > DRWAV_SIZE_MAX)
  {
    drwav_uninit(&wav);

    return false;
  }

  (*data) = allocwav(samples * channels);
  (*size) = samples * channels;

  if (drwav_read_pcm_frames_f32(&wav, samples, (*data)) != samples)
  {
    drwav_uninit(&wav);

    freewav(*data);

    (*data) = NULL;
    (*size) = 0;

    return false;
  }

  if (channels > 1)
  {
    for (size_t i = 0; i < samples; ++i)
    {
      (*data)[i] = (*data)[i * channels];

      for (size_t j = 1; j < channels; ++j)
      {
        (*data)[i] += (*data)[i * channels + j];
      }

      (*data)[i] /= channels;
    }

    (*size) = samples;
  }

  (*samplerate) = wav.sampleRate;

  drwav_uninit(&wav);

  return true;
}

bool writewav(const char* path, const float* data, const size_t size, const size_t samplerate)
{
  const size_t samples = size;
  const size_t channels = 1;
  const size_t bytes = samples * channels * sizeof(drwav_int32);

  if (bytes > DRWAV_SIZE_MAX)
  {
    return false;
  }

  drwav wav;
  drwav_data_format format;

  format.container = drwav_container_riff;
  format.format = DR_WAVE_FORMAT_PCM;
  format.bitsPerSample = sizeof(drwav_uint32) * 8;
  format.channels = (drwav_uint16)channels;
  format.sampleRate = (drwav_uint32)samplerate;

  if (drwav_init_file_write(&wav, path, &format, NULL) != DRWAV_TRUE)
  {
    return false;
  }

  if (drwav_write_pcm_frames_f32_to_s32(&wav, samples, data) != samples)
  {
    drwav_uninit(&wav);

    return false;
  }

  drwav_uninit(&wav);

  return true;
}
