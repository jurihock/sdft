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

bool readwav(const char* path, float** data, size_t* size, double* samplerate)
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
