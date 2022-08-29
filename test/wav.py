import numpy
import sys
import wave


def readwav(path):

    with wave.open(path, 'rb') as file:
        sr = file.getframerate()
        bytes = file.getsampwidth()
        channels = file.getnchannels()
        data = file.readframes(file.getnframes())

    assert bytes in [1, 2, 3, 4]
    bits = bytes * 8
    scaler = 2 ** (bits - 1) - 1

    data = numpy.frombuffer(data, dtype=numpy.uint8).reshape(-1, bytes)
    data = numpy.asarray([
        int.from_bytes(frame, signed=(bits != 8), byteorder=sys.byteorder)
        for frame in data])
    data = data.astype(float).reshape(-1, channels)

    data -= 128 if bits == 8 else 0
    data = (data + 0.5) / (scaler + 0.5)
    data = data.clip(-1, +1)

    data = data.mean(axis=-1)

    return data, sr


def writewav(path, data, sr, bits=32):

    data = numpy.asarray(data)

    assert data.dtype == float
    assert data.ndim == 1
    assert data.size > 0

    assert bits in [8, 16, 24, 32]
    bytes = bits // 8
    scaler = 2 ** (bits - 1) - 1

    data = data.clip(-1, +1)
    data = (data * (scaler + 0.5)) - 0.5
    data += 128 if bits == 8 else 0

    data = b''.join([
        int(frame).to_bytes(length=bytes, signed=(bits != 8), byteorder=sys.byteorder)
        for frame in data])

    with wave.open(path, 'wb') as file:
        file.setframerate(sr)
        file.setsampwidth(bytes)
        file.setnchannels(1)
        file.writeframes(data)
