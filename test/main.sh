#!/bin/bash

DFTSIZE=1000
HOPSIZE=100
WINDOW=hann
LATENCY=1

BUILD=build

cmake -DCMAKE_BUILD_TYPE=Release -B ${BUILD} || exit $?
cmake --build ${BUILD} || exit $?

rm -f "test.*.dft"
rm -f "test.*.wav"

${BUILD}/sdft-test-c   ${DFTSIZE} ${HOPSIZE} ${WINDOW} ${LATENCY} "test.wav" "test.c.wav"   "test.c.dft"
${BUILD}/sdft-test-cpp ${DFTSIZE} ${HOPSIZE} ${WINDOW} ${LATENCY} "test.wav" "test.cpp.wav" "test.cpp.dft"
uv run test.py         ${DFTSIZE} ${HOPSIZE} ${WINDOW} ${LATENCY} "test.wav" "test.py.wav"  "test.py.dft"
uv run main.py         ${DFTSIZE} ${HOPSIZE} ${WINDOW} ${LATENCY} "test.wav" "test.{}.wav"  "test.{}.dft"
