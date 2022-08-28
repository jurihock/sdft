#!/bin/bash

DFTSIZE=512
HOPSIZE=1000

BUILD=../build-release

mkdir -p ${BUILD}

pushd ${BUILD} >/dev/null 2>&1
cmake -DCMAKE_BUILD_TYPE=Release .. || exit $?
cmake --build . || exit $?
popd >/dev/null 2>&1

rm -f "test.*.dft"
rm -f "test.*.wav"

${BUILD}/sdft-test-c   ${DFTSIZE} ${HOPSIZE} "test.wav" "test.c.wav" "test.c.dft"
${BUILD}/sdft-test-cpp ${DFTSIZE} ${HOPSIZE} "test.wav" "test.cpp.wav" "test.cpp.dft"
python3 test.py        ${DFTSIZE} ${HOPSIZE} "test.wav" "test.py.wav" "test.py.dft"
python3 main.py        ${DFTSIZE} ${HOPSIZE} "test.wav" "test.{}.wav" "test.{}.dft"
