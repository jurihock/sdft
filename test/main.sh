#!/bin/bash

DFTSIZE=512
HOPSIZE=1000

BUILD=../build-release

mkdir -p ${BUILD}

pushd ${BUILD} >/dev/null 2>&1
cmake -DCMAKE_BUILD_TYPE=Release .. || exit $?
cmake --build . || exit $?
popd >/dev/null 2>&1

rm -f "test.*.dfts"

${BUILD}/sdft-test-c   ${DFTSIZE} ${HOPSIZE} "test.wav" "test.c.dfts"
${BUILD}/sdft-test-cpp ${DFTSIZE} ${HOPSIZE} "test.wav" "test.cpp.dfts"
python3 test.py        ${DFTSIZE} ${HOPSIZE} "test.wav" "test.py.dfts"
python3 main.py        ${DFTSIZE} ${HOPSIZE} "test.wav" "test.{}.dfts"
