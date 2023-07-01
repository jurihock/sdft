@echo off

set DFTSIZE=1000
set HOPSIZE=100
set WINDOW=hann
set LATENCY=1

set BUILD=build

mkdir %BUILD%

pushd %BUILD%
cmake -A x64 .. || exit /b
cmake --build . --config Release || exit /b
popd

del /f "test.*.dft"
del /f "test.*.wav"

%BUILD%\Release\sdft-test-c.exe   %DFTSIZE% %HOPSIZE% %WINDOW% %LATENCY% "test.wav" "test.c.wav"   "test.c.dft"
%BUILD%\Release\sdft-test-cpp.exe %DFTSIZE% %HOPSIZE% %WINDOW% %LATENCY% "test.wav" "test.cpp.wav" "test.cpp.dft"
python test.py                    %DFTSIZE% %HOPSIZE% %WINDOW% %LATENCY% "test.wav" "test.py.wav"  "test.py.dft"
python main.py                    %DFTSIZE% %HOPSIZE% %WINDOW% %LATENCY% "test.wav" "test.{}.wav"  "test.{}.dft"
