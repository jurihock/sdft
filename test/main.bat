@echo off

set DFTSIZE=1000
set HOPSIZE=100
set WINDOW=hann
set LATENCY=1

set BUILD=build

cmake -A x64 -B %BUILD% || exit /b
cmake --build %BUILD% --config Release || exit /b

del /f "test.*.dft"
del /f "test.*.wav"

%BUILD%\Release\sdft-test-c.exe   %DFTSIZE% %HOPSIZE% %WINDOW% %LATENCY% "test.wav" "test.c.wav"   "test.c.dft"
%BUILD%\Release\sdft-test-cpp.exe %DFTSIZE% %HOPSIZE% %WINDOW% %LATENCY% "test.wav" "test.cpp.wav" "test.cpp.dft"
uv run test.py                    %DFTSIZE% %HOPSIZE% %WINDOW% %LATENCY% "test.wav" "test.py.wav"  "test.py.dft"
uv run main.py                    %DFTSIZE% %HOPSIZE% %WINDOW% %LATENCY% "test.wav" "test.{}.wav"  "test.{}.dft"
