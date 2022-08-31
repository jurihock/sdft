@echo off

set DFTSIZE=1000
set HOPSIZE=100

set BUILD=..\build-release

mkdir %BUILD%

pushd %BUILD%
cmake -A x64 .. || exit /b
cmake --build . --config Release || exit /b
popd

del /f "test.*.dft"
del /f "test.*.wav"

%BUILD%\Release\sdft-test-c.exe   %DFTSIZE% %HOPSIZE% "test.wav" "test.c.wav"   "test.c.dft"
%BUILD%\Release\sdft-test-cpp.exe %DFTSIZE% %HOPSIZE% "test.wav" "test.cpp.wav" "test.cpp.dft"
python test.py                    %DFTSIZE% %HOPSIZE% "test.wav" "test.py.wav"  "test.py.dft"
python main.py                    %DFTSIZE% %HOPSIZE% "test.wav" "test.{}.wav"  "test.{}.dft"
