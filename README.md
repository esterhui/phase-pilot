# phase-pilot

Fast signal detection via FFT and tracking of tones with FLL/PLL and PRN codes with DLL

```sh
conan install . --build=missing       # install project dependencies
cmake --preset conan-release          # prepare build system
cmake --build --preset conan-release  # build
```

This will build & install debug versions of all the dependencies and set up a preset that targets the Debug cmake build.

```sh
conan install . --build=missing -s build_type=Debug
cmake --preset conan-debug
cmake --build --preset conan-debug
```
