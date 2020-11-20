# projectCS
This is a repository for our 15 cc project in CS

# Prerequisites
The project requires cmake >= `3.10` and a compiler supporting c++17.
For graphs to be shown, you also need `gnuplot` installed and in your path.

# Building
Unix:
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

Windows:
```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

This will output the executable as `dist/projectcs` (unix), or `dist/Release/projectcs.exe` (windows).

Compiler support for c++17 (specifically `std::filesystem`) on mac seems to vary, so YMMV building this on mac.