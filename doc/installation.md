# Install

This page summarizes the build and setup information from the repository
README.

## System Requirements

- Linux, with Ubuntu 16.04 or newer recommended
- macOS 10.10 or newer
- CMake 3.5 or newer
- GNU `gcc` and `gfortran` 14.2 or newer
- GNU make or make, used to build SIMPLE's pinned OpenBLAS
- FFTW 3.3 or newer, with double, single, and threaded builds
- libcurl 7 or newer
- libTIFF 4 or newer
- jbigkit 2 or newer, optional
- Python 3.7 or newer for the core build; NICE currently requires Python 3.10-3.13

SIMPLE always uses its pinned OpenBLAS release, which provides both BLAS and LAPACK. Its build scratch lives
in the build folder; the verified source archive and the installed static library persist outside it, the
archive under the gitignored `.deps/downloads/` (checked against its pinned SHA-256), the library under
`.deps/numerics/`. System BLAS and LAPACK libraries, including macOS Accelerate, are not selected.

With GNU C and Fortran compilers and OpenMP enabled, SIMPLE builds the OpenMP flavour of OpenBLAS: a product
called from serial code runs on SIMPLE's OpenMP threads, a product called inside one of SIMPLE's OpenMP teams
stays on the calling thread. CMake picks the GNU C compiler that sits beside the chosen gfortran
(`cmake/AutoDetectGCC.cmake`), so Homebrew GCC on macOS takes this path too. Only a C compiler set by hand that
cannot share the Fortran OpenMP runtime (AppleClang) gets the sequential flavour with its locked allocator:
every product runs on the calling thread. The dependency summary printed by CMake reports the version and
flavour. The cache key includes the platform,
architecture, compiler versions and flavour; the recipe record beside the library additionally includes the
compiler paths and build flags. Consequently `compile_clean.sh`, `compile_debug.sh` and other clean rebuilds
reuse a compatible library without compiling it again, and a changed recipe rebuilds it. Remove
`.deps/numerics/` only to force every cached configuration to rebuild.

The library is built with `NUM_THREADS=128`, independent of the machine that builds it. That caps a
serial-region product at 128 threads and sizes the library's buffer table, which every thread inside a BLAS
call holds a slot of: SIMPLE's teams stay within the core count, so any node with up to 128 cores is covered.
OpenBLAS does not wait at that limit: it warns on stderr ("precompiled NUM_THREADS exceeded"), then corrupts
its table and crashes. A larger node needs a larger pin in `cmake/Dependencies.cmake`.

## Obtain The Source

Download and unpack the stable release:

```shell
gunzip SIMPLE-3.0.0.tar.gz
tar -xvf SIMPLE-3.0.0.tar
cd SIMPLE-3.0.0
```

Or clone the repository to install the version under active development:

```shell
git clone https://github.com/hael/SIMPLE.git
cd SIMPLE
```

## Standard Build

```shell
mkdir build
cd build
cmake ..
make -j install
```

## Environment Setup

### Bash Shell

BASH users should append the contents of `add2.bashrc` to the bash configuration file `~/.bashrc` with the following command:

```shell
cat add2.bashrc >> ~/.bashrc
```

### C shell

TCSH/CSH users can use the following command instead:

```shell
cat add2.tcshrc >> ~/.tcshrc
```

### Manual Configuration

Alternatively, the following changes can also be made manually in the default configuration file for the shell such as `~/.bashrc`,  `~/.tcshrc` or `~/.zshrc` to setup the relevant environment variables.

Add these paths to the `PATH` environemnt variable:

```text
<simple_path>/build/bin
<simple_path>/build/scripts
```

Set the `SIMPLE_PATH` environment vairable:

```text
SIMPLE_PATH=<simple_path>
```

## Custom Installation

To change the installation directory, use the following Cmake option before running make install:

```shell
cmake -DCMAKE_INSTALL_PREFIX=<install_dir> ..
make -j install
```

In case of non-standard location for gcc/gfortran compilers or FFTW, provide the paths to Cmake as follows:

```shell
FC=<gfortran_path> CC=<gcc_path> FFTW_DIR=<fftw_path> cmake ..
make -j install
```

For instance, on MacOS such locations could be:

- Macports:   FC=/opt/local/bin/gfortran FFTW_DIR=/opt/local
- Fink users: FC=/sw/bin/gfortran FFTW_DIR=/sw/
- Homebrew:   FC=/usr/local/bin/gfortran FFTW_DIR=/usr/local/

For more advanced FFTW installations, see *Installation of FFTW* below.

## Developer Installation and Updates

Get the developemnt version of SIMPLE by cloning the git repository:

```shell
git clone https://github.com/hael/SIMPLE.git
cd SIMPLE
```

Run the following commands to install the debug build type:

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=debug ..
make -j install
```

Or use one of the helper build script such as `compile_clean.sh` provided in the SIMPLE directory for various installation cases. Install SIMPLE with the `compile_debug.sh` shell script

```bash
bash compile_debug.sh
```

## GUI Build With NICE (New Interface for CryoEM)

Use the `compile_gui.sh` script to install SIMPLE with the NICE GUI support:

```bash
bash compile_gui.sh
```

NICE currently supports Python 3.10 through 3.13 because its pinned NumPy
2.2.6 provides wheels for those versions. The script prefers Python 3.12 and
passes the selected interpreter explicitly to CMake. Set
`NICE_PYTHON3=/path/to/python3.12` to override the automatic selection.

Or run the following commands:

```shell
mkdir build
cd build
cmake .. -DNICE=ON -DPython3_EXECUTABLE="$(command -v python3.12)"
make -j install
nice_local
```

## Testing the Build

To check a build, compile it with the tests (the default: `./compile_clean.sh` also runs the fast test gate at the end of the build; `--exclude-tests` skips both) and run the unit tests. Execute the following in a separate terminal to ensure the environment variables have been correctly set:

```shell
simple_test_exec test=units
```

It runs every fast unit-test suite in one process and creates its own folder SIMPLE_TEST_units_*date*, which holds the report (simple_test_units_report.txt) and any temporary files. It ends with a summary of the checks; if any check fails, the summary lists it and the program exits with a non-zero status. The folder can be safely removed. If you detect an error, please carefully check the SIMPLE and FFTW installations and the gfortran version. If you still have issues, please file a help ticket on the webpage.

## FFTW Installation

SIMPLE requires three versions of the FFTW library. A double-precision, a single-precision and a threaded-single precision build.

Copy the source from the FFTW homepage (<www.fftw.org>) and unpack it:

```shell
wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.8.tar.gz
tar zxvf fftw-3.3.8.tar.gz
cd fftw-3.3.8
```

By default, the bootstrap script configures and builds the double precision and threaded libraries:

```shell
./bootstrap.sh --prefix=/usr/local
make -j
sudo make install
make distclean
```

Now build the single-precision libraries:

```shell
./bootstrap.sh --prefix=/usr/local --enable-single
make -j
sudo make install
make distclean
```

## New Mac Installation

1. Install homebrew

        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    Add brew to your path

2. Install git, gcc, make, libtiff, jbigkit and cmake. Install Python 3.12 for NICE. SIMPLE builds its pinned OpenBLAS itself.

        brew install git
        brew install cmake
        brew install gcc
        brew install make
        brew install libtiff
        brew install jbigkit
        brew install python@3.12

3. Install fftw and libraries as described above

4. Install SIMPLE

        git clone https://github.com/hael/SIMPLE.git
        cd SIMPLE
        mkdir build
        cd build
        cmake -D NICE=on ..
        make -j install
        cat add2.bashrc >> ~/.zshrc
