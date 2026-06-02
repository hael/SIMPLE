# Install

This page summarizes the build and setup information from the repository
README.

## System Requirements

- Linux, with Ubuntu 16.04 or newer recommended
- macOS 10.10 or newer
- CMake 3.5 or newer
- GNU `gcc` and `gfortran` 14.2 or newer
- FFTW 3.3 or newer, with double, single, and threaded builds
- libcurl 7 or newer
- libTIFF 4 or newer
- jbigkit 2 or newer, optional
- Python 3.10 or newer

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
Or use run the following commands:
```shell
mkdir build
cd build
cmake .. -D NICE=on
make -j install
nice_local
```

## Testing the Build

To ensure that SIMPLE has been correctly installed, we recommend running the application simple_test_install. It will perform elementary tests of the base components in the SIMPLE library. Execute the following in a separate terminal to ensure the environment variables have been correctly set:
```shell
simple_test_install
```
The program will create its own folder SIMPLE_TEST_INSTALL*date* where temporary files are stored. Upon succesful completion you should see

    $ **** SIMPLE_TEST_INSTALL NORMAL STOP ****

simple_test_install can be executed anywhere and the folder created can be safely removed. If any of the individual tests fail an error message will be displayed. If you detect an error, please carefully check the SIMPLE and FFTW installations and the gfortran version. If you still have issues, please file a help ticket on the webpage.

## FFTW Installation

SIMPLE requires three versions of the FFTW library. A double-precision, a single-precision and a threaded-single precision build.

Copy the source from the FFTW homepage (www.fftw.org) and unpack it:
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
1.  Install homebrew

        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    Add brew to your path

3.  Install the dependencies git, gcc, libtiff, jbigkit and cmake – make sure you have a python version higher than 3.10

        brew install git
        brew install cmake
        brew install gcc
        brew install libtiff
        brew install jbigkit
        brew install python@3.10

3.  Install fftw and libraries as described above

4.  Install SIMPLE

        git clone https://github.com/hael/SIMPLE.git
        cd SIMPLE
        mkdir build
        cd build
        cmake -D NICE=on ..
        make -j install
        cat add2.bashrc >> ~/.zshrc
