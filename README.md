           _______ _____ _______  _____         _______
          |______   |   |  |  | |_____] |      |______
          ______| __|__ |  |  | |       |_____ |______  v3b1

# SIMPLE

![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey)
![Build](https://img.shields.io/badge/build-CMake-blue)
![Compiler](https://img.shields.io/badge/compiler-gcc%20%7C%20gfortran-orange)
![FFTW](https://img.shields.io/badge/dependency-FFTW%203.x-green)
![Version](https://img.shields.io/badge/version-3.0.0-brightgreen)

**SIMPLE** is a high‑performance, end‑to‑end cryo‑EM image processing
and reconstruction platform.

It transforms raw electron microscopy data --- movies, micrographs, and
particle images --- into scientifically meaningful outputs including:

-   2D class averages
-   High-resolution 3D reconstructions
-   Symmetry analysis
-   Atomic-level structural models

SIMPLE supports interactive batch workflows and fully automated
streaming pipelines, scaling seamlessly from single workstations to
distributed HPC environments.

------------------------------------------------------------------------

## 🏗️  Build Status
[![Build SIMPLE & Run Tests (Linux/MacOS)](https://github.com/hael/SIMPLE/actions/workflows/ci_build_and_test.yml/badge.svg)](https://github.com/hael/SIMPLE/actions/workflows/ci_build_and_test.yml) 

------------------------------------------------------------------------

## 🚀 Release

Stable release (v3.0.0):

https://github.com/hael/SIMPLE/releases/tag/v3.0.0

Public unit tests:

https://zenodo.org/records/18789533

------------------------------------------------------------------------

# 🧬 Features

-   End-to-end cryo‑EM workflow
-   Streaming data processing
-   HPC scalability
-   Modular architecture
-   Test-driven development
-   Production-ready scientific computing environment

------------------------------------------------------------------------

# 📦 Installation

## System Requirements

-   Linux (Ubuntu 16.04+ recommended)
-   macOS 10.10+
-   CMake ≥ 3.5
-   GNU gcc & gfortran ≥ 14.2
-   FFTW ≥ 3.3 (double + single + threaded builds)
-   libcurl ≥ 7 
-   libTIFF ≥ 4 
-   jbigkit ≥ 2 (optional)
-   Python ≥ 3.10

------------------------------------------------------------------------

## Standard Build

### 1️⃣ Obtain source

Stable release:

``` bash
gunzip SIMPLE-3.0.0.tar.gz
tar -xvf SIMPLE-3.0.0.tar
cd SIMPLE-3.0.0
```

Or clone:

``` bash
git clone https://github.com/hael/SIMPLE.git
cd SIMPLE
```

------------------------------------------------------------------------

### 2️⃣ Build

``` bash
mkdir build
cd build
cmake ..
make -j install
```

------------------------------------------------------------------------

### 3️⃣ Environment Setup

For bash:

``` bash
cat add2.bashrc >> ~/.bashrc
```

For tcsh/csh:

``` bash
cat add2.tcshrc >> ~/.tcshrc
```

Minimal manual setup:

Add to `PATH`:

    <simple_path>/build/bin
    <simple_path>/build/scripts

Define:

    SIMPLE_PATH=<simple_path>

------------------------------------------------------------------------

# ⚙ Advanced Configuration

Custom install directory:

``` bash
cmake -DCMAKE_INSTALL_PREFIX=<install_dir> ..
make -j install
```

Custom compiler / FFTW paths:

``` bash
FC=<gfortran_path> CC=<gcc_path> FFTW_DIR=<fftw_path> cmake ..
make -j install
```

------------------------------------------------------------------------

# 🔁 Updating

``` bash
cd ~/SIMPLE
git stash
git stash clear
git pull
cd build
make -j install
source add2.bashrc
nice_local
```

------------------------------------------------------------------------

# 📜 License

SIMPLE is distributed under the **GNU General Public License v3
(GPLv3)** or later.

This software is provided **without warranty**, including without
implied warranty of merchantability or fitness for a particular purpose.

See the GNU General Public License for details.

------------------------------------------------------------------------

# 🧠 Philosophy

SIMPLE combines:

-   High-performance numerical kernels
-   Cryo‑EM domain abstractions
-   Workflow orchestration
-   Engineering discipline

into a scalable scientific application suite built for real experimental
workflows.
