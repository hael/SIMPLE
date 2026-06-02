# Install FFTW

## Mac
```shell
sudo fink install fftw3 fftw3-shlibs
```
Then set link by: `-L/sw/lib/ -lfftw3f -lfftw3f_threads`

# Linux
```shell
make clean
./configure prefix=/path/to/fftw/ --enable-float --enable-threads --with-openmp
make
make install
```
Then set link by: `-L/path/to/fftw/ -lfftw3f -lfftw3f_threads`
