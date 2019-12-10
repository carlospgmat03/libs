# qc-cpp
C++ toolbox for quantum computation

## Getting Started

### Prerequisites
* [GNU gcc compiler][gcc]
* [Clang][clang]
* [GNU fortran compiler][fortran]
* [TCLAP][tclap]
* [IT++][IT++]
* [FFTW][fftw]
* [Doxygen][doxygen]
* [GNU Make][make]


### Ubuntu Build Guide
To install the dependencies, run the following commands.

```bash
sudo apt update
```

GNU gcc compiler

```bash
sudo apt install g++ clang ibc++-helpers pentium-builder
```
GNU fortran compiler

```bash
sudo apt install gfortran
```

IT++

```bash
sudo apt install libblas-dev librtlsdr-dev libncurses5-dev libitpp-dev -y
```

Doxygen

```bash
sudo apt install doxygen graphviz
```
[gcc]: https://www.gnu.org/software/gcc/
[fortran]: https://gcc.gnu.org/fortran/
[clang]: https://clang.llvm.org/
[tclap]: https://sourceforge.net/projects/tclap/
[IT++]: https://sourceforge.net/projects/itpp/
[fftw]: http://fftw.org/
[doxygen]: http://www.doxygen.org/
[make]: https://www.gnu.org/software/make/
