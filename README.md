
This set of files contain a toolbox in multiple languages, including c++, fortran 77, 
fortran 90, mathematica, matlab (and possibly others) to play with random matrices, 
quantum computation and tensor networks. 


To be able to use de c++ toolbox, one must install tclap and itpp. In ubuntu it corresponds to the packages libitpp-dev libtclap-dev


Maxima -  Files provided by Thomas Gorin to calculate integral of polynomials of uniatry matrices with the Haar measure
cpp - C++ toolbox
matlab_libraries - Matlab toolbox
modules - Fortran 90 toolbox
.gitignore - configuration file
Carlos.m - Various tools to plot an process information in mathematica
Quantum.m - Quantum information tools
README - This file
SpinChain.m - Spin chain files
TheoreticalRMT.m - Files to plot RMT results in mathematica
test.nb - File that allows testing of the mathematica libraries. 


## Installation on Linux-based Systems

1. Clone or download the [`libs`](https://github.com/carlospgmat03/libs) repository into your home folder:

    ```bash
    git clone https://github.com/carlospgmat03/libs ~/libs
    ```
    or if you prefer, download the specific [v3.1.0](https://github.com/carlospgmat03/libs/releases/tag/v3.1.0) version.

2. Add the following line to your `~/.Wolfram/Kernel/init.m` file:

    ```mathematica
    AppendTo[$Path, FileNameJoin[{"/home/" <> ToString[$Username], "libs"}]];
    ```

> **Note**: If `init.m` does not exist, create the file and the necessary directory path.
