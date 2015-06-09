This file is a readme for a shared `libStraticounter` library that can be used to interface with the scibox software.

# Prerequisites

1. MATLAB Runtime v8.4 (for R2014b) - [Found Here](http://www.mathworks.com/products/compiler/mcr/)
2. Your compiler (g++) must be version 4.7x (may involve using update-alternatives to switch default version on your computer)

# Building the Library

To create a shared library, run the following command from the root of StratiCounter:

```
mcc -W cpplib:libStraticounter -T link:lib  -d MCR_Library/ -I Subroutines/ -I Subroutines/algorithm/ -I Subroutines/batchresults/ -I Subroutines/checkinput/ -I Subroutines/layertemplates/ -I Subroutines/matchmaker/ -I Subroutines/preprocessing/ -I Subroutines/showresults/ -I Subroutines/syntheticdata/ -I Subroutines/utilities/  straticounter_scibox.m
```

This will create a shared object file along with a C++ source and C++ header in the `MCR_Library/` directory.

# Approach

A small interface function called `straticounter_scibox.m` was written to interface with the `straticounter.m` function. This is an attempt to decouple the interface of scibox from that of the MATLAB tool. The main funciton `straticounter.m` needed to be modified in two key ways:

1. Some code in the file should not run when used as a library (like plots) - Specifically for plots it turns out the easiest/best way to address this was to force the Runtype.plotlevel variable to be zero when the code is deployed. In general the `isdeployed` variable was used to make these conditions. Example:

    ```
    if ~isdeployed
        addpath(genpath('./Subroutines'))
        addpath(genpath('./Settings'))
    end
    ```

2. The interface to scibox requires the ability to specify an arbitrary output directory - This was addressed by extending the `Runtype` structure to include a variable called `outdir` and modifying the `makeoutputfolder.m` to use the `Runtype.outdir` variable for the path of the directory if it is not empty.
3. In order to accomplish items 1 and 2 the function `straticounter.m` now takes a variable number of inputs. When given 1 parameter the function works the same as before. When called with 2 parameters (like `straticounter_scibox.m` does) the inputs are a path to the data file to load and the path for the desired output directory.

#Terminal interface for testing

In order to create an executable (that can be run from the terminal):

```mbuild MCR_Library/Straticounter.cpp MCR_Library/libStraticounter.so```
