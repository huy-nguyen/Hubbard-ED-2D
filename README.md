# Exact Diagonalization of the Hubbard Model in 2-D

This repository contains the `MATLAB` code to perform exact calculations of the *imaginary-time correlation functions* of the *Hubbard* model in two dimensions.

The Hubbard model is widely believed to be the model that describes [high-temperature superconductivity](https://en.wikipedia.org/wiki/High-temperature_superconductivity). I presented the theory behind this model in a manner accessible to senior-year physics majors in Chapter 2 of my [undergraduate thesis](https://www.huynguyen.io/Reed-thesis).

This repo is a branch and an extension of the [one-dimensional version](https://github.com/huy-nguyen/Hubbard-ED-1D).

The code is designed to run in parallel within a single compute node on multi-core Linux high-performance clusters using `MATLAB`'s [Parallel Computing Toolbox](http://www.mathworks.com/products/parallel-computing/). Furthermore, `MATLAB` will automatically distributes the work load over multiple compute nodes if you have a license for their [Distributed Computing Server](http://www.mathworks.com/products/distriben/).

Please keep in mind that the computational cost of this code increases very quickly. The largest system I've been able to run is a 4x4 square lattice with 5 spin-up and 5 spin-down electrons. This calculation took more than a day on an 8-core compute node with 192GB of memory. Peak memory usage was about 180GB.

To run the code:

1. Modify the parameters of the model in `sample.m`
2. Open the Linux command line and `cd` to the directory containing that `m` file.
3. Execute this command `matlab -nodisplay -nosplash -nojvm <sample.m> sample.out`.

The output is a series of `mat` files, one for each value of `tau`. Each file contains 2 square matrices whose entries are the values of the correlation functions for the spin-up and spin-down sectors for that particular value of `tau`.

This repository includes a number of unit tests. All files whose names start with *test_* are unit test files that can be run run with [xUnit Test Framework](http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework/content/matlab_xunit_3_1_1/xunit/assertEqual.m).


