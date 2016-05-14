% The number of cores for parallel running:
NUM_CORE = 4;

matlabpool('open', NUM_CORE);

% Physical configuration of the Hubbard model:
t = 1;
U = 4;
Lx = 3;
Ly = 3;
noOfUp = 2;
noOfDn = 2;

% The code will calculate the imaginary-time correlations for all tau values
% between tau_start and tau_end, in steps of tau_step:
tau_start = 1;
tau_end = 3;
tau_step = 1;

% How many eigenvalues to request from MATLAB's sparse matrix diagonalization.
% The higher the number eigenvalues requtested, the longer the code takes to run.
NUM_OF_EIGEN_VALUES = 17;

% Indicate whether we want the correlation functions for the spin up or spin down sector:
sector = 'up';

% Do not change these during production run:
method = 'long_tau';
commit_number = 'testtesttest';
need_profiling = 'No';

list_of_generated_files = unequalTimeGF_long_tau_parallel_2D( t, U, tau_start, tau_end, tau_step, Lx, Ly, noOfUp, noOfDn, NUM_OF_EIGEN_VALUES, sector, method, commit_number, need_profiling, NUM_CORES );
