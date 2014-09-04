NUM_CORE = 4;

matlabpool('open', NUM_CORE);
t = 1;
U = 4;
Lx = 3;
Ly = 3;
noOfUp = 2;
noOfDn = 2;

file_name = strcat('Test_',datestr(now,'_yymmdd_HHMMSS'),'.mat')

[ totalHamiltonian, kineticHamiltonian,  potentialHamiltonian] = hubbardHamiltonian_2D( t, U, Lx, Ly, noOfUp, noOfDn, NUM_CORE );

[v, d] = eigs(totalHamiltonian, 1, 'sa');

save(file_name, 'v', 'd', '-v7.3');

matlabpool('close');