t = 1;
U = 4;
Lx = 3;
Ly = 3;
noOfUp = 5;
noOfDn = 5;
NUM_CORES = 4;

hamiltonian = hubbardHamiltonian_2D( t, U, Lx, Ly, noOfUp, noOfDn, NUM_CORES );

save('4x4_6u5d_Hamiltonian.mat', 'firstHamiltonian', '-v7.3');