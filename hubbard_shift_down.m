function index_out = hubbard_shift_down( index_in, Lx, Ly )

% switch to base 0
row = mod( (index_in-1), Ly );
col = floor( (index_in-1)/Ly );

row = mod( row+1, Ly);

% now back to base 1:
index_out = col * Ly + row + 1;
end