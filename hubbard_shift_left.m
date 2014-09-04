function index_out = hubbard_shift_left( index_in, Lx, Ly )

% switch to base 0
row = mod( (index_in-1), Ly );
col = floor( (index_in-1)/Ly );

col = mod( col - 1, Lx );

% now back to base 1:
index_out = col * Ly + row + 1;
end