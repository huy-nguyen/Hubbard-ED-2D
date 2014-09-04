function totalHamiltonian = hubbardHamiltonian_2D( t, U, Lx, Ly, noOfUp, noOfDn, NUM_CORES )

if (Lx==2) || (Ly==2)
   error('Lx and Ly must be both greater than 2.'); 
   % If we want to deal with the case of Lx = 2 or Ly = 2, put this statement:
   % kinetic = unique( kinetic, 'rows', 'stable'); 
   % before the definition
   % kineticHamiltonian = sparse( ...)
   % at the end of this function
end
noOfSites = Lx * Ly;

totalNoOfPossiblestates = nchoosek( noOfSites, noOfUp) * nchoosek( noOfSites, noOfDn);

%% POTENTIAL HAMILTONIAN:
aux_file_names_potential = {};

fprintf('    Generating the potential Hamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
for i_core =1:NUM_CORES    
aux_file_names_potential{i_core} = strcat('aux_potential_',num2str(noOfSites, '%02d'),...
                                    '_sites_',num2str(noOfUp, '%02d'),...
                                    'u',num2str(noOfDn, '%02d'),...
                                    '_U_',num2str(U, '%4.2f'),...
                                    '_t_',num2str(t),...                                       
                                    '_num_', num2str(i_core, '%02d'),...
                                    ' ',datestr(now,'_yymmdd_HHMMSS'),'.mat');
save(aux_file_names_potential{i_core}, 't', 'U', 'noOfSites', 'noOfUp', 'noOfDn', 'i_core');
end

parfor core_counter_potential=1:NUM_CORES
    fprintf('        Worker %d: Begin.\n', core_counter_potential)
    
    [ combinedBasis_inside_parfor, num_of_states_inside_parfor,dymmy2, dummy3, dummy4, dummy5 ] = generateBasis( noOfSites, noOfUp, noOfDn );
    splitsize = 1 / NUM_CORES * num_of_states_inside_parfor;    
    start_index = floor(round((core_counter_potential-1)*splitsize)) + 1;
    stop_index = floor(round((core_counter_potential)*splitsize));
    j_to_work_on = start_index:stop_index;    
    up_states_to_work_on = combinedBasis_inside_parfor(j_to_work_on, 2);
    dn_states_to_work_on = combinedBasis_inside_parfor(j_to_work_on, 3);        
    results_in_core_loop = zeros( length(j_to_work_on), 3);
    for j = 1:length(j_to_work_on)
        upSectorDec= up_states_to_work_on(j);
        dnSectorDec= dn_states_to_work_on(j);
        upSector = de2bi_modified(upSectorDec, noOfSites);
        dnSector= de2bi_modified(dnSectorDec, noOfSites); 
        doubleOccupancy=bitand(upSector,dnSector); % find all doubly-occupied sites
        potentialEnergy=sum(doubleOccupancy)*U; % sum up the number of doubly occupied sites and multiply by U
        results_in_core_loop(j, 3) = potentialEnergy;
        results_in_core_loop(j, 1) = j_to_work_on(j);
        results_in_core_loop(j, 2) = j_to_work_on(j);
    end
    save_potential_Hamiltonian_segment(aux_file_names_potential{core_counter_potential}, results_in_core_loop);
    
    fprintf('        Worker %d: End.\n', core_counter_potential)
end

fprintf('    Assembling the potential Hamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
potential_sparse_input = zeros( totalNoOfPossiblestates, 3);
last_non_zero_elem_in_potential = 0;
for i_core =1:NUM_CORES 
    current_aux_file_object = matfile( aux_file_names_potential{i_core} );
    [nrows, dummy] = size(current_aux_file_object, 'results_in_core_loop');
    potential_sparse_input( (last_non_zero_elem_in_potential+1): (last_non_zero_elem_in_potential+nrows) , :) ...
                                            = current_aux_file_object.results_in_core_loop;
    last_non_zero_elem_in_potential = last_non_zero_elem_in_potential + nrows;                           
end

potentialHamiltonian = sparse( potential_sparse_input(:, 1), potential_sparse_input(:, 2), potential_sparse_input(:, 3), totalNoOfPossiblestates, totalNoOfPossiblestates);
fprintf('    Done with assembling the potential Hamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
clearvars upSectorDec dnSectorDec upSector dnSector potential_sparse_input;

%% KINETIC HAMILTONIAN:

aux_file_names_kinetic = {};
fprintf('    Generating the kinetic Hamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
for i_core =1:NUM_CORES    
aux_file_names_kinetic{i_core} = strcat('aux_kinetic_',num2str(noOfSites, '%02d'),...
                                    '_sites_',num2str(noOfUp, '%02d'),...
                                    'u',num2str(noOfDn, '%02d'),...
                                    '_U_',num2str(U, '%4.2f'),...
                                    '_t_',num2str(t),...                                       
                                    '_num_', num2str(i_core, '%02d'),...
                                    ' ',datestr(now,'_yymmdd_HHMMSS'),'.mat');
save(aux_file_names_kinetic{i_core}, 't', 'U', 'noOfSites', 'noOfUp', 'noOfDn', 'i_core');
end

max_kinetic_num_non_zero_per_iteration = totalNoOfPossiblestates * 4 * (noOfUp + noOfDn);
actual_num_non_zero_elems_kinetic = 0;
for core_counter_kinetic = 1:NUM_CORES  % will be parfor   
    fprintf('        Worker %d: Begin.\n', core_counter_kinetic)
    [ combinedBasis_inside_parfor, num_of_states_inside_parfor,dummy2, totalNoOfDnStates, upStates, dnStates ] = generateBasis( noOfSites, noOfUp, noOfDn );
    splitsize = 1 / NUM_CORES * num_of_states_inside_parfor;
    start_index = floor(round((core_counter_kinetic-1)*splitsize)) + 1;
    stop_index = floor(round((core_counter_kinetic)*splitsize));
    m_to_work_on = start_index:stop_index;   
    up_states_to_work_on = combinedBasis_inside_parfor( m_to_work_on, 2);
    dn_states_to_work_on = combinedBasis_inside_parfor( m_to_work_on, 3);
    
    KINETIC_COUNTER = 0;
    kinetic_core_rows = zeros(ceil(max_kinetic_num_non_zero_per_iteration / NUM_CORES), 1);
    kinetic_core_cols = zeros( ceil( max_kinetic_num_non_zero_per_iteration / NUM_CORES), 1);
    kinetic_core_elems = zeros( ceil(max_kinetic_num_non_zero_per_iteration / NUM_CORES), 1);
    
    
    for internal_index = 1:length(m_to_work_on)    
        m = m_to_work_on(internal_index);
        % save the unshifted spin up and spin down sectors:
        upSectorDec = up_states_to_work_on(internal_index);
        dnSectorDec = dn_states_to_work_on(internal_index);        
        upSector= de2bi_modified(upSectorDec, noOfSites);
        dnSector= de2bi_modified(dnSectorDec, noOfSites);                  
        % find the occupied lattice sites:    
        upNonZero=find(upSector);
        dnNonZero=find(dnSector);
        % shift for spin up:
        for n=upNonZero % for each occupied site
            % left shift
            leftShiftResult=upSector;
            leftShiftedIndex = hubbard_shift_left( n, Lx, Ly);
            if upSector(leftShiftedIndex)~= 1 % if destination is not occupied
               % perform the shift:
               leftShiftResult(n)=0;
               leftShiftResult(leftShiftedIndex)=1;
                % figure out where in the basis this shifted state is                
               upIndexOfLeftShiftedResult=  binaraysearchasc(upStates, bi2de_modified(leftShiftResult) );
               dnIndexOfLeftShiftedResult=mod( m-1,totalNoOfDnStates)+1;
               basisIndexOfLeftShiftedResult=(upIndexOfLeftShiftedResult-1)*totalNoOfDnStates+dnIndexOfLeftShiftedResult;
               % figure out the sign
               small_index = min(n, leftShiftedIndex) + 1;
               large_index = max(n, leftShiftedIndex) - 1;
               num_electrons_in_between = sum( leftShiftResult( small_index : large_index ) );
               if mod( num_electrons_in_between, 2) % odd # electrons in between: sign change.
                   matrix_element = +t;
               else % even # electrons in between: no sign change
                   matrix_element = -t;
               end
               % assign the matrix element               
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
               kinetic_core_cols(KINETIC_COUNTER) = m;
               kinetic_core_elems(KINETIC_COUNTER) = matrix_element;
            end             
            leftShiftResult=[];leftShiftedIndex=[];upIndexOfLeftShiftedResult=[];dnIndexOfLeftShiftedResult=[];basisIndexOfLeftShiftedResult=[];matrix_element=[];small_index=[];large_index=[];num_electrons_in_between=[];
            
            % right shift            
            rightShiftResult=upSector;
            rightShiftedIndex = hubbard_shift_right( n, Lx, Ly );
            if upSector(rightShiftedIndex)~= 1 % if destination is not occupied
               % perform the shift:
               rightShiftResult(n)=0;
               rightShiftResult(rightShiftedIndex)=1;
                % figure out where in the basis this shifted state is                
               upIndexOfRightShiftedResult=  binaraysearchasc(upStates, bi2de_modified(rightShiftResult) );
               dnIndexOfRightShiftedResult=mod( m-1,totalNoOfDnStates)+1;
               basisIndexOfRightShiftedResult=(upIndexOfRightShiftedResult-1)*totalNoOfDnStates+dnIndexOfRightShiftedResult;
               % figure out the sign
               small_index = min(n, rightShiftedIndex) + 1;
               large_index = max(n, rightShiftedIndex) - 1;
               num_electrons_in_between = sum( rightShiftResult( small_index : large_index ) );
               if mod( num_electrons_in_between, 2) % odd # electrons in between: sign change.
                   matrix_element = +t;
               else % even # electrons in between: no sign change
                   matrix_element = -t;
               end
               % assign the matrix element               
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
               kinetic_core_cols(KINETIC_COUNTER) = m;
               kinetic_core_elems(KINETIC_COUNTER) = matrix_element;
            end
            rightShiftResult=[];rightShiftedIndex=[];upIndexOfRightShiftedResult=[];dnIndexOfRightShiftedResult=[];basisIndexOfRightShiftedResult=[];matrix_element=[];small_index=[];large_index=[];num_electrons_in_between=[];
            
            % up shift          
            upShiftResult=upSector;            
            upShiftedIndex = hubbard_shift_up( n, Lx, Ly);
            if upSector(upShiftedIndex)~= 1 % if destination is not occupied
               % perform the shift:
               upShiftResult(n)=0;
               upShiftResult(upShiftedIndex)=1;
                % figure out where in the basis this shifted state is                
               upIndexOfUpShiftedResult=  binaraysearchasc(upStates, bi2de_modified(upShiftResult) );
               dnIndexOfUpShiftedResult=mod( m-1,totalNoOfDnStates)+1;
               basisIndexOfUpShiftedResult=(upIndexOfUpShiftedResult-1)*totalNoOfDnStates+dnIndexOfUpShiftedResult;
               % figure out the sign
               small_index = min(n, upShiftedIndex) + 1;
               large_index = max(n, upShiftedIndex) - 1;
               num_electrons_in_between = sum( upShiftResult( small_index : large_index ) );
               if mod( num_electrons_in_between, 2) % odd # electrons in between: sign change.
                   matrix_element = +t;
               else % even # electrons in between: no sign change
                   matrix_element = -t;
               end
               % assign the matrix element               
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfUpShiftedResult;
               kinetic_core_cols(KINETIC_COUNTER) = m;
               kinetic_core_elems(KINETIC_COUNTER) = matrix_element;
            end            
            upShiftResult=[];upShiftedIndex=[];upIndexOfUpShiftedResult=[];dnIndexOfUpShiftedResult=[];basisIndexOfUpShiftedResult=[];matrix_element=[];small_index=[];large_index=[];num_electrons_in_between=[];
            
            % down shift            
            downShiftResult=upSector;
            downShiftedIndex = hubbard_shift_down( n, Lx, Ly );
            if upSector(downShiftedIndex)~= 1 % if destination is not occupied
               % perform the shift:
               downShiftResult(n)=0;
               downShiftResult(downShiftedIndex)=1;
                % figure out where in the basis this shifted state is                
               upIndexOfDownShiftedResult=  binaraysearchasc(upStates, bi2de_modified(downShiftResult) );
               dnIndexOfDownShiftedResult=mod( m-1,totalNoOfDnStates)+1;
               basisIndexOfDownShiftedResult=(upIndexOfDownShiftedResult-1)*totalNoOfDnStates+dnIndexOfDownShiftedResult;
               % figure out the sign
               small_index = min(n, downShiftedIndex) + 1;
               large_index = max(n, downShiftedIndex) - 1;
               num_electrons_in_between = sum( downShiftResult( small_index : large_index ) );
               if mod( num_electrons_in_between, 2) % odd # electrons in between: sign change.
                   matrix_element = +t;
               else % even # electrons in between: no sign change
                   matrix_element = -t;
               end
               % assign the matrix element               
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfDownShiftedResult;
               kinetic_core_cols(KINETIC_COUNTER) = m;
               kinetic_core_elems(KINETIC_COUNTER) = matrix_element;
            end            
            downShiftResult=[];downShiftedIndex=[];upIndexOfDownShiftedResult=[];dnIndexOfDownShiftedResult=[];basisIndexOfDownShiftedResult=[];matrix_element=[];small_index=[];large_index=[];num_electrons_in_between=[];
        end
        n = [];
        % shift for spin down:
        for n=dnNonZero % for each occupied site
            % left shift
            leftShiftResult=dnSector;
            leftShiftedIndex = hubbard_shift_left( n, Lx, Ly);
            if dnSector(leftShiftedIndex)~= 1 % if destination is not occupied
               % perform the shift:
               leftShiftResult(n)=0;
               leftShiftResult(leftShiftedIndex)=1;
                % figure out where in the basis this shifted state is   
               dnIndexOfLeftShiftedResult=  binaraysearchasc(dnStates, bi2de_modified(leftShiftResult) );
               upIndexOfLeftShiftedResult= floor(( m - 1 )/totalNoOfDnStates)+1;
               basisIndexOfLeftShiftedResult=(upIndexOfLeftShiftedResult-1)*totalNoOfDnStates+dnIndexOfLeftShiftedResult;  
               small_index = min(n, leftShiftedIndex) + 1;
               large_index = max(n, leftShiftedIndex) - 1;
               num_electrons_in_between = sum( leftShiftResult( small_index : large_index ) );
               if mod( num_electrons_in_between, 2) % odd # electrons in between: sign change.
                   matrix_element = +t;
               else % even # electrons in between: no sign change
                   matrix_element = -t;
               end
               % assign the matrix element               
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfLeftShiftedResult;
               kinetic_core_cols(KINETIC_COUNTER) = m;
               kinetic_core_elems(KINETIC_COUNTER) = matrix_element;
            end             
            leftShiftResult=[];leftShiftedIndex=[];upIndexOfLeftShiftedResult=[];dnIndexOfLeftShiftedResult=[];basisIndexOfLeftShiftedResult=[];matrix_element=[];small_index=[];large_index=[];num_electrons_in_between=[];
            
            % right shift
            rightShiftResult=dnSector;
            rightShiftedIndex = hubbard_shift_right( n, Lx, Ly);
            if dnSector(rightShiftedIndex)~= 1 % if destination is not occupied
               % perform the shift:
               rightShiftResult(n)=0;
               rightShiftResult(rightShiftedIndex)=1;
                % figure out where in the basis this shifted state is   
               dnIndexOfRightShiftedResult=  binaraysearchasc(dnStates, bi2de_modified(rightShiftResult) );
               upIndexOfRightShiftedResult= floor(( m - 1 )/totalNoOfDnStates)+1;
               basisIndexOfRightShiftedResult=(upIndexOfRightShiftedResult-1)*totalNoOfDnStates+dnIndexOfRightShiftedResult;  
               small_index = min(n, rightShiftedIndex) + 1;
               large_index = max(n, rightShiftedIndex) - 1;
               num_electrons_in_between = sum( rightShiftResult( small_index : large_index ) );
               if mod( num_electrons_in_between, 2) % odd # electrons in between: sign change.
                   matrix_element = +t;
               else % even # electrons in between: no sign change
                   matrix_element = -t;
               end
               % assign the matrix element               
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfRightShiftedResult;
               kinetic_core_cols(KINETIC_COUNTER) = m;
               kinetic_core_elems(KINETIC_COUNTER) = matrix_element;
            end             
            rightShiftResult=[];rightShiftedIndex=[];upIndexOfRightShiftedResult=[];dnIndexOfRightShiftedResult=[];basisIndexOfRightShiftedResult=[];matrix_element=[];small_index=[];large_index=[];num_electrons_in_between=[];         
            
            % up shift
            upShiftResult=dnSector;     
            upShiftedIndex = hubbard_shift_up( n, Lx, Ly);
            if dnSector(upShiftedIndex)~= 1 % if destination is not occupied
               % perform the shift:
               upShiftResult(n)=0;
               upShiftResult(upShiftedIndex)=1;
                % figure out where in the basis this shifted state is   
               dnIndexOfUpShiftedResult=  binaraysearchasc(dnStates, bi2de_modified(upShiftResult) );
               upIndexOfUpShiftedResult= floor(( m - 1 )/totalNoOfDnStates)+1;
               basisIndexOfUpShiftedResult=(upIndexOfUpShiftedResult-1)*totalNoOfDnStates+dnIndexOfUpShiftedResult;  
               small_index = min(n, upShiftedIndex) + 1;
               large_index = max(n, upShiftedIndex) - 1;
               num_electrons_in_between = sum( upShiftResult( small_index : large_index ) );
               if mod( num_electrons_in_between, 2) % odd # electrons in between: sign change.
                   matrix_element = +t;
               else % even # electrons in between: no sign change
                   matrix_element = -t;
               end
               % assign the matrix element               
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfUpShiftedResult;
               kinetic_core_cols(KINETIC_COUNTER) = m;
               kinetic_core_elems(KINETIC_COUNTER) = matrix_element;
            end             
            upShiftResult=[];upShiftedIndex=[];upIndexOfUpShiftedResult=[];dnIndexOfUpShiftedResult=[];basisIndexOfUpShiftedResult=[];matrix_element=[];small_index=[];large_index=[];num_electrons_in_between=[];         
            
            
            % down shift
            downShiftResult=dnSector;
            downShiftedIndex = hubbard_shift_down( n, Lx, Ly); 
            if dnSector(downShiftedIndex)~= 1 % if destination is not occupied
               % perform the shift:
               downShiftResult(n)=0;
               downShiftResult(downShiftedIndex)=1;
                % figure out where in the basis this shifted state is   
               dnIndexOfDownShiftedResult=  binaraysearchasc(dnStates, bi2de_modified(downShiftResult) );
               upIndexOfDownShiftedResult= floor(( m - 1 )/totalNoOfDnStates)+1;
               basisIndexOfDownShiftedResult=(upIndexOfDownShiftedResult-1)*totalNoOfDnStates+dnIndexOfDownShiftedResult;  
               small_index = min(n, downShiftedIndex) + 1;
               large_index = max(n, downShiftedIndex) - 1;
               num_electrons_in_between = sum( downShiftResult( small_index : large_index ) );
               if mod( num_electrons_in_between, 2) % odd # electrons in between: sign change.
                   matrix_element = +t;
               else % even # electrons in between: no sign change
                   matrix_element = -t;
               end
               % assign the matrix element
               KINETIC_COUNTER = KINETIC_COUNTER + 1;
               kinetic_core_rows(KINETIC_COUNTER) = basisIndexOfDownShiftedResult;
               kinetic_core_cols(KINETIC_COUNTER) = m;
               kinetic_core_elems(KINETIC_COUNTER) = matrix_element;
            end             
            downShiftResult=[];downShiftedIndex=[];upIndexOfDownShiftedResult=[];dnIndexOfDownShiftedResult=[];basisIndexOfDownShiftedResult=[];matrix_element=[];small_index=[];large_index=[];num_electrons_in_between=[];
            
            
        end
    end
    kinetic_core_rows = kinetic_core_rows( kinetic_core_rows ~= 0);
    kinetic_core_cols = kinetic_core_cols( 1:length(kinetic_core_rows));
    kinetic_core_elems = kinetic_core_elems( 1:length(kinetic_core_rows));    
    kinetic_per_core = horzcat(kinetic_core_rows, kinetic_core_cols, kinetic_core_elems); 
    save_kinetic_Hamiltonian_segment( aux_file_names_kinetic{core_counter_kinetic}, kinetic_per_core);
    actual_num_non_zero_elems_kinetic = actual_num_non_zero_elems_kinetic + length(kinetic_core_rows);
    fprintf('        Worker %d: End.\n', core_counter_kinetic)
end


fprintf('    Assembling the kinetic Hamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
kinetic = zeros( actual_num_non_zero_elems_kinetic, 3);
last_non_zero_elem_in_kinetic = 0;
for i_core =1:NUM_CORES 
    current_aux_file_object = matfile( aux_file_names_kinetic{i_core} );
    [nrows, dummy] = size(current_aux_file_object, 'kinetic_per_core');
    kinetic( (last_non_zero_elem_in_kinetic+1): (last_non_zero_elem_in_kinetic+nrows) , :) ...
                                            = current_aux_file_object.kinetic_per_core;
    last_non_zero_elem_in_kinetic = last_non_zero_elem_in_kinetic + nrows;
end
    
kineticHamiltonian = sparse( kinetic(:,1), kinetic(:,2), kinetic(:,3), totalNoOfPossiblestates, totalNoOfPossiblestates);
fprintf('    Done with assembling the kinetic Hamiltonian at time %s.\n', datestr(now,'yymmdd_HHMMSS'))
%% TOTAL HAMILTONIAN:
totalHamiltonian=kineticHamiltonian+potentialHamiltonian;

for i_core =1:NUM_CORES    
    delete( aux_file_names_potential{i_core}, aux_file_names_kinetic{i_core});
end

end

function save_kinetic_Hamiltonian_segment(aux_file_name, kinetic_per_core)
save(aux_file_name, '-append', 'kinetic_per_core', '-v7.3');
end

function save_potential_Hamiltonian_segment(aux_file_name, results_in_core_loop)
save(aux_file_name, '-append', 'results_in_core_loop', '-v7.3');
end
