clear all
close all


% ====== system setup ============
alpha = 0.05;
do_plot = true;
ion_soft = 1.0;                       % soft constant for v_ext potential
nspin  = 1;                           % number of spin
tsmear = 0.01/27.2114;                % according hartree-fock assumptions
natom   = 12;
atom_Z  = ones(natom,1)*1.2;          % 12 H atoms
atom_Z_real  = ones(natom,1)*1.0;     % 12 H atoms
q_total   = natom*1; % total electron number in system
nbox = 6;



% === define alternating H chain ===
d = 2.0;
delta = 0.0;



coord = zeros(natom,1);

%Set up first coordinate as d/2 + delta / 2
coord(1) = d/2+delta/2;


% shift coordinates according to the bond and distance between atoms; d and delta/2 for odd and even number
%atoms
for i=2:natom
    %Odd numbered atoms
    if (mod(i,2)==1)
        coord(i)= (i-1) * d + delta/2 + d/2;
        
    %Even numbered atoms
    else
        coord(i)=(i-1) * d - delta/2 + d/2;
    end
end


%>>>>>>>>  get box/cell length <<<<<<<<<<<<<<<<<<<
box_len = coord(end)+d/2+delta/2;  % box starts at zero, in bohr
ngrid   = floor(box_len*4);          % include the last point length of cell
norb    = ngrid;                    % number of orbitals to solve.


ksdft

