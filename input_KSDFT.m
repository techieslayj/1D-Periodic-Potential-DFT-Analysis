clear all
close all


% ====== system setup ============
do_plot = true;
ion_soft = 1.0;                       % soft constant for v_ext potential
nspin  = 1;                           % number of spin
tsmear = 0.01/27.2114;                % according hartree-fock assumptions
natom   = 12;
atom_Z  = ones(natom,1)*1.2;          % 14 H atoms
atom_Z_real  = ones(natom,1)*1.0;     % 14 H atoms
q_total   = natom*1;                  % total electron number in system


% === define alternating H chain ===
eq_dist = 2.2;
bond = 1.8;   % equilibrium distance (for LDA energy) for eq_dist=2.2
dist = (eq_dist*(natom-1) - bond*natom/2)/(natom/2-1);
coord = zeros(natom,1);
coord(1) = 0.5;
% shift coordinates according to the bond and distance between atoms
for i=1:natom-1
    if (mod(i,2)==1)
        coord(i+1)=coord(i) + bond;
    else
        coord(i+1)=coord(i) + dist;
    end
end


%>>>>>>>>  get box/cell length <<<<<<<<<<<<<<<<<<<
box_len = coord(end)+0.5;            % box starts at zero, in bohr
ngrid   = floor(box_len*3)          % include the last point length of cell
norb    = ngrid;                    % number of orbitals to solve.


ksdft

