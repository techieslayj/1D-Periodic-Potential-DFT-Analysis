%clear all
%close all

% ps = parallel.Settings;
% ps.Pool.AutoCreate = false;

fprintf('\n    >>> KS-DFT <<< \n\n');

%>>>>>>>>>>>>>>>>>> parameters below should not be modified <<<<<<<<<<<<<<<<<

%--------------------------------------------------------------
% Soft Coulumb potential is from
% "Reference electronic structure calculations in one dimension"
% Phys. Chem. Chem. Phys., 2012, 14, 8581–8590
%--------------------------------------------------------------

% check
if length(coord)~=length(atom_Z)
    fprintf('length(coord)~=length(atom_Z)\n');
    stop
end

% set up grid.
natom = length(atom_Z);
h = box_len/(ngrid-1);
x = zeros(ngrid,1);
for i=1:ngrid
    x(i) = (i-1)*h;
end

% fourier transformation z/r and fourier transformation of 1/r is 1/q^2
% V atom = z / r + delta(v)
%=========== currenty creating an algorithm to account for 1D system external
% potential==================
% because if use current system with fourier transformation will
% cause system to diverge because since we are not dealing with a 3D system
% r^2 is absent and the denominator will cause system to diverge.

[vext]  = make_external_pot(ion_soft,x,natom,atom_Z,coord,ngrid);
[e_ion] = make_ionic_energy(natom, coord, atom_Z);


%================
%     KS-DFT
%================
print_occ = false;
[ ee,ev,rho,occ,vks,mu,vxc,vhart ] = scf( x, vext, q_total, norb, ...
    box_len, tsmear, e_ion, print_occ );



% count the number of occupied KS orbitals
noccupied = 0;
for i=1:ngrid
    if abs(occ(i))>1e-4
        noccupied = noccupied + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% plots results of DFT 
%%%%%%%%%%%%%%%%%%%%%%%%%
if do_plot
    figure;
    plot(x,vext);
    hold on; plot(x,vks);
    hold on; plot(x,rho);
    tmp2 = zeros(size(coord));
    plot(coord,tmp2,'ro');
    legend('external potential','KS potential','density');
    fprintf('press any key to continue ....\n');
    pause;
else
    fprintf('writing eps files ...\n');
    f=figure('Visible','off');
    plot(x,vext);
    hold on; plot(x,vks);
    hold on; plot(x,rho);
    tmp2 = zeros(size(coord));
    plot(coord,tmp2,'ro');
    legend('external potential','KS potential','density');
    print('-depsc2','-r300','ks_dft.eps');
    
    % canonical orbitals ...
    f=figure('Visible','off');
    for i=1:noccupied
        hold on;
        plot(x,ev(:,i));
    end
    tmp2 = zeros(size(coord));
    plot(coord,tmp2,'rx');
    title('MOs')
    print('-depsc2','-r300','MOs.eps');
end


return




%%%%%%%%%%% below are not related to KS-DFT procedure %%%%%%%%%%%%%%










% =============================
%   make total density matrix
% =============================
dm_tot = zeros(ngrid);
for ib=1:norb
    dm_tot = dm_tot + occ(ib)*ev(:,ib)*ev(:,ib)';
end
ref_rho = rho;
ref_vks = vks;


%%
%%%%%% localization using boys method %%%%%

norb_localize = 4;
[outwfr] = localization_MOs(x,ev(:,1:norb_localize));
%
figure
sum_rhor = 0.0;
for q=1:norb_localize
    if q==1 || q==3
        continue
    end
    oo = outwfr(:,q).^2*2;
    sum_rhor = sum_rhor+oo;
    hold on
    plot(x,oo);
end
plot(x,sum_rhor,'-k');
hold on
plot(x,rho,'ro');

fprintf('continue DM partitioning?\n');
pause;

omega_max = 10;
nfreq = 20;
rpa = compute_rpa_energy(omega_max,nfreq,x,ev,ee,occ);
exx = calc_exx_energy(ev,ee,occ,x)


% >>>>>>>>>> DM partitioning related <<<<<<<<

comm_chempot  = true;                % equlibrate cluster and env chemical potential?
q_clu         = 4.0;                 % electron number in cluster
q_env         = q_total-q_clu;  % electron number in env.
atom_list_clu = [1 2 3 4];                 % atoms grouped to form cluster
atom_list_env = [4 5 6 7 8];         % atoms grouped to form environment

% % dmfet can partition
% atom_Z = ones(1,14)*2;               % six H atoms
% dx=2.0;
% comm_chempot = false;                % equlibrate cluster and env chemical potential?
% q_clu   = 7.0;                       % electron number in cluster
% atom_list_clu = [1:7];               % atoms grouped to form cluster
% atom_list_env = [8:14];              % atoms grouped to form environment
% q_total = 14;                        % total electron number in system

% %%% high temperature mimics the metal system, fractional occupation.
% tsmear = 0.5/27.2114;                % in hartree
% ngrid = 201;                         % include the last point
% box_len = 60;                        % box stars at zero, in bohr
% atom_Z = ones(1,14)*2;               % six H atoms
% dx=2.5;
% coord = [10:dx:10+dx*13];            % cooridates
% norb = 30;                           % number of orbitals to solve.
% tol_ksdft = 1e-4;                    % tolerance for stopping KSDFT
% pen_coeff = 1e-4;                    % pen_coeff for regularizing vemb.
% q_total = 24;                        % total electron number in system
% atom_list_clu = [1 2 3 4 5 6 7];               % atoms grouped to form cluster
% atom_list_env = [6 7 8 9 10 11 12 13 14];


%%%%%%%%%%%%%%%%%%%%%%%%%%
%  density partitioning
%%%%%%%%%%%%%%%%%%%%%%%%%%

density_partition;
% plot
figure
plot(x,ref_rho);
hold on;
plot(x,rho_clu);
plot(x,rho_env);
plot(x,rho_clu+rho_env,'ro');
plot(x,vemb);
title ('density partitioning and embedding potential')
stop






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% density matrix partitioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\ndensity matrix partitioning...\n\n')
fprintf('box-length: %f \nngrid: %d  \ndx: %f\n',box_len,ngrid,x(2)-x(1));
fprintf('q_tot: %f  \nq_clu: %f  \nq_env: %f\n\n',q_total,q_clu,q_env);

% Initialize KS potential for cluster and env.
vext_clu = zeros(length(vext),1);
vext_env = zeros(length(vext),1);
for j=1:ngrid
    for i=1:natom
        dist = coord(i) - x(j);
        pot = -atom_Z(i)/sqrt(dist*dist+1.0);
        location = find(atom_list_clu==i);
        if  length(location)==1
            % this atom belongs to cluster
            vext_clu(j) = vext_clu(j) + pot;
        else
            % this atom belong to env.
            vext_env(j) = vext_env(j) + pot;
        end
    end
end
vks_clu = vext_clu;
vks_env = vext_env;

% f_handler = @(vemb)wrap_W_dm(vemb,nspin,x,comm_chempot,box_len,ev_len,ngrid,norb,...
%         tsmear,q_clu,q_env,vks_clu,vks_env,dm_tot);
% options = optimset('Display','iter');
% [vemb2,fval,exitflag,output] = lbfgs(f_handler, reshape(vemb,1,ngrid^2), options);
% vemb = reshape(vemb,ngrid,ngrid);

vemb = zeros( ngrid );
tsmear = 0.2/27.2114;

for iter=1:10
    
    % ----- BFGS minimize (-W) functional -----
    fprintf('solve for vemb (L-BFGS) ...\n\n');
    f_handler = @(vemb)(wrap_W_dm(vemb,nspin,x,comm_chempot,box_len,ngrid,norb,...
        tsmear,q_clu,q_env,vks_clu,vks_env,dm_tot));
    options = optimset('GradObj','on','Display','Iter','MaxIter',10);
    [vemb,final_W] = fminlbfgs(f_handler,reshape(vemb,1,ngrid^2),options);
    vemb = reshape(vemb,ngrid,ngrid);
    
    [W,grad,rho_clu,rho_env]=W_dm(vemb,nspin,x,comm_chempot,box_len,ngrid,norb,...
        tsmear,q_clu,q_env,vks_clu,vks_env,dm_tot);
    
    vks_clu_old = vks_clu;
    vks_env_old = vks_env;
    
    % update vks_clu and vks_env
    fprintf('update vks of cluster and env.\n');
    vks_clu =0.d0;
    vks_env =0.d0;
    % exchange
    [ex,vks_clu] = cal_ex(x, rho_clu);
    [ex,vks_env] = cal_ex(x, rho_env);
    % hartree
    [ex,vtmp] = cal_hartree(x, rho_clu); vks_clu = vks_clu + vtmp;
    [ex,vtmp] = cal_hartree(x, rho_env); vks_env = vks_env + vtmp;
    % external
    vks_clu = vks_clu + vext_clu;
    vks_env = vks_env + vext_env;
    
    % mixing
    vks_clu = (vks_clu_old+vks_clu)/2.0;
    vks_env = (vks_env_old+vks_env)/2.0;
end

% plot
figure
plot(x,ref_rho);
hold on;
plot(x,rho_clu);
plot(x,rho_env);
plot(x,rho_clu+rho_env,'ro');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstruct cluster KS system using penalty function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Reconstruct cluster with penalty technique?\n');
pause;



%
% localization with SCMD method
%
% ne = 7;
% [Q, R, piv] = qr(ev',0);
% Pc = ev*(ev(piv(1:ne),:)');
% S = Pc(piv(1:ne),:);
% Rchol = chol(S);
% Phi = Pc/(Rchol);

%%%%%%%%%% another coding of SCMD %%%%%%%%%
norb = 7;
occ_matrix = eye(norb,norb);
for q=1:norb
    occ_matrix(q,q) = sqrt(occ(q));
end
phi = ev(:,1:norb)*occ_matrix;
[Q,R] = qr( transpose(phi) );
loc_phi = phi*Q;

figure;
hold on;
for q=1:ne
    plot(x,ev(:,q));
end


figure
hold on
for q=1:ne
    plot(x,loc_phi(:,q).^2,'-')
end
plot(x,rho,'ro');
plot(x,sum(loc_phi(:,:).^2,2))


% make the density matrix from orbitals
%[dm] = wfr2dm(norb,ngrid,ev,occ);

