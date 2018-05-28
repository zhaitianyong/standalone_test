function PRASD(protein_pathway,protein_name,in_min_res,in_max_res,output_pathway)
% mac; linux

if isdeployed
  %fid = fopen(fullfile(ctfroot,'myfolder','myfile.dat'))
  mydir = ctfroot ;
%end
addpath([mydir '\code']);
addpath([mydir '\hanso']);
addpath([mydir '\tt']);
addpath([mydir '\matrix_completion']);
in_min_res=str2num(in_min_res);
in_max_res=str2num(in_max_res);
end

load 'ainfo.mat'

% Number of BFGS iterations
gd_tol  = 10^-9;
cg_iter = 500;
f = [10 10 10 10 10]; % [f_hb, f_tau, f_tal, f_vdw, f_tas]



fprintf('==========================================================================\n')


% protein_name_cell={'1b4r','1cn7','1g6j','2gjy','2k49','2k62','2k7h','2kte','2yt0','2l7b','2e8o','2l3o'};

%protein_name='1g6j';

tstart=tic;
   

in_hbond_file = '';
% in_max_res = [];
% in_min_res = [];

% switch protein_name
%     case '1g6j'
%         in_seq_file   = '1g6j.seq';
%         in_upl_file   = {'1g6j.upl'};
%         in_ang_file   = '1g6j.aco';
%     case '1cn7'
%         in_seq_file = '1cn7.seq';
%         in_upl_file = {'1cn7.upl'};
%         in_ang_file = '1cn7.aco';
%         in_min_res = 2;
%         in_max_res = 105;        
%     case '2k7h'
%         in_seq_file = '2k7h.seq';
%         in_upl_file = {'2k7h2.upl','2k7h.upl'};
%         in_ang_file = '2k7h.aco';
%     case '2kte'
%         in_seq_file = '2kte.seq';
%         in_upl_file = {'2kte.upl'};
%         in_ang_file = '2kte.aco';
%     case '1b4r'
%         in_seq_file = '1b4r.seq';
%         in_upl_file = {'1b4r.upl'};
%         in_ang_file = '1b4r.aco';
%         in_min_res = 8;
%         in_max_res = 87;        
%     case '2k49'
%         in_seq_file = '2k49.seq';
%         in_upl_file = {'2k49.upl'};
%         in_ang_file = '2k49.aco';
%     case '2k62'
%         in_seq_file = '2k62.seq';
%         in_upl_file = {'2k62.upl'};
%         in_ang_file = '2k62.aco';
%     case '2gjy'
%         in_seq_file = '2gjy.seq';
%         in_upl_file = {'2gjy.upl'};
%         in_ang_file = '2gjy.aco';        
%     case '2yt0'
%         in_seq_file = '2yt0.seq';
%         in_upl_file = {'2yt0.upl'};
%         in_ang_file = '2yt0.aco';
%     case '2l7b'
%         in_seq_file = '2l7b.seq';
%         in_upl_file = {'2l7b.upl'};
%         in_ang_file = '2l7b.aco';
%     case '2e8o'
%         in_seq_file = '2e8o.seq';
%         in_upl_file = {'2e8o.upl'};
%         in_ang_file = '2e8o.aco';        
%     case '2l3o'
%         in_seq_file = '2l3o.seq';
%         in_upl_file = {'2l3o.upl'};
%         in_ang_file = '2l3o.aco';
%         in_min_res = 33;
%         in_max_res = 156;
%     case '2mz7'
%         in_seq_file = '2mz7.seq';
%         in_upl_file = {'2mz7.upl'};
%         in_ang_file = '2mz7.aco';
%         in_min_res = 267;
%         in_max_res =312;
% end
in_seq_file   = [protein_name '.seq'];
in_upl_file   = {[protein_name '.upl']};
in_ang_file   = [protein_name '.aco'];

protein_path = [protein_pathway '\' protein_name '\'];
fprintf('*************************************************************************\n');
fprintf('Protein: %s\n', protein_name);
fprintf('-Reading input files...\n')
% Reading input data
%==========================================================================
seq_file = [protein_path in_seq_file];
[seq, num] = seq_reader(seq_file);
max_res = max(num);
min_res = min(num);
    
num_upl = size(in_upl_file, 2);
upl_file = cell(1, num_upl);
for i = 1:num_upl
    upl_file{i} = [protein_path in_upl_file{i}];
end
if ~isempty(in_hbond_file)
    hbond_file = [protein_path in_hbond_file];
    hbond_write_file = [protein_path protein_name '_hbo.upl'];
    hbond_reader(hbond_file,hbond_write_file);
    num_upl = num_upl + 1;
    upl_file{num_upl} = hbond_write_file;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
raw_up    = cell(1, num_upl);
raw_up_ho = cell(1, num_upl);
temp_min_res = +inf;
temp_max_res = -inf;
for i = 1:num_upl
    raw_up{i} = dist_reader(upl_file{i}, A);
    temp_min_res = min(temp_min_res, min([raw_up{i}.tres raw_up{i}.sres]));
    temp_max_res = max(temp_max_res, max([raw_up{i}.tres raw_up{i}.sres]));
end
if isempty(in_min_res)
    min_res = max(temp_min_res-1, min_res);
else
    min_res = in_min_res;
end
if isempty(in_max_res)
    max_res = min(temp_max_res+1, max_res);
else
    max_res = in_max_res;
end
    
% Remove informationless (w/o any constraints)
% parts from and N- and C-terminus
ind_del_N = num < min_res;
ind_del_C = num > max_res;
ind_del   = ind_del_N | ind_del_C;
seq(ind_del) = [];
num(ind_del) = [];
    
% dihedral angle constraints
ang_file = [protein_path in_ang_file];
[phi_cons, psi_cons] = ang_reader(ang_file, num);
%==========================================================================
trend = toc(tstart);
fprintf('\tdone: %4.1f sec\n', trend)
fprintf('-Sampling a random molecule...\n')
% Generating a random structure
%==========================================================================
[phi, psi] = ang_sampler(seq,phi_cons,psi_cons);

 [rand_X, Comp] = ibuildprot(seq,num,phi,psi,A);
 
if exist(ang_file, 'file')
    [ang_lo_cons, ang_up_cons] = ang_dist_conmaker(phi_cons,psi_cons,Comp);
end
%==========================================================================
trand = toc(tstart);
fprintf('\tdone: %4.1f sec\n',trand - trend)
fprintf('-Forming constraints...\n')
% Generating upper and lower bounds constraints
%==========================================================================
start = 0;
up_bounds = nan(50000,4);
for i = 1:num_upl
    temp_upl = upper_maker(raw_up{i}, Comp);
    up_bounds(start+1:start+size(temp_upl,1),:) = temp_upl;
    start = start + size(temp_upl,1);
end
up_bounds(isnan(up_bounds(:,1)), :) = [];
% adding torsion-angle constraints
if exist(ang_file, 'file')
   up_bounds = [up_bounds; ang_up_cons];
end
    
 % equality cons
 eq_cons = equality_cons_former(rand_X,Comp);
 % vdw bounds
 vdw_bounds = vdw_bound_maker(Comp);
    
 lo_bounds     = vdw_bounds;
 if exist(ang_file, 'file')
    lo_bounds    = [lo_bounds; ang_lo_cons];
end
%=========================================================================
% %
num_atom=size(Comp.atom_names,1);
cons_column=[eq_cons;up_bounds(:,1:3)];
ori_d=sparse(cons_column(:,1),cons_column(:,2),cons_column(:,3),num_atom,num_atom);
ori_d=full(ori_d);
d=ori_d+ori_d';
tbounds = toc(tstart);
fprintf('\tdone: %4.1f sec\n',tbounds - trand)
fprintf('Forming Triangle inequality constraint...\n')
fprintf('\n');
Dtri= tri3stri(d);
ttriang = toc(tstart);
fprintf('\tdone: %4.1f sec\n',ttriang - tbounds)
disp('*********************************************************************')
sample_rate = 10*5/num_atom;
[Dret,Dretsample ] = triChange( Dtri, d,num_atom );
Drecover = randmsamplesym(Dretsample,sample_rate);
DrecoverF = Drecover+d;
DrecoverF2=DrecoverF.^2;
% %
% scaledASD
ttriang_then = toc(tstart);
fprintf('solving the MC problem by ASD...\n')
fprintf('\n');
opts = default_opts();
[m,n]=size(DrecoverF2);
r=5;
[II,JJ]=find(DrecoverF2);
Omega2=sub2ind([m,n],II,JJ);
data2=DrecoverF2(Omega2);
% start2 = make_start_x('ASD',m,n,r,Omega2,data2);
% [Mout2, Out2] = ASD(m,n,r,Omega2,data2,start2,opts);
start2 = make_start_x('ScaledASD',m,n,r,Omega2,data2);
[Mout2, Out2] =ScaledASD(m,n,r,Omega2,data2,start2,opts);
tcompu_end = toc(tstart);
fprintf('\tdone: %4.1f sec\n',tcompu_end - ttriang_then)
disp('*********************************************************************')
%==========================================================================
% Post-processing
%==========================================================================
% ASD_post_processing
fprintf('Post processing...\n')
fprintf('\n');
H=eye(num_atom)-ones(num_atom,num_atom)/num_atom;
G2=-1/2*H*Mout2*H;
[VZ2,LZ2]=eigb(G2);
LZ2=real(LZ2);
LZ2(LZ2<0)=0;
Z2=LZ2^0.5*VZ2';
rawX=Z2(1:3,:); 
disp('*********************************************************************')
% Analysis of the output
%==========================================================================
fprintf('Violations (raw)\n')
check_eq_cons =eq_cons;
report = protchecker(rawX(1:3,:),Comp,check_eq_cons,lo_bounds,up_bounds,1);
disp('*********************************************************************')
%==========================================================================
% Post-processing
%==========================================================================
fprintf('Post-Processing...\n')
fprintf('\n');
f = [10 10 10 10 10];
W = [2 1 1 -1];

% Phase I - GD
% Refinement by HANSO
[pX, hinfo] = hanso_post_processing(rand_X, rawX, Comp, lo_bounds, up_bounds, W, f);
fprintf('Violations (GD-I)\n')
p_report = protchecker(pX,Comp,check_eq_cons,lo_bounds,up_bounds,1);
%---------------------------------------------------------------------------
% simple fix using the fact that most residues lie on the left half of
% Ramachandran plot
if sum(p_report.phi(~isnan(p_report.phi)) > 0) > 0.5*length(p_report.phi)
    pX(1,:) = -pX(1,:);
end
p_report = protchecker(pX,Comp,check_eq_cons,lo_bounds,up_bounds,0);
fprintf('\nCorrecting chiralities:\n')
pXc = chirality_correction(pX,Comp,p_report.chiral);
fprintf('Violations (after fixing chiralities)\n')
p_report = protchecker(pXc,Comp,check_eq_cons,lo_bounds,up_bounds,1);
disp('*********************************************************************')

[pXc, c_info] = hanso_post_processing(rand_X,pXc,Comp,lo_bounds,up_bounds,W,f);
fprintf('Violations (GD-II)\n')
pc_report = protchecker(pXc,Comp,check_eq_cons,lo_bounds,up_bounds,1);
fprintf('\nCorrecting chiralities:\n')
pXc = chirality_correction(pXc,Comp,pc_report.chiral);
disp('*********************************************************************')
[fX, wh_info] = hanso_post_processing(rand_X,pXc,Comp,lo_bounds,up_bounds,W,f);
 fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
 fprintf('\nViolations (FINAL)\n')
 fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
 final_report = protchecker(fX,Comp,check_eq_cons,lo_bounds,up_bounds,1);
 disp('*********************************************************************')
%==========================================================================


tpost_end = toc(tstart);
fprintf('\tdone: %4.1f sec\n',tpost_end-tcompu_end)
disp('*********************************************************************')
fprintf('\n Writing pdb\n')
pdb_writer_file=[output_pathway '\' protein_name '.pdb'];
pdb_writer(pdb_writer_file,fX,Comp)
twrite_end = toc(tstart);
fprintf('\tdone: %4.1f sec\n',twrite_end-tpost_end)
fprintf('==========================================================================\n')
fprintf('\tOverall time: %4.1f sec\n',toc(tstart))
fprintf('done!\n')
fprintf('==========================================================================\n')
