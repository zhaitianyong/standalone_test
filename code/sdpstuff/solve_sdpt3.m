% This function solves the SDP problem
%
% Babak Alipanahi
% University of Waterloo
% July 19, 2010
% edit: August 31, 2010
% edit: March 2, 2011 (SDPT3 customization)

function [X eigens slacks] = solve_sdpt3(U,equality_cons,upper_bounds,lower_bounds,vdw_bounds,f)

% finding equivalent U
% Find V such that V'*U'*e = 0
UTe = sum(U)';
[V ,~] = qr(UTe);
U = U*V(:,2:end);

% f_hb  = f[0];
% f_tau = f[1];
% f_tal = f[2];
% f_vdw = f[3];
% f_tas = f[4];

hydrogen_factor = f(1);
torsionu_factor = f(2);
torsionl_factor = f(5);
vdw_factor      = f(4);

num_upper_bounds    = size(upper_bounds,1);
num_equality_cons   = size(equality_cons,1);
num_lower_bounds    = size(lower_bounds,1);
num_vdw_bounds      = size(vdw_bounds,1);

[num_atoms sdp_dim] = size(U);

num_up_slacks  = 2*num_upper_bounds;
num_lo_slacks  = 2*num_lower_bounds;
num_vdw_slacks = 2*num_vdw_bounds;

gamma  = 50*(num_atoms)/num_upper_bounds;
lambda = 0*(num_atoms)/num_upper_bounds;
num_slacks = num_up_slacks + num_lo_slacks + num_vdw_slacks;

%num_cons = num_equality_cons + num_upper_bounds + num_lower_bounds + num_vdw_bounds + 1;
num_cons = num_equality_cons + num_upper_bounds + num_lower_bounds + num_vdw_bounds;


b = zeros(num_cons,1);
Vec_cons = nan(sdp_dim,num_cons);

S_slacks = nan(num_slacks,3);
C_slacks = nan(num_slacks,1);


cons_cnt = 1;
for i = 1:num_equality_cons
    ti = equality_cons(i,1);
    tj = equality_cons(i,2);
    tmpU = (U(ti,:) - U(tj,:))';
    Vec_cons(:,cons_cnt) = tmpU;
    b(cons_cnt) = equality_cons(i,3)^2;
    cons_cnt = cons_cnt + 1;
end

slack_cnt = 1;
for i = 1:num_upper_bounds
    ti = upper_bounds(i,1);
    tj = upper_bounds(i,2);
    tmpU = (U(ti,:) - U(tj,:))';
    Vec_cons(:,cons_cnt) = tmpU;
    b(cons_cnt) = upper_bounds(i,3)^2;
    S_slacks(slack_cnt,:) = [slack_cnt cons_cnt +1];
    C_slacks(slack_cnt)   = lambda;
    slack_cnt = slack_cnt + 1;
    S_slacks(slack_cnt,:) = [slack_cnt cons_cnt -1];
    if upper_bounds(i,4) >= 0
        C_slacks(slack_cnt)   = gamma;
    % hydrogen bonds
    elseif upper_bounds(i,4) == -1
        C_slacks(slack_cnt)   = hydrogen_factor*gamma;
    % torsion angles
    elseif upper_bounds(i,4) == -2
        C_slacks(slack_cnt)   = torsionl_factor*gamma;
    end
    slack_cnt = slack_cnt + 1;
    cons_cnt  = cons_cnt + 1;
end


for i = 1:num_lower_bounds
    ti = lower_bounds(i,1);
    tj = lower_bounds(i,2);
    tmpU = (U(ti,:) - U(tj,:))';
    Vec_cons(:,cons_cnt) = tmpU;
    b(cons_cnt) = lower_bounds(i,3)^2;
    S_slacks(slack_cnt,:) = [slack_cnt cons_cnt -1];
    C_slacks(slack_cnt)   = 0;
    slack_cnt = slack_cnt + 1;
    S_slacks(slack_cnt,:) = [slack_cnt cons_cnt +1];
    if lower_bounds(i,4) == -2
        C_slacks(slack_cnt) = torsionu_factor*gamma;
    end
    slack_cnt = slack_cnt + 1;
    cons_cnt  = cons_cnt + 1;    
end

for i = 1:num_vdw_bounds
    ti = vdw_bounds(i,1);
    tj = vdw_bounds(i,2);
    tmpU = (U(ti,:) - U(tj,:))';
    Vec_cons(:,cons_cnt) = tmpU;
    b(cons_cnt) = vdw_bounds(i,3)^2;
    S_slacks(slack_cnt,:) = [slack_cnt cons_cnt -1];
    C_slacks(slack_cnt)   = 0;
    slack_cnt = slack_cnt + 1;
    S_slacks(slack_cnt,:) = [slack_cnt cons_cnt +1];
    C_slacks(slack_cnt)   = vdw_factor*gamma;
    slack_cnt = slack_cnt + 1;
    cons_cnt  = cons_cnt + 1;    
end

% cent_vec = U'*ones(num_atoms,1);
% Vec_cons(:,cons_cnt) = cent_vec;
% b(cons_cnt) = 0;

num_indep_cons = length(b);

Cs = cell(2,1);
Cs{1} = C_slacks;
Cs{2} = -speye(sdp_dim);



blk = cell(2,3);
At  = cell(2,3);

blk{1,1} = 'l';
blk{1,2} = num_slacks;
blk{2,1} = 's';
blk{2,2} = sdp_dim;
blk{2,3} = ones(1,num_indep_cons);

At{1,1} = sparse(S_slacks(:,1),S_slacks(:,2),S_slacks(:,3),num_slacks,num_cons,num_slacks);
At{2,1} = [];
At{2,2} = Vec_cons;
At{2,3} = ones(num_indep_cons,1);

OPTIONS.printlevel = 3;

% [~,X] = sdpnalplus(blk,At,Cs,b,[],[],[],[],[],OPTIONS);
[~, X] =sqlp(blk,At,Cs,b,OPTIONS);

slacks = X{1};
Z = X{2};

[VZ LZ] = eigb(Z);
LZ = real(LZ);
LZ(LZ < 0) = 0;
X = LZ^0.5*(U*VZ)';
eigens = diag(LZ);

dLZ = diag(LZ);
ind = dLZ/max(dLZ); 
fprintf('\tRank: %d\n',sum(ind > 1e-6));

