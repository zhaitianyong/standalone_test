% This function builds a random protein
% based on the given sequence and phi and
% psi ranges
%
% Babak Alipanahi
% University of Waterloo
% May 10, 2010
% edit: June 17,2010

function [X, Comp, new_X, hComp] = ibuildprot(seq, num, phi, psi, A)

phi = pi*phi/180;
psi = pi*psi/180;

% add a dummy ALA at the end
seq = [seq 1];
phi(end+1) = -60;


seq_len = length(seq);
num_atoms = 0;
for i = 1:seq_len-1
    num_atoms = num_atoms + A(seq(i)).last_meta - A(seq(i)).first_meta + 1;
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Computing the atom omission map
kept_atom_cnt = 0;
atom_cnt = 0;
atoms_to_keep = true(1,num_atoms);
atoms_map = nan(num_atoms,2);
atoms_map(:,1) = 1:num_atoms;
hresidue_start = nan(1,seq_len-1);
for i = 1:seq_len-1
    hresidue_start(i) = kept_atom_cnt;
    cur_num_atoms = A(seq(i)).last_meta - A(seq(i)).first_meta + 1;
    atoms_to_keep(atom_cnt+1:atom_cnt+cur_num_atoms) = A(seq(i)).omission_map;
    atom_cnt = atom_cnt + cur_num_atoms;
    kept_atom_cnt = kept_atom_cnt + sum(A(seq(i)).omission_map);
end
index_kept = find(atoms_to_keep);
num_atoms_kept = length(index_kept);
atoms_map(index_kept,2) = 1:num_atoms_kept;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


X        = zeros(3,num_atoms);
bonds    = nan(4*num_atoms,2);
info     = nan(6,num_atoms);
atom_names = cell(num_atoms,1);
atom_types = zeros(num_atoms,1);
planes  = zeros(seq_len-1,6);
sides   = cell(seq_len-1,1);
schains = cell(seq_len-1,1);
new_schains = cell(seq_len-1,1);
residue     = nan(num_atoms,1);
residue_start = nan(1,seq_len);

% First residue only has valid PSI angle (no PHI)
res = seq(1);
curX = A(res).X;
max_j = size(A(res).dihedral,2);
for j = 1:max_j
    if strcmp(A(res).dihedral(j).name, 'PSI')
        psi_index = A(res).dihedral(j).atoms(1:3);
        break;
    end
end
psiX = curX(:, psi_index);

start = 0;
bond_cnt = 0;
for i = 1:seq_len-1
    res     = seq(i);
    nextres = seq(i+1);
    last    = A(res).last_meta;
    first   = A(res).first_meta;
    cur_num_atoms = last - first + 1;
    temp_X  = zeros(3,cur_num_atoms);
    residue_start(i) = start;
    %----------------------------------------------------
    %----------------------------------------------------
    % find the next N by psi
    next_N = ditocord(psiX,res,nextres,psi(i),'PSI');
    %----------------------------------------------------
    cur_N  = psiX(:,1);
    cur_CA = psiX(:,2);
    cur_C  = psiX(:,3);
    %----------------------------------------------------
    % Find the atoms in a plane
    % atoms in outX: CA O C N HN CA (res ~= PRO)
    %                CA O C N CD CA (res == PRO)
    outX = PlaneFinder([cur_CA cur_C next_N],res,nextres);
    cur_O = outX(:,2);
    %----------------------------------------------------
    % For all residues other than PRO
    if res ~= 15
        % atoms in plane > i: CA C O | i+1: N HN CA
        planes(i,:) = start + [3 (last-first)+(0:1) cur_num_atoms+(1:3)];
        if i == 1
            cur_HN = A(res).X(:,4);
        else
            cur_HN = next_HN;
        end
        cur_backbone = [cur_N cur_HN cur_CA cur_C cur_O];
        index_backbone = [1:3 (last - first) + (0:1)];
        index_sidechain = 4:(last-first)-1;
        atom_types(start + index_backbone)  = 1;    % backbone  = 1
        atom_types(start + index_sidechain) = 2;    % sidechain = 2
        temp_X(:,index_backbone) = cur_backbone;
        %----------------------------------------------------
        % atoms in sides{i} are N CA C and side chain atoms
        sides{i} = start + [1 3 (last-first) index_sidechain];
        num_cliq = size(A(res).cliqs,1);
        schains{i} = cell(num_cliq,1);
        for q = 1:num_cliq
            schains{i}{q} = A(res).cliqs{q} + residue_start(i);
        end
        %++++++++++++++++++++++++++++++++++++++++++++++++++++
        num_new_cliq = size(A(res).cliqs_woh,1);
        new_schains{i} = cell(num_new_cliq,1);
        for q = 1:num_new_cliq
           new_schains{i}{q} = A(res).cliqs_woh{q} + residue_start(i); 
        end
        %++++++++++++++++++++++++++++++++++++++++++++++++++++
        %----------------------------------------------------
        % atoms used for side chain alignement [cur_N cur_CA cur_C]
        cur_sidechain = SidechainFinder(cur_backbone(:,[1 3 4]),index_backbone([1 3 4]),A(res));
        cur_sidechain = cur_sidechain(:,index_sidechain);
        temp_X(:,index_sidechain) = cur_sidechain;
        % For PRO
    else
        % atoms in plane > i: CA C O | i+1: N CD CA
        planes(i,:) = start+[3 (last-first)+(0:1) cur_num_atoms+(1:3)];
        if i == 1
            cur_CD = A(res).X(:,4);
        else
            cur_CD = next_CD;
        end
        % no HN in PRO
        cur_backbone = [cur_N cur_CD cur_CA cur_C cur_O];
        index_backbone = [1:3 (last - first) + (0:1)];
        index_sidechain = 4:(last-first)-1;
        atom_types(start + index_backbone)  = 1;
        atom_types(start + index_sidechain) = 2;
        temp_X(:,index_backbone) = cur_backbone;
        %----------------------------------------------------
        % atoms in sides{i} are N CA C CD and side chain atoms
        sides{i} = start + [1 3 (last-first) 2 index_sidechain];
        num_cliq = size(A(res).cliqs,1);
        schains{i} = cell(num_cliq,1);
        for q = 1:num_cliq
            schains{i}{q} = A(res).cliqs{q} + residue_start(i);
        end
        %++++++++++++++++++++++++++++++++++++++++++++++++++++
        num_new_cliq = size(A(res).cliqs_woh,1);
        new_schains{i} = cell(num_new_cliq,1);
        for q = 1:num_new_cliq
           new_schains{i}{q} = A(res).cliqs_woh{q} + residue_start(i); 
        end
        %++++++++++++++++++++++++++++++++++++++++++++++++++++
        %----------------------------------------------------
        cur_sidechain = SidechainFinder(cur_backbone(:,[1 3 4]),index_backbone([1 3 4]),A(res));
        cur_sidechain = cur_sidechain(:,index_sidechain);
        temp_X(:,index_sidechain) = cur_sidechain;
    end
    % Put the computed residue into the protein
    X(:,start+1:start+cur_num_atoms) = temp_X;
    residue(start+1:start+cur_num_atoms) = num(i);
    for j = 1:cur_num_atoms
        atom_names{start+j} = A(res).atom(j+first-1).name;
        % pseudo atoms have type 3
        if regexp(atom_names{start+j}, 'PSEUD')
            atom_types(start+j) = 3;
        end
        info(:,start+j) = A(res).info(:,j+first-1);
        % info(6,j) is the pseudo atom corresponding to atom j
        if info(6,start+j)
            info(6,start+j) = info(6,start+j) + residue_start(i);
        end
    end
    
    bond_increment = size(A(res).meta_bonds,1);
    bonds(bond_cnt+(1:bond_increment),:) = start + A(res).meta_bonds;
    bond_cnt = bond_cnt + bond_increment;
    %----------------------------------------------------
    %----------------------------------------------------
    % this is the peptide bond
    % between C(i) and N(i+1)
    if i < seq_len - 1
        bond_cnt = bond_cnt + 1;
        if res ~= 15
            bonds(bond_cnt,:) = start + [index_backbone(4) cur_num_atoms+1];
        else
            bonds(bond_cnt,:) = start + [index_backbone(3) cur_num_atoms+1];
        end
    end
    %----------------------------------------------------
    next_CA = outX(:,6);
    next_HN = outX(:,5);
    next_CD = outX(:,5);
    % atoms required for determing next C' by PHI angle
    phiX = [cur_C next_N next_CA];
    %----------------------------------------------------
    % find the next C' by PHI
    next_C = ditocord(phiX,res,nextres,phi(i+1),'PHI');
    %----------------------------------------------------
    % This is for next residue's PSI angle
    psiX = [next_N next_CA next_C];
    %----------------------------------------------------
    start = start + cur_num_atoms;
end
bonds(isnan(bonds(:,1)),:) = [];
info(:,isnan(info(1,:))) = [];
residue_start(isnan(residue_start)) = [];

Comp.info = info;
Comp.bonds = bonds;

Comp.atom_names = atom_names;
Comp.atom_types = atom_types;

Comp.planes  = planes;
Comp.sides   = sides;
Comp.residue = residue;
Comp.schains = schains;
Comp.residue_bias = residue_start;

Comp.seq = seq;
Comp.num_seq = num;

residue_type = nan(length(residue),1);
for i = 1:length(residue_type)
   temp_res = residue(i);
   ind = num == temp_res;
   residue_type(i) = seq(ind);
end

Comp.residue_type = residue_type;

% We form the Cliques sparse matrix
num_atoms  = size(X,2);
num_planes = size(Comp.planes,1);
num_res    = size(Comp.schains,1);
num_sides  = 0;
for i = 1:num_res
    cliqs_per_side = size(Comp.schains{i},1);
    num_sides = num_sides + cliqs_per_side;
end
raw_estimate = num_planes*6 + num_sides*30;
I = nan(raw_estimate,1);
J = nan(raw_estimate,1);
V = ones(raw_estimate,1);

%---------------------------------------------------
%---------------------------------------------------
I(1:3) = 1:3;
J(1:3) = 1;
start = 3;
cliq_cnt = 1;
for i = 1:num_planes
    cliq_cnt = cliq_cnt + 1;
    num_atoms_clique = sum(Comp.planes(i,:) <= num_atoms);
    I(start+1:start+num_atoms_clique) = Comp.planes(i,Comp.planes(i,:) <= num_atoms);
    J(start+1:start+num_atoms_clique) = cliq_cnt;
    start = start + num_atoms_clique;
end

for i = 1:num_res
    cliqs_per_side = size(Comp.schains{i},1);
    for j = 1:cliqs_per_side
        cliq_cnt = cliq_cnt + 1;
        num_atoms_clique = sum(Comp.schains{i}{j} <= num_atoms);
        I(start+1:start+num_atoms_clique) = Comp.schains{i}{j}(Comp.schains{i}{j} <= num_atoms);
        J(start+1:start+num_atoms_clique) = cliq_cnt;
        start = start + num_atoms_clique;
    end
end

index_bad = isnan(I);
I(index_bad) = [];
J(index_bad) = [];
V(index_bad) = [];

Comp.Cq = sparse(I,J,V,num_atoms,cliq_cnt,2*sum(V));
%---------------------------------------------------
%---------------------------------------------------

% Dropping some of the Hydrogen atoms
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
new_X = X(:,atoms_to_keep);
hComp.info          = Comp.info(:,atoms_to_keep);
hComp.atom_names    = Comp.atom_names(atoms_to_keep);
hComp.atom_types    = Comp.atom_types(atoms_to_keep);
hComp.residue       = Comp.residue(atoms_to_keep);
hComp.seq           = Comp.seq;
hComp.num_seq       = Comp.num_seq;
hComp.residue_bias  = hresidue_start;


residue_type = nan(length(hComp.residue),1);
for i = 1:length(residue_type)
   temp_res = hComp.residue(i);
   ind = hComp.num_seq == temp_res;
   residue_type(i) = hComp.seq(ind);
end

hComp.residue_type = residue_type;

% new bonds
num_bonds = size(bonds,1);
new_bonds = nan(num_bonds,2);
for i = 1:num_bonds
    ti = bonds(i,1);
    tj = bonds(i,2);
    if ~isnan(atoms_map(ti,2)*atoms_map(tj,2))
        new_bonds(i,:) = [atoms_map(ti,2) atoms_map(tj,2)];
    end
end
new_bonds(isnan(new_bonds(:,1)),:) = [];
hComp.bonds = new_bonds;

% new planes
num_planes = size(planes,1);
new_planes = nan(num_planes,6);
for i = 1:num_planes
    for j = 1:6
        if planes(i,j) <= num_atoms
            new_planes(i,j) = atoms_map(planes(i,j),2);
            if isnan(new_planes(i,j))
                error('There is an error in the mapping atoms');
            end
        end
    end
end
hComp.planes = new_planes;

% new side chains
new_sides = sides;
for i = 1:num_res
    num_atoms_in_side = length(sides{i});
    index_bad = false(1,num_atoms_in_side);
    for j = 1:num_atoms_in_side
       if isnan(atoms_map(sides{i}(j),2))
          index_bad(j) = true; 
       end
    end
    new_sides{i}(index_bad) = [];
    new_sides{i} = atoms_map(new_sides{i},2)';
end
hComp.sides = new_sides;


% new side chain cliques
for i = 1:num_res
    num_cliq = size(new_schains{i},1); 
    for q = 1:num_cliq
        cliq_len = length(new_schains{i}{q});
        for j = 1:cliq_len
            new_schains{i}{q}(j) = atoms_map(new_schains{i}{q}(j),2);
            if isnan(new_schains{i}{q}(j))
                 error('There is an error in the mapping atoms');
            end
        end
    end
end
hComp.schains = new_schains;

% Now, form the cliq matrix

num_atoms  = size(new_X,2);
num_planes = size(hComp.planes,1);
num_res    = size(hComp.schains,1);
num_sides  = 0;
for i = 1:num_res
    cliqs_per_side = size(hComp.schains{i},1);
    num_sides = num_sides + cliqs_per_side;
end
raw_estimate = num_planes*6 + num_sides*30;
I = nan(raw_estimate,1);
J = nan(raw_estimate,1);
V = ones(raw_estimate,1);

%---------------------------------------------------
%---------------------------------------------------
I(1:3) = 1:3;
J(1:3) = 1;
start = 3;
cliq_cnt = 1;
for i = 1:num_planes
    cliq_cnt = cliq_cnt + 1;
    num_atoms_clique = sum(hComp.planes(i,:) <= num_atoms);
    I(start+1:start+num_atoms_clique) = hComp.planes(i,hComp.planes(i,:) <= num_atoms);
    J(start+1:start+num_atoms_clique) = cliq_cnt;
    start = start + num_atoms_clique;
end

for i = 1:num_res
    cliqs_per_side = size(hComp.schains{i},1);
    for j = 1:cliqs_per_side
        cliq_cnt = cliq_cnt + 1;
        num_atoms_clique = sum(hComp.schains{i}{j} <= num_atoms);
        I(start+1:start+num_atoms_clique) = hComp.schains{i}{j}(hComp.schains{i}{j} <= num_atoms);
        J(start+1:start+num_atoms_clique) = cliq_cnt;
        start = start + num_atoms_clique;
    end
end

index_bad = isnan(I);
I(index_bad) = [];
J(index_bad) = [];
V(index_bad) = [];

hComp.Cq = sparse(I,J,V,num_atoms,cliq_cnt,2*sum(V));
hComp.atoms_map = atoms_map;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


