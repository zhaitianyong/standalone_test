% This function adds pseudo atoms
%
% Babak Alipanahi
% University of Waterloo
% Aug 15, 2010

function   [Comp X xplor_flag] = pseudo_adder(Comp_raw,X_raw,A)


amino_3L = {'ALA', 'ARG', 'ASN','ASP','CYS', 'GLU', 'GLN','GLY', 'HIS', 'ILE',...
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL','CYSS'};

num_models = numel(Comp_raw);
Comp = cell(num_models,1);
X    = cell(num_models,1);

% checking xplor atom nomenoclature
xplor_flag = 0;
for i = 1:size(X_raw{1},2)
    if ~isempty(regexp(Comp_raw{1}.atom_names{i},'HB1', 'once')) ||  ~isempty(regexp(Comp_raw{1}.atom_names{i},'HA1', 'once'))
        tres = Comp_raw{1}.residue(i);
        tres = Comp_raw{1}.seq(find(Comp_raw{1}.num_seq == tres,1));
        if tres == 1
            continue
        end
        % since ALA has HB1 we need to exlcude it
        if tres ~= 1
            xplor_flag = 1;
        end
        if tres == 8
            xplor_flag = 1;
        end
        break
    end
end

for m = 1:num_models
   atom_cnt = 1;
   tComp = Comp_raw{m};
   tX    = X_raw{m};
   tseq  = tComp.seq;
   tnum_seq = tComp.num_seq;
   tresidue = tComp.residue;
   num_res  = numel(tseq);
   num_atom = numel(tComp.atom_names);
   % each resisdue has at most 5 pseudo atoms
   raw_size = num_atom + num_res*5;
   X{m} = nan(3,raw_size);
   if xplor_flag
       for i = 1:num_atom
           temp_name = tComp.atom_names{i};
           tres = tComp.residue(i);
           tres = find(tnum_seq == tres,1,'first');
           tres_3L = amino_3L{tseq(tres)};
           new_name = atom_nom(tres_3L,temp_name,'xplor');
           if ~isempty(new_name)
               tComp.atom_names{i} = new_name;
           end
       end
   end
   Comp{m}.atom_names = cell(raw_size,1);
   Comp{m}.atom_types = nan(raw_size,1);
   Comp{m}.residue    = nan(raw_size,1);
   for i = 1:num_res              
       index_res = tresidue == tnum_seq(i);
       rX = tX(:,index_res);
       res_atoms = size(rX,2);
       % first we put in the known atoms
       res_atom_range = atom_cnt + (1:res_atoms) - 1;
       X{m}(:,res_atom_range) = rX;
       Comp{m}.atom_names(res_atom_range) = tComp.atom_names(index_res);
       Comp{m}.atom_types(res_atom_range) = tComp.atom_types(index_res);
       Comp{m}.residue(res_atom_range) = tComp.residue(index_res);
       atom_cnt = atom_cnt + res_atoms;
       % now time too put the pseudo atoms in
       tres = tseq(i);
       tA   = A(tres); 
       tres_atom_names = tComp.atom_names(index_res);       
       tres_orig_names = {tA.atom.name};
       tres_orig_names = tres_orig_names(tA.first_meta:tA.last_meta);
       num_pseudos = numel(tA.pseudos); 
       for j = 1:num_pseudos
           pX = zeros(3,1);
           num_atoms_in_pseudo = numel(tA.pseudos{j});
           atom_ind = index_finder(tres_atom_names,tres_orig_names(tA.pseudos{j}));
           if any(isnan(atom_ind))
               keyboard
           end
           for k = 1:num_atoms_in_pseudo
               pX = pX + rX(:,atom_ind(k));
           end
           pX = pX/num_atoms_in_pseudo;
           X{m}(:,atom_cnt) = pX;
           Comp{m}.atom_names{atom_cnt} = tA.pseudo_names{j};
           Comp{m}.atom_types(atom_cnt) = -1;
           Comp{m}.residue(atom_cnt) = tnum_seq(i);
           atom_cnt = atom_cnt + 1;
       end     
   end
   index_bad = isnan(X{m}(1,:));
   X{m}(:,index_bad) = [];
   Comp{m}.atom_names(index_bad) = [];
   Comp{m}.atom_types(index_bad) = [];
   Comp{m}.residue(index_bad)    = [];
   
end


function ind = index_finder(atoms,names)

if iscell(names)
    ind = nan(1,numel(names));
    for j = 1:numel(names)
        for i = 1:numel(atoms)
            if strcmp(atoms{i},names{j})
                ind(j) = i;
                break
            end
        end
    end
else
    ind = nan(1);    
    for i = 1:numel(atoms)
        if strcmp(atoms{i},names)
            ind = i;
            break
        end
    end    
end

