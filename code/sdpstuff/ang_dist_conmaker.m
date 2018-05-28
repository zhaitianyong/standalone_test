% This function translates dihedral angle constraints into cooresponding
% distance constraints: d(Ci-1,Ci) for phi and d(Ni,Ni+1) for psi
% 
%
% Babak Alipanahi
% University of Waterloo
% March 8, 2011

function [lo_cons up_cons] = ang_dist_conmaker(phi_cons,psi_cons,Comp)

atom_names = Comp.atom_names;
residue = Comp.residue;
seq_num = Comp.num_seq;
bias    = Comp.residue_bias;

num_phi = size(phi_cons,1);
num_psi = size(psi_cons,1);
num_all = num_phi + num_psi;

lo_cons = nan(num_all,4);
up_cons = nan(num_all,4);

for i = 1:num_phi
    if phi_cons(i,2) - phi_cons(i,1) < 360 && i > 1
        dmin = -inf;
        d1 = DCC(phi_cons(i,1));
        d2 = DCC(phi_cons(i,2));
        if phi_cons(i,2) * phi_cons(i,1) > 0
            % both on the same side, lower and
            % upper bounds are required
            dmax = max(d1,d2);
            dmin = min(d1,d2);
        else
            % now lower bound is necessary
            dmax = max(d1,d2);
        end
        atoms_cur = atom_names(residue == seq_num(i));
        atoms_per = atom_names(residue == seq_num(i-1));        
        cur_index = index_finder(atoms_cur,'C') + bias(i);
        per_index = index_finder(atoms_per,'C') + bias(i-1);        
        up_cons(i,:) = [per_index cur_index dmax -2];        
        if ~isinf(dmin)
            lo_cons(i,:) = [per_index cur_index dmin -2];
        end          
   end
end

for i = 1:num_psi
    if psi_cons(i,2) - psi_cons(i,1) < 360 && i < num_psi
        dmin = -inf;
        d1 = DNN(psi_cons(i,1));
        d2 = DNN(psi_cons(i,2));
        if psi_cons(i,2) * psi_cons(i,1) > 0
            % both on the same side, lower and
            % upper bounds are required
            dmax = max(d1,d2);
            dmin = min(d1,d2);
        else
            % now lower bound is necessary
            dmax = max(d1,d2);
        end
        atoms_cur = atom_names(residue == seq_num(i));
        atoms_nex = atom_names(residue == seq_num(i+1));        
        cur_index = index_finder(atoms_cur,'N') + bias(i);
        nex_index = index_finder(atoms_nex,'N') + bias(i+1);        
        up_cons(num_phi + i,:) = [cur_index nex_index dmax -2];        
        if ~isinf(dmin)
            lo_cons(num_phi + i,:) = [cur_index nex_index dmin -2];
        end          
   end
end
lo_cons(isnan(lo_cons(:,1)),:) = [];
up_cons(isnan(up_cons(:,1)),:) = [];

%**********************************************************************************

%**********************************************************************************
function out = DCC(phi)

phi = phi * pi/180;
out = 3.221 - 0.4886*cos(1.044*phi);
%**********************************************************************************
function out = DNN(psi)

psi = psi * pi/180;
out = 3.157 - 0.5022*cos(1.046*psi);

%**********************************************************************************
function ind = index_finder(list,query)

ind = -1;
for i = 1:numel(list)
    if strcmp(list{i},query)
        ind = i;
        break
    end
end
%**********************************************************************************



