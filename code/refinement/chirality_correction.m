% This function checks correct chirality
%
% Babak Alipanahi
% University of Waterloo
% February 22, 2011

function outX = chirality_correction(X,Comp,chiral_err)

seq = Comp.seq;
num_seq = Comp.num_seq;
atom_names = Comp.atom_names;
atom_types = Comp.atom_types;
residue_bias = Comp.residue_bias;

num_res = length(seq) - 1;
outX = X;

for i = 1:num_res
    res = seq(i);
    if ~chiral_err(i) && res ~= 8
        chiral_atoms = {'CA','N','C','CB','HA'};
        chiral_atoms_index = nan(1,5);        
        temp_atoms = Comp.residue == num_seq(i);
        temp_X = X(:,temp_atoms);
        temp_atom_names = atom_names(temp_atoms);
        temp_atom_types = atom_types(temp_atoms);
        for j = 1:5
            for k = 1:length(temp_atom_names)
                if strcmp(chiral_atoms{j},temp_atom_names{k})
                    chiral_atoms_index(j) = k;
                    break;
                end
            end
        end
        cX = temp_X(:,chiral_atoms_index);
        
        
        v1 = cX(:,2) - cX(:,1);     % CA-N vector
        v2 = cX(:,3) - cX(:,1);     % CA-C vector
        n  = cross(v1,v2);          % CA-N-C plane normal vector    
        n  = n/norm(n);             % normalized
        
        side_chain = find(temp_atom_types > 1);
                
        sChain = temp_X(:,side_chain);
        tSChain = nan(3,length(side_chain));
        for j = 1:length(side_chain)
            tSChain(:,j) = plane_reflector(cX(:,1),n,sChain(:,j));
        end
        outX(:,side_chain + residue_bias(i)) = tSChain;
    end
end


% This function computes the reflection of P w.r.t plane
% identified by the normal vector n and a point O
function out = plane_reflector(O,n,P)

v  = P - O;
vn = (v'*n)*n;
vp = v - vn;
vr = vp - vn;
out = O + vr;




