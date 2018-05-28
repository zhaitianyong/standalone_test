% This function writes the output into a
% PDB file
%
% Babak Alipanahi
% University of Waterloo
% July 26, 2010


%function lo_bounds = vdw_bound_maker(Comp,up_bounds)
function lo_bounds = vdw_bound_maker(Comp)


permitted_peneteration = 0.9;

num_atoms = size(Comp.info,2);
G = ones(num_atoms);

Cq = Comp.Cq;
num_cliques = size(Cq,2);

pseudo_atoms = Comp.info(5,:) == 0;
G(pseudo_atoms,:) = 0;
G(:,pseudo_atoms) = 0;

for q = 1:num_cliques
   index_cliq = find(Cq(:,q));
   cliq_size  = numel(index_cliq);
   for i = 1:cliq_size
       for j = i+1:cliq_size
            G(index_cliq(i),index_cliq(j)) = 0;
       end
   end
end

% num_bonds = size(Comp.bonds,1);
% for i = 1:num_bonds
%    G(Comp.bonds(i,1),Comp.bonds(i,2)) = 0; 
% end

G = min(G,G');
G = triu(G,1);

[I J] = find(G);
max_i = length(I);
lo_bounds = zeros(max_i,4);
for i = 1:max_i
    lo_bounds(i,1:3) = [I(i) J(i) Comp.info(2,I(i)) + Comp.info(2,J(i))];
end
lo_bounds(:,3) = lo_bounds(:,3)*permitted_peneteration;
lo_bounds(:,4) = 1;

% for i = 1:size(up_bounds,3)
%    s = up_bounds(i,1);
%    t = up_bounds(i,2);
%    d = up_bounds(i,3);
%    ind = (lo_bounds(:,1) == s & lo_bounds(:,2) == t) | (lo_bounds(:,2) == s & lo_bounds(:,1) == t);
%    if sum(ind) > 1
%        error('there should be one such atom pair')
%    end
%    if lo_bounds(ind,3) >= d
%       lo_bounds(ind,3) = -1; 
%    end
% end
% lo_bounds(lo_bounds(:,3)==-1,:) = [];




% function lo_bounds = vdw_bound_maker(Comp)
% 
% permitted_peneteration = 0.9;
% 
% num_atoms = size(Comp.info,2);
% G = ones(num_atoms);
% 
% num_planes = size(Comp.planes,1);
% num_sides  = size(Comp.sides,1);
% 
% pseudo_atoms = Comp.info(5,:) == 0;
% G(pseudo_atoms,:) = 0;
% G(:,pseudo_atoms) = 0;
% 
% for i = 1:num_planes
%     tatoms = Comp.planes(i,:);
%     tatoms(tatoms > num_atoms) = [];
%     tatoms(isnan(tatoms)) = [];
%     max_j = length(tatoms);
%     for j = 1:max_j
%        for k = j+1:max_j
%           G(tatoms(j),tatoms(k)) = 0; 
%        end
%     end
% end
% 
% for i = 1:num_sides
%     tatoms = Comp.sides{i};
%     max_j = length(tatoms);
%     for j = 1:max_j
%        for k = j+1:max_j
%            if tatoms(j) == 4 && tatoms(k) == 6
%                keyboard
%            end
%           G(tatoms(j),tatoms(k)) = 0; 
%        end
%     end
% end
% 
% num_bonds = size(Comp.bonds,1);
% for i = 1:num_bonds
%    G(Comp.bonds(i,1),Comp.bonds(i,2)) = 0; 
% end
% 
% G = min(G,G');
% %G = G - diag(diag(G));
% G = triu(G,1);
% 
% [I J] = find(G);
% max_i = length(I);
% lo_bounds = zeros(max_i,3);
% for i = 1:max_i
%     lo_bounds(i,:) = [I(i) J(i) Comp.info(2,I(i)) + Comp.info(2,J(i))];
% end
% lo_bounds(:,3) = lo_bounds(:,3)*permitted_peneteration;


