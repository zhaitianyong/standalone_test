function equality_cons = equality_cons_former(X,Comp,type)

if ~exist('type','var')
    type = 1;
end

Cq = Comp.Cq;
[num_atoms,num_cliq] = size(Cq);

G = zeros(num_atoms);
% add a dummy constraint
% between HA and H of the
% first residue
G(2,4) = 1;


for i = 1:num_cliq
    index_cliq = find(Cq(:,i));
    num = size(index_cliq,1);
    base = index_cliq;
    %end
    for j = 1:num
        for k = j+1:num
            tj = base(j);
            tk = base(k);
            G(tj,tk) = 1;
        end
    end
    if type > 1
        num_points_cliq = numel(index_cliq);
        for j = 1:num_points_cliq
            for k = j+1:num_points_cliq
                tj = index_cliq(j);
                tk = index_cliq(k);
                G(tj,tk) = 1;
            end
        end
    end
    if type < 0
        index_cliq = find(Cq(:,i));
        num = Comp.cliq_dims(i);
        base = BaseFormer(X(:,index_cliq),num);
        base = index_cliq(base);        
        for j = 1:num+1
            for k = j+1:num+1
                tj = base(j);
                tk = base(k);
                G(tj,tk) = 1;
            end
        end
    end
end

%G = max(G,G');
%G = triu(G);
[I J] = find(G);
    

D = distmex(X);
Vals = (D((J-1)*num_atoms+I)).^0.5;

equality_cons(:,1) = I;
equality_cons(:,2) = J;
equality_cons(:,3) = Vals;