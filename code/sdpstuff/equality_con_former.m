% This function extracts the equality constraints
%
% Babak Alipanahi
% University of Waterloo
% July 16, 2010
% revise: August 19, 2010
% 
%     type = 1 > cons for preserving cliques
%     type = 2 > all equality cons, for post-processing    

function equality_cons = equality_con_former(X,Comp,type)

if ~exist('type','var')
    type = 1;
end

Cq = Comp.Cq;
[num_atoms num_cliq] = size(Cq);

G = zeros(num_atoms);
% add a dummy constraint
% between HA and H of the
% first residue
G(2,4) = 1;


for i = 1:num_cliq
    index_cliq = find(Cq(:,i));
    dim = Comp.cliq_dims(i);
    %base = index_cliq(1:dim+1);
    %if pcarank(X(:,base)) < dim
    base = BaseFormer(X(:,index_cliq),dim);
    base = index_cliq(base);
    %end
    for j = 1:dim+1
        for k = j+1:dim+1
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
        dim = Comp.cliq_dims(i);
        base = BaseFormer(X(:,index_cliq),dim);
        base = index_cliq(base);        
        for j = 1:dim+1
            for k = j+1:dim+1
                tj = base(j);
                tk = base(k);
                G(tj,tk) = 1;
            end
        end
        index_nonbase = setdiff(index_cliq,base);
        if ~isempty(index_nonbase)
           for j = 1:dim+1
              for k = 1:numel(index_nonbase)
                 tj = base(j);
                 tk = index_nonbase(k);
                 G(tj,tk) = 1;
              end
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


function out = BaseFormer(X,d)

num = size(X,2);
if num < d + 1
    out = 1:num;
else
    rankX = pcarank(X);
    rand_index = randperm(num);
    test_points = X(:,rand_index(1:rankX + 1));
    while pcarank(test_points) < rankX
        rand_index = randperm(num);
        test_points = X(:,rand_index(1:rankX + 1));
    end
    out = rand_index(1:rankX + 1);
end


