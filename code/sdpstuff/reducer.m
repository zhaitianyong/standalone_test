% This function reduces the cliques and computes
% U in K = UZU'
%
% Babak Alipanahi
% University of Waterloo
% August 18, 2010

function [U, cliq_dims] = reducer(X, Comp, type)

if ~exist('type', 'var')
    type = 0;
end

Cq  = Comp.Cq;
[num_atoms num_cliq] = size(Cq);

if type > 0
    U = [];
    cliq_dims = zeros(num_cliq,1);    
    for i = 1:num_cliq
        cliq_index = logical(Cq(:,i));
        qX = X(:,cliq_index);
        dim = pcarank(qX);
        cliq_dims(i) = dim;
    end
else    
    % Compute U for each clique
    cliq_dims = zeros(num_cliq,1);
    qU = cell(num_cliq,1);
    for i = 1:num_cliq
        cliq_index = logical(Cq(:,i));
        qX = X(:,cliq_index);
        qD = distmex(qX);
        qK = kdagger(qD);
        len_cliq = sum(cliq_index);
        [tU, ~] = eigb(qK);
        dim = pcarank(qX);
        cliq_dims(i) = dim;
        qU{i} = tU(:,1:dim);
        qU{i} = [qU{i} ones(len_cliq,1)/sqrt(len_cliq)];
    end
    
    % intersection of cliques
    % 1-intersection of cliques row-wise
    iSq = Cq'*Cq;
    rowU = cell(num_cliq,1);
    for i = 1:num_cliq
        num_intersect = find(iSq(i,:));
        num_intersect(num_intersect < i) = [];
        len_intersect = length(num_intersect);
        if len_intersect < 1
            error('At least one clique in the intersection')
        end
        rowU{i} = qU{i};
        for j = 2:len_intersect
            tj = num_intersect(j);
            [rowU{i} Cq(:,i)] = intersecter(rowU{i},qU{tj},Cq(:,i),Cq(:,tj));
        end
    end
    % 2-intersection of row-wise intersections
    % merged = Cq(:,1);
    % U = rowU{1};
    % for i = 2:num_cliq
    %     [U merged] = intersecter(U,rowU{i},merged,Cq(:,i));
    % end
    
    num_cliq_level = num_cliq;
    old_num_cliq_level = num_cliq_level;
    num_levels = floor(log2(num_cliq)) + 1;
    
    oldU = rowU;
    oldCq = Cq;
    for i = 1:num_levels
        num_cliq_level = ceil(num_cliq_level/2);
        tempU = cell(num_cliq_level,1);
        tCq = zeros(num_atoms,num_cliq_level);
        for j = 1:num_cliq_level
            qU1 = oldU{(j-1)*2+1};
            q1  = oldCq(:,(j-1)*2+1);
            tempU{j} = qU1;
            tCq(:,j) = q1;
            if (j-1)*2 + 2 <= old_num_cliq_level
                qU2 = oldU{(j-1)*2+2};
                q2  = oldCq(:,(j-1)*2+2);
                [tempU{j} tCq(:,j)] = intersecter(qU1,qU2,q1,q2);
            end
        end
        oldCq = tCq;
        oldU  = tempU;
        old_num_cliq_level = num_cliq_level;
    end    
    U = oldU{1};
end

if size(U,1) ~= size(X,2) && type == 0
    error('U is incomplete');
end

%*************************************************************************
% Subfunction written by Nathan Krislock
%
%
%*************************************************************************
function [U,s] = subspace_intersection(tU1,tU2,k)
% [U,s] = subspace_intersection(U1,U2)
%
% Computes U satisfying
%
%       range(U) = range(U1) intersect range(U2).
%
% s : a vector containing the cosines of the principle angles between the
% subspaces range(U1) and range(U2).
%
% This uses the subspace intersection from Golub / Van Loan.



len_1 = size(tU1,1);
len_2 = size(tU2,1);

U1 = blkdiag(tU1,eye(len_2-k));
U2 = blkdiag(eye(len_1-k),tU2);

% U1 = orth(U1);
% U2 = orth(U2);

if size(U1'*U1,1) ~= size(eye(size(U1,2)),1)
    error('Size of U1 is incorrect')
end

if norm(U1'*U1 - eye(size(U1,2))) > size(U1,2)*0.001
    error('Faces must be orthogonal')
end

C = U1'*U2;
[Uc,Sc,Vc]= svd(C); %#ok<NASGU>

s = diag(Sc);


U = U1*Uc;
U = U(:,s > 1 - 1e-3);
%*************************************************************************
%*************************************************************************
%
% Subfunction for changing the permutation
%
%*************************************************************************
function [U w] = intersecter(U1,U2,v1,v2)

num = length(v1);
len_1 = sum(v1);   ind_1 = find(v1);
len_2 = sum(v2);   ind_2 = find(v2);

w1 = sparse(ind_1,ones(len_1,1),1:len_1,num,1);
w2 = sparse(ind_2,ones(len_2,1),1:len_2,num,1);

w  = v1 + v2; w(w > 0) = 1;
c  = find(v1.*v2);          % intersection of two cliques
k = sum(v1.*v2);            % size of intersection

% b1: atoms only in clique 1
b1 = max(v1 - v2,0);    b1 = find(b1);  ind_top_1 = w1(b1); %#ok<*FNDSB>
% b2: atoms only in clique 2
b2 = max(v2 - v1,0);    b2 = find(b2);  ind_bot_2 = w2(b2);

ind_bot_1 = w1(c);
ind_top_2 = w2(c);

tU1 = [U1(ind_top_1,:); U1(ind_bot_1,:)];
tU2 = [U2(ind_top_2,:); U2(ind_bot_2,:)];

%pU1 = blkdiag(tU1,eye(len_2-k));
%pU2 = blkdiag(eye(len_1-k),tU2);

%U = face_intersect(tU1,tU2,k);
U = subspace_intersection(tU1,tU2,k);


raw_index = [b1; c; b2];
[~, index_sorted] = sort(raw_index);
U = U(index_sorted,:);
%*************************************************************************


