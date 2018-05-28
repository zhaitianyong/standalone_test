% This program puts back hydrogens
%
% Babak Alipanahi
% University of Waterloo
% February 14, 2011

function oX = hydrogen_mapper(X,orig_rX,Comp,orig_Comp,A)

seq = Comp.seq;
num_res = length(seq);
num_atoms = size(orig_rX,2);
oX = nan(3,num_atoms);



for i = 1:num_res-1
    
    index_orig = find(orig_Comp.residue == Comp.num_seq(i));
    index_reduced = Comp.residue == Comp.num_seq(i);
    
    cur_num_atoms = length(index_orig);
    
    toX = nan(3,cur_num_atoms);
    tX = X(:,index_reduced);
    ind = Comp.atoms_map(index_orig,:);
    ind = ~isnan(ind(:,2));
    toX(:,ind) = tX;
    
    sX = orig_rX(:,index_orig);
    
    num_anchors = size(A(seq(i)).anchors,1);
    for j = 1:num_anchors
        temp_missing   = A(seq(i)).anchors{j,1};
        temp_available = A(seq(i)).anchors{j,2};
        
        refX = toX(:,temp_available);
        mapX = sX(:,temp_available);
        
        %         [err, ~, R] = procrustes(refX',mapX','scaling',false);
        %         %[err, ~, R] = procrustes(refX',mapX','scaling',false,'reflection',false);
        %         if err > 1e-2
        %             fprintf('Warning: error is too large %5.2f.\n',err);
        %         end
        %         tvec = repmat(R.c(1,:),size(sX(:,[temp_missing temp_available]),2),1);
        %         temp_toX = (sX(:,[temp_missing temp_available])'*R.T + tvec)';
        
        reflection = 1;
        [R t err]  = procrustesb(refX,mapX,reflection);
        
        if err > 5e-1
            fprintf('Warning: error is too large %5.2f.\n',err);
        end
        
        temp_toX = R*sX(:,[temp_missing temp_available]);
        temp_toX = bsxfun(@plus,temp_toX,t);
        
        toX(:,temp_missing) = temp_toX(:,1:length(temp_missing));
    end
    oX(:,index_orig) = toX;
end


