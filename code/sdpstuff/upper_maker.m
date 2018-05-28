% This function reads the hydrogen bond file in
% .hbo format and writes restraints in .upl format
%
% Babak Alipanahi
% University of Waterloo
% June 23, 2010
% -edit: February 15, 2011

function upper_bounds = upper_maker(dist_con,Comp)

max_i = numel(dist_con);
upper_bounds = nan(max_i,4);

max_res = max(Comp.residue);
min_res = min(Comp.residue);

for i =1:max_i
    sres  = dist_con(i).sres;
    tres  = dist_con(i).tres;
    
    if sres <= max_res && sres >= min_res && tres <= max_res && tres >= min_res        
        satom = dist_con(i).satom;
        tatom = dist_con(i).tatom;
        
        index_s = Comp.residue == sres;
        index_t = Comp.residue == tres;
        satom_names = Comp.atom_names(index_s);
        tatom_names = Comp.atom_names(index_t);
        
        local_s = index_finder(satom_names,satom);
        local_t = index_finder(tatom_names,tatom);
        
        global_s = find(index_s,local_s);
        global_t = find(index_t,local_t);
        
        upper_bounds(i,1) = global_s(local_s);
        upper_bounds(i,2) = global_t(local_t);
        
        if isnan(upper_bounds(i,1)*upper_bounds(i,2))
            error('unkown atom name');
        end
        upper_bounds(i,3) = dist_con(i).dist;
        upper_bounds(i,4) = dist_con(i).peak;
    end
end

index_bad = isnan(sum(upper_bounds,2));
upper_bounds(index_bad,:) = [];


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

%     max_j = size(satom_names,1);
%     max_k = size(tatom_names,1);
%     for j = 1:max_j
%         if strcmp(satom_names{j},satom)
%             upper_bounds(i,1) = index_s(j);
%             break
%         end
%     end
%     for k = 1:max_k
%         if strcmp(tatom_names{k},tatom)
%             upper_bounds(i,2) = index_t(k);
%             break
%         end
%     end