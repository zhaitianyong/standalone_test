% This function adjusts distance constraints
%
% Babak Alipanahi
% University of Waterloo
% February 23, 2011
% -edit: March 3, 2011

function oC = upperBound_adjuster_homitted(C,A)

oC = C;

num = size(C,2);
for i = 1:num
    satom  = C(i).satom;
    tatom  = C(i).tatom;
    stype  = C(i).stype;
    ttype  = C(i).ttype;
    dist   = C(i).dist;
    
    shoswaps = A(stype).hoswaps;
    thoswaps = A(ttype).hoswaps;
    max_s  = size(shoswaps,2);
    max_t  = size(thoswaps,2);
    
    delta_s = 0;
    for s = 1:max_s
        if strcmp(shoswaps{s}{1},satom)
            oC(i).satom = shoswaps{s}{2};
            inds = index_finder(A(stype).atom,{shoswaps{s}{1},shoswaps{s}{2}});
            delta_s = norm(A(stype).X(:,inds(1)) - A(stype).X(:,inds(2)));
            break
        end
    end
    delta_t = 0;
    for t = 1:max_t
        if strcmp(thoswaps{t}{1},tatom)            
            oC(i).tatom = thoswaps{t}{2};
            inds = index_finder(A(ttype).atom,{thoswaps{t}{1},thoswaps{t}{2}});
            delta_t = norm(A(ttype).X(:,inds(1)) - A(ttype).X(:,inds(2)));
            break
        end
    end    
    oC(i).dist = dist + delta_s + delta_t;
end


function ind = index_finder(atoms,names)

ind = nan(1,numel(names));
for j = 1:numel(names)
    for i = 1:numel(atoms)
        if strcmp(atoms(i).name,names{j})
            ind(j) = i;
            break
        end
    end
end