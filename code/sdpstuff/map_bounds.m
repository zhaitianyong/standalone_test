% This program convert constraints from complete
% protein to reduced-atom set protein
%
% Babak Alipanahi
% University of Waterloo
% February 10, 2011

function bounds = mapBounds(in_bounds,atoms_map)

[num_bounds num_cols] = size(in_bounds);
bounds = nan(num_bounds,num_cols);
for i = 1:num_bounds
    
    if in_bounds(i,1)*in_bounds(i,2) == 0
        error('There is a error in mapping atoms');
    end
    bounds(i,[1 2]) = atoms_map(in_bounds(i,[1 2]),2);    
    bounds(i,3:end) = in_bounds(i,3:end);
end

if any(isnan(bounds(:)))  
   error('There is an error in mapping atoms'); 
end