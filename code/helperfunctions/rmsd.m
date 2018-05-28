%===RMSD Calculator=======================
% Syntax: [rmsd_value]=rmsd(cord_1,cord_2)

function [rmsd_value] = rmsd(cord_1,cord_2)

diff_vector = (cord_1-cord_2).^2;
rmsd_value = sqrt(mean(diff_vector));