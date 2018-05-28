% This function randomly samples phi and psi values
%
% Babak Alipanahi
% University of Waterloo
% July 16, 2010
% -edit: July 29, 2010

function [phi psi] = ang_sampler(seq,phi_con,psi_con)

len = size(phi_con,1);
phi = zeros(len,1);
psi = zeros(len,1);
for i = 1:len
    if seq(i) == 15
        phi(i) = -60;
        psi(i) = -50;
    elseif i < len
        if seq(i+1) == 15
            phi(i) = -120;
            psi(i) = 120;            
        else
            phi_width = phi_con(i,2) - phi_con(i,1);
            phi(i) = phi_width*rand + phi_con(i,1);
            
            psi_width = psi_con(i,2) - psi_con(i,1);
            psi(i) = psi_width*rand + psi_con(i,1);
        end
    end
end
