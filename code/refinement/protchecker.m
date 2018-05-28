% This function verifies the quality of the protein
%
% Babak Alipanahi
% University of Waterloo
% March 10, 2011

function report = protchecker(X,Comp,eq_cons,lo_bounds,up_bounds,print_flag)

if ~exist('print_flag','var')
    print_flag = 0;
end

report.eq_err = bondcheckmex(X,eq_cons);
report.lo_err = loboundmex(X,lo_bounds);
report.up_err = upboundmex(X,up_bounds);
report.chiral = chirality_check(X,Comp,0);
[report.phi report.psi] = ang_checker(X,Comp,0);

report.hbond = report.up_err(up_bounds(:,4) == -1);
report.dihed = report.up_err(up_bounds(:,4) == -2);


if print_flag > 0    
    fprintf('Distance  violations >\tpercent: %5.2f%%\t\tmean: %4.2f\n',100*sum(report.eq_err ~= 0)/numel(report.eq_err),mean(report.eq_err(report.eq_err > 0)));     
    fprintf('Lo bound  violations >\tpercent: %5.2f%%\t\tmean: %4.2f\n',100*sum(report.lo_err ~= 0)/numel(report.lo_err),mean(report.lo_err(report.lo_err > 0)));
    fprintf('Up bound  violations >\tpercent: %5.2f%%\t\tmean: %4.2f\n',100*sum(report.up_err ~= 0)/numel(report.up_err),mean(report.up_err(report.up_err > 0)));
    fprintf('Chirality violations >\tpercent: %5.2f%%\n',100*sum(report.chiral < 1)/numel(report.chiral));
end
