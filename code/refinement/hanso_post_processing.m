% This function performs post processing on the
% SDP output using the HANSO method
%
% Babak Alipanahi
% University of Waterloo
%
% August 19, 2011


function [X info] = hanso_post_processing(ref_X,X0,Comp,lo_bounds,up_bounds,w,f)

equality_cons = equality_cons_former(ref_X,Comp,2);

pars.nvar = numel(X0);
pars.fgname = 'fgcalc';
% 
% f_hb  = f[0];
% f_ta  = f[1];
% f_vdw = f[2];
%f = [10 10 1];

ow = w;
w(4) = 0.0;

pars.equality_cons = equality_cons;
pars.up_bounds = up_bounds;
pars.lo_bounds = lo_bounds;
pars.dim = size(X0,1);

obj  = objfunmex(X0,equality_cons,lo_bounds,up_bounds,w,f);

w(4) = ow(4)*obj/(25*trace(X0'*X0));
pars.w = w;
pars.f = f;

options.x0 = X0(:);
options.prtlevel = 0;
options.maxit = 250;
options.normtol = 10^-9;
%options.quitLSfail = 0;
%[x f] = bfgs(pars,options);
[x, ~, ~, ~, ~, ~, ~, ~, ~, info.obj] = bfgs(pars,options);

max_i = numel(info.obj);
temp = nan(max_i,1);
for i = 1:max_i
   temp(i) = mean(info.obj{i}); 
end
info.obj = temp;

     
X = reshape(x,size(X0,1),size(X0,2));

