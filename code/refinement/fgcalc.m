% BFGS objective and gradient finder for
% HANSO
%
% Babak Alipanahi
% University of Waterloo
%
% August 19, 2011

function [o g] = fgcalc(x,pars)

dim = pars.dim;
w = pars.w;
E = pars.equality_cons;
L = pars.lo_bounds;
U = pars.up_bounds;
f = pars.f;

X = reshape(x,dim,numel(x)/dim);

o = objfunmex(X,E,L,U,w,f);
g = gradfunmex(X,E,L,U,w,f);
g = g(:);