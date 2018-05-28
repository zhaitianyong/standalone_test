% This function finds PCA Rank
%
% Babak Alipanahi
% University of Waterloo
% July 29, 2010

function out = pcarank(X)

[~, sv] = pcab(X);
out = numel(sv);