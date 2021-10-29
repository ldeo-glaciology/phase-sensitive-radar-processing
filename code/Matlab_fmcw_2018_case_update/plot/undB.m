function x = undB(dBx)

% Convert dBx from decibells to ratio
%
% Craig Stewart
% 2014/6/2

x = 10.^(dBx/20);