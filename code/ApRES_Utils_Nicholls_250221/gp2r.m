function r = gp2r;

% Quick tool for getting range difference from fmcw phase with mouse click
% input
%
% Craig Stewart
% 2014/10/29

lambda = 1/sqrt(3.1); % = 0.568m fmcw wavelength
[~,y] = ginput(2);
r = lambda*diff(y)/(4*pi);