function fmcw_showpro(vdat)

% fmcw_showpro(vdat)
%
% Display processing history for FMCW data in vdat

% Craig Stewart
% 2014/5/30

disp(' ')
disp('Processing record')
for ii = 1:length(vdat.processing)
    disp(vdat.processing{ii})
end
