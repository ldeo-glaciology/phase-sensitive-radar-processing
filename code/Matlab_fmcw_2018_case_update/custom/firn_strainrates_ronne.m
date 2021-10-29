
clear all;
close all;
%% extract file names
fd = dir('/Users/elizabeth/Google Drive/home_research/projects/ApRES/data/Ronne_ApRES_KN/Site5cApRESdata 3/*.DAT');

%% choose earliest date as baseline (filename1) for config file

% 1
% 2
f1 = 'DATA2015-01-10-1448.DAT';

%% compare each file to baseline and record time, vsr, vsre, cfg, f (or g) .rangecoarse, f, g

for i = 2:length(fd)
    fn = fd(i).name;
    site = fmcw_dmelt(fn);

    ronne_vsr(i).t1 = site.t1; % time
    ronne_vsr(i).t2 = site.t2; % time
    ronne_vsr(i).dt = site.dt;
    ronne_vsr(i).file1 = site.file1;
    ronne_vsr(i).file2 = site.file2;
    ronne_vsr(i).fit = site.fit; % /year
    ronne_vsr(i).vsr = site.vsr; % /year
    ronne_vsr(i).vsre = site.vsre; % /year
    ronne_vsr(i).resid = site.resid;
    ronne_vsr(i).cfg = site.cfg;
    ronne_vsr(i).rangef = site.rangef;
    ronne_vsr(i).rangeg = site.rangeg;
    %ronne_vsr.f = f %will be the same for all
    %ronne_vsr.g = g
    
    clearvars -except fd ronne_vsr f1
    close all
end