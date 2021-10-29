%function fmcw_demo_clean(filename)

% Clean up a file plotting along the way
%
% Craig Stewart
% 2014/5/20

%if nargin ==0
%    filename = uigetfile({'*.dat;*.DAT;*.000;*.mat','Radar files: .dat, .DAT, .000 and .mat'},'Choose a radar files to demo clean');
%end

filename = ('a02_2014-01-12_213419.dat');

% Load
vdat = fmcw_load(filename,1);

% Cull f range
frange = [2.5987456e8 3.230498e8];
disp(['Culling frequency range to: ' num2str(frange(1)/1e6) '-' num2str(frange(2)/1e6) 'MHz'])
vdat2 = fmcw_cull_freq(vdat,frange);
%vdat2.processing

% Cull gross contaminated chirps
noisePowerLimit = 0.01;
disp(['Removing grossly contaminated shots'])
vdat3 = fmcw_cull_bad(vdat2,noisePowerLimit);

% Cull noisey remaining
n = 10;
disp(['Culling noisest ' int2str(n) ' shots'])
vdat4 = fmcw_cull_noisey(vdat3,n);

% Display processing
fmcw_showpro(vdat4)

% Plot
figure
a = plot(vdat.t,transpose(vdat.vif),'b');
hold on
b = plot(vdat2.t,transpose(vdat2.vif),'r');
c = plot(vdat3.t,transpose(vdat3.vif),'c');
d = plot(vdat3.t,transpose(vdat4.vif),'g');
legend([a(1) b(1) c(1) d(1)],{'raw','frequency cropped','bad shots removed',['remaining noiseyest ' int2str(n) ' shots removed']})