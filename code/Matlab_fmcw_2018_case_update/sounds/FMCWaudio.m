%%%%%%% sonify your apres data %%%%%%%%
clear; close all;



%% load the data
% filename= '/Users/elizabeth/Documents/projects/ApRES/phase-sensitive-radar-processing/data/DATA2015-10-27-1531.dat';
% filename = '
% burstnum=1;
vdat = fmcw_load('/Users/elizabeth/Documents/projects/juneau_icefield/2019/apres/Jonny_093/DIR2018-07-28-0817/DATA2018-07-28-0927.DAT');
disp(vdat.filename)

Fs = vdat.fs;

% split by attenuator
vdats = fmcw_burst_split_by_att(vdat);

% averages chirp
attnset=1; %only works for 1 setting for now (setting 1-4)
vdat = fmcw_burst_mean(vdats(attnset));

% phase process data
%[rc,rf,sr,su] = fmcw_range(vdat,p,maxrange,win);
%amp = 20*log10(abs(su));

t = vdat.t;
voltage = vdat.vif(attnset,:);

%normalizes voltage around 0 
voltage = voltage - mean(voltage);
voltage = voltage/max(abs(voltage));

rsv = resample(voltage,3,1);

soundsc(voltage, Fs)

audiowrite("JuneauSonify.wav",voltage,Fs);

%%
% set a sampling frequency
 % (= the default)
% 
% % make a time array
% sample = 1:42000;
% 
% % make a sin wave with the pitch of middle C
% C = sin(2*pi*261/Fs*sample);
% 
% % play middle C 
% soundsc(C,Fs)
% % 
% % % play a gong sound
% load gong.mat;
% soundsc(y)