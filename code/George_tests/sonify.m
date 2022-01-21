%% load the data
filename= '../../data/DATA2021-12-13-1417.dat';
burstnum=1;
vdat = fmcw_load(filename, burstnum);

% split by attenuator
vdats = fmcw_burst_split_by_att(vdat);

% averages chirp
attnset=1; %only works for 1 setting for now (setting 1-4)
vdat = fmcw_burst_mean(vdats(attnset));

% phase process data
%[rc,rf,sr,su] = fmcw_range(vdat,p,maxrange,win);
%amp = 20*log10(abs(su));

%t = vdat.t;
voltage = vdat.vif(attnset,:);

%normalizes voltage around 0 
voltage = voltage - mean(voltage);
voltage = voltage/max(abs(voltage));
%%
Fs = 40000;
soundsc(voltage,Fs)
audiowrite('wall_sonify.wav',voltage,Fs);