function vdat = fmcw_burst_subset(vdat,chirplist)

% vdat = fmcw_burst_subset(vdat,chirplist)
%
% Keep only requested chirps from bursts
%
% Craig Stewart
% 2014/5/2

%chirplist = sort(chirplist);

% crop all chirp specific values to this selection
vdat.vif = vdat.vif(chirplist,:); % data
vdat.chirpNum = vdat.chirpNum(chirplist);
vdat.chirpAtt = vdat.chirpAtt(chirplist);
vdat.chirpTime = vdat.chirpTime(chirplist);
vdat.ChirpsInBurst = size(vdat.vif,1); % should this be left as a record of the original burst length?
attSetList = unique(vdat.chirpAtt,'stable');
vdat.Attenuator_1 = real(attSetList);
vdat.Attenuator_2 = imag(attSetList);
vdat.processing = [vdat.processing {[mfilename ': keeping chirps ' mat2str(vdat.chirpNum)]}];