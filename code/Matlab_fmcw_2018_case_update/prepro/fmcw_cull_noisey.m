function vdat = fmcw_cull_noisey(vdat,n)

% Cull n noisest chirps from a FMCW radar burst
%
% args: vdat = fmcw data structure
% n = number of noisest chirps to cull
%
% Craig Stewart
% 2014/5/20

% Check if there are more than one attenuator settings
if length(unique(vdat.chirpAtt)) > 1
    error('Multiple attenuator settings present in burst - subset into attenuator settings before culling noisy chirps')
end

% Check we're not killing all the data
v = vdat.vif;
nchirps = size(v,1);
if n >= nchirps
    disp([mfilename ' cancelled  - to few chirps to cull'])
    vdat.processing = [vdat.processing {[mfilename ': cancelled - to few chirps to cull']}];
    return
end

% Find noisy shots
mv = repmat(mean(v,1),size(v,1),1); % mean of all chirps
p = sqrt(mean((v-mv).^2,2)); % rms difference from mean for each shot
[~,ii] = sort(p,1,'ascend'); % ii gives the row numbers of v sorted by rms difference from mean (descending)


chirpsToKeep = ii(1:nchirps-n);
chirpsToCull = ii(1+nchirps-n:end);
vdat.processing = [vdat.processing {[mfilename ': removed noisest ' int2str(n) ' chirps: ' mat2str(vdat.chirpNum(chirpsToCull))]}];
chirpsToKeep = sort(chirpsToKeep);
vdat = fmcw_burst_subset(vdat,chirpsToKeep);