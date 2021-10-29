function vdat = fmcw_burst_mean(vdat)

% fmcw_burst_mean(vdat)
%
% Takes mean of fmcw_burst
%
% Craig Stewart
% 2014/5/2
% 2014/5/20 - changed mean to first dimension fixing bug for single chirp
% burst

% Check if there are more than one attenuator settings
if length(unique(vdat.chirpAtt)) > 1
    error('Trying to average across multiple attenuator settings')
end

% Only do mean if there is more than one chirp in the burst
if size(vdat.vif,1)>1
    vdat.vif = mean(vdat.vif,1);
    vdat.ChirpsInBurst = 1;
    %vdat.chirpNum = nan;
    vdat.chirpTime = mean(vdat.chirpTime);
    vdat.chirpAtt = mean(vdat.chirpAtt);
    vdat.processing = [vdat.processing {'burst mean'}];
end
    