function vdat = fmcw_cal(vdat,calArray)

% fmcw_cal(vdat,ca)
% Calibrate shot using cal file
%
% System attenuation of the signal at high frequencies causes phase
% gradients. Correction of the attenuation in the raw signal removes this.
%
% args:
% vdat fmcw data array 
% ca: callibration array
%
% Craig Stewart
% 2014/10/15

% Define default frequency dependent amplitude correction
[nShots,nSamples] = size(vdat.vif);
if nargin==1
    calArray = [1 4];
end

% Make calibration array
cav = linspace(calArray(1),calArray(2),nSamples); % calibration vector
cam = repmat(cav,nShots,1); % calibration matrix (size of vdat.vif)

% Apply cal
vdat.vifOrig = vdat.vif; % keep copy of the raw data
vifMean = mean(vdat.vif(:));
vdat.vif = vifMean + cam.*(vdat.vif-vifMean); % apply cal
vdat.processing = [vdat.processing {[mfilename ': calibrated for system RF attenuation (scale range ' num2str(cav(1)) '-' num2str(cav(1)) ')']}];

% give a warning for now so we don't forget we're doing this
disp('Calibrating raw timeseries!')
