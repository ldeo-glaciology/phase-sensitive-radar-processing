function vdat = fmcw_cull_freq(vdat,frange)

% Crop data from a FMCW radar burst to specified frequency range
% (i.e. a specific time range)
%
% args: vdat = fmcw data structure
% frange = 2 element vector or frequency range to crop to (Hz)
% 
% Craig Stewart
% 2014/5/20

% Calculate the original frequencies at each sample point
sampleNo = 1:vdat.SamplesPerChirp;

% Check cutoff frequencies are sensible
if (min(frange)<vdat.f(1)) || (max(frange)>vdat.f(end))
    if frange(1) < vdat.f(1)
        frange(1) = vdat.f(1);
        disp('Warning: low frequency cutoff below start frequency - not cropping low')
    end
    if frange(2) > vdat.f(end)
        frange(2) = vdat.f(end);
        disp('Warning: high frequency cutoff above stop frequency - not cropping high')
    end
end

% Crop frequency range
n1 = round(interp1(vdat.f,sampleNo,frange(1)));
n2 = round(interp1(vdat.f,sampleNo,frange(2)));
if n1 == n2
    error('No data left in frequency range')
end
vdat.vif = vdat.vif(:,n1:n2); % signal cropped
vdat.f0 = vdat.f(n1); % overwrite vdat.f0 in vdat
vdat.f1 = vdat.f(n2); % overwrite vdat.f0 in vdat
vdat.t = vdat.t(n1:n2); % keep original sampling times
vdat.f = vdat.f(n1:n2); % keep original sampling freqs

% Recalculate derived parameters from new start/end frequencies
vdat = fmcw_derive_parameters(vdat);

% Update processing record
vdat.processing = [vdat.processing {[mfilename ': cropped timeseries to frequency range: ' mat2str(frange/1e6) ' MHz']}];