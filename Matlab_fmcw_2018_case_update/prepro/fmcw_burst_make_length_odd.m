function vdat = fmcw_burst_make_length_odd(vdat)

% Drops last sample from each chirp
%
% Craig Stewart
% 2014/8/21

% Ensure that number of samples is odd (so that we can shift the data to
% get the phase centre at the beginning of the sequence prior to fft)
n = size(vdat.vif,2); % number of samples per chirp
if ~mod(n,2)
    vdat.vif = vdat.vif(:,1:n-1); % remove last sample
    %disp('Cropping data to odd number of points per chirp')
    vdat.processing = {[mfilename ': cropped chirps to ' int2str(n-1) ' points']};
    vdat = fmcw_derive_parameters(vdat);
    vdat.t = vdat.t(1:n-1);
    vdat.f = vdat.f(1:n-1);
end
