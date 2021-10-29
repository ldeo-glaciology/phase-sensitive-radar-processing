function vdat = fmcw_derive_parameters_update(vdat)

% vdat = fmcw_derive_parameters(vdat)
%
% Calculate fmcw radar parameters from fundamental parameters
% > size(vdat.vif)
% > fs            (sampling frequency)
% > f0            (start frequency)
% > K             (chirp gradient)
% > er            (material permitivity)

% Craig Stewart
% 2013-10-23
% 2014/5/20 removed vdat.t (to fmcw_load)
% 2014/5/21 changed how radar chirp is defined (now using chirp gradient K as
% fundamental parameter)

%vdat.ChirpsInBurst = size(vdat.vif,1); %already set in LoadBurstRMB5 and 6
%vdat.SamplesPerChirp = size(vdat.vif,2); %already set in LoadBurstRMB5 and 6

    vdats.dt = 1/vdats.fs; % sample interval (s)
    vdat.T = (size(vdat.vif,2)-1)/vdat.fs; % period between first and last sample
    %vdat.T = size(vdat.vif,2)/vdat.fs; % period of sampling (cls test 26 aug 2014)
    % - this makes the amplitude of teh fft centred at the right range, but phase wrong

    vdat.f1 = vdat.f0 + vdat.T*vdat.K/(2*pi); % stop frequency
    %vdat.f1 = vdat.f0 + vdat.dt*(vdat.SamplesPerChirp-1)*vdat.K/(2*pi); % stop frequency

    %vdat.B = vdat.f1-vdat.f0; % bandwidth (hz)
    %vdat.B = vdat.T*(vdat.K/(2*pi)); % bandwidth (hz)
    vdat.B = (size(vdat.vif,2)/vdat.fs)*(vdat.K/(2*pi)); % bandwidth (hz)

    vdat.fc = mean([vdat.f0 vdat.f1]); % Centre frequency
    %vdat.fc = vdat.f0 + vdat.B/2; % Centre frequency
    vdat.ci = 3e8/sqrt(vdat.er); % velocity in material
    vdat.lambdac = vdat.ci/vdat.fc; % Centre wavelength
