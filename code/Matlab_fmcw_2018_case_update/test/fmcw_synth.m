function vdat = fmcw_synth(R,A,N,S,adcBitDepth)

% Create synthetic FMCW radar data
%
% args:
% R = range vector (m)
% A = amplitude of reflectors at range R (dB)
% N = noise floor level (dB)
% S = scale factor (pre scale data to this fraction of the adc range to avoid clipping)
% adcBitDepth = bit depth of adc (default 16)
% 
% Craig Stewart
% 2013 April 24
% 2014/5/30 finished!

if nargin == 1
    A = ones(size(R));
end
if size(R) ~= size(A)
    error('Inputs must be the same size')
end
R = R(:);
A = A(:);
if nargin <3
    N = nan; % no noise
end
if nargin <4
    S = 0; % no scaling
end
if nargin <5
    adcBitDepth = 16;
end

% Constants
vdat.fs = 4e4; % sampling frequency
vdat.f0 = 2e8; % start frequency
vdat.K = 2*pi*2e8; % chirp gradient in rad/s/s (200MHz/s)
%vdat.f0 = vdat.f0 + (vdat.K/(4*pi))/vdat.fs; % start frequency
vdat.SamplesPerChirp = 40000;
vdat.er = 3.1;
%T = (SamplesPerChirp-1)*1/fs; % period between first and last sample
%T = SamplesPerChirp/fs; % period of sampling
%f1 = f0 + T*K/(2*pi); % stop frequency
%
%ci = 3e8/sqrt(er); % velocity in material
%fc = mean([f0 f1]); % Centre frequency
%t = repmat([0:1/fs:T-1/fs],numel(R),1); % sample times
vdat.vif = ones(1,vdat.SamplesPerChirp); % dummy so we can use fmcw_derive_parameters
vdat = fmcw_derive_parameters(vdat);
vdat.t = vdat.dt*(0:size(vdat.vif,2)-1); % sampling times (rel to first)
vdat.f = vdat.f0 + vdat.t.*vdat.K/(2*pi);


% fs = 4e4;
% T = 1-1/fs;
% f0 = 2e8;
% f1 = 4e8;
% B = f1-f0; % bandwidth (hz)
% K = 2*pi*B/T; % chirp gradient (rad/s/s)
% ci = 3e8/sqrt(3.1); % c ice
% fc = mean([f0 f1]);
% t = repmat([0:1/fs:T],numel(R),1); % sample times

% Create timeseries for each reflector (one row per reflector)
amp = repmat(10.^(A/20),1,size(vdat.t,2)); % volts
tau = repmat(2*R/vdat.ci,1,size(vdat.t,2)); % 2-way travel time
t = repmat(vdat.t,size(R,1),1);
VIF = amp.*cos(2*pi*vdat.fc*tau + vdat.K*tau.*(t-vdat.T/2) - vdat.K*tau.^2/2); % voltage at IF
% (Brennan et al 2013 equation 13.)
vif = sum(VIF,1); % Sum all the seperate reflectors into a single timeseries

% Add noise
if ~isnan(N)
    noise = randn(size(vif));
    noise = 100 * (noise/rms(noise)) * 10.^(N/20);
    %noise = (noise/rms(noise)) * 10.^(N/20);
    % Note: the factor of 100 was experimental but matches levels ok
    % this is because we're then taking fft of 40K samples which cancells a
    % lot of noise - should be related to sqrt(num samples)?
    vif = vif + noise;
end

if S
    % Scale data
    vif = 0.5*vif/max(abs(vif)); % first to within +/- 0.5
    vif = S*vif; % further scale used defined to force clipping etc
    
    % Clip data to range
    vif(vif>0.5) = 0.5;
    vif(vif<-0.5) = -0.5;
end
vdat.Attenuator_1 = 0;
vdat.Attenuator_2 = 0;

% Digitise to 16 bits
if adcBitDepth < inf
    levels = 2^adcBitDepth; % ADC resolution
    vif = round(levels*vif)/levels; % vif digitized to 16 bit
end
% Rescale to ADCoutput voltage range
adcVRange = 2.5; % ADC input range
vif = adcVRange*(vif + 0.5); % vif no 0-adcVRange

% Create fmcw data structure with esential metadata
vdat.v = vif;
vdat.Startind = 1;
vdat.Endind = numel(vif);
vdat.ChirpsInBurst = 1;
vdat.TimeStamp = now;
vdat.Code = 10;

%eval(['save testdata_shot' int2str(ind) ' vdat '])