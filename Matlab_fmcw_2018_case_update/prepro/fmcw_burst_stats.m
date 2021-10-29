function S = fmcw_burst_stats(spec)

% Calculate statistical properties of the burst (ensemble of chirps)
% starting with the complex spectrum.
% 
% args:
% spec: complex spectrum size(nChirps,nSamples), each row representing a
%   chirp, each column a range.
%
% output: structure S with fields:
% .mean: mean spectrum (over all chirps in burst)
% .std: standard deviation
% .stdPhase: phase standard deviaion
% .stdSNR: signal to noise power ratio (chirp)
%
% The following fields a relevant to the burst mean
% .n: number of chirps in burst
% .ste: standard error
% .stePhase: phase standard error
% .steSNR: signal to noise power ratio

% 2014/10/24
% Craig Stewart

% Statistical properties of the chirps
S.n = size(spec,1); % number of chirps in burst
S.mean = mean(spec,1); % coherent average of chirps
S.std = std(spec,[],1); % RMS noise
S.stdPhase = S.std./(sqrt(2)*abs(S.mean)); % phase noise (valid for noise <~15deg)
S.stdSNR = (abs(S.mean)./S.std).^2; % signal to noise power ratio

% Statistical properties of the burst mean
% based on standard error with denominator sqrt(S.n-1)
S.ste = S.std./sqrt(S.n-1); % Standard error (expected rms noise on S.mean)
S.stePhase = S.ste./(sqrt(2)*abs(S.mean)); % standard error of phase
S.steSNR = (abs(S.mean)./S.ste).^2; % signal to noise power ratio