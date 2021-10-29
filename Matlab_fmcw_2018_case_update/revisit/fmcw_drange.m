function c = fmcw_drange(vdat,p,maxrange,winfun,Bf)

% c = fmcw_drange(vdat,maxrange,winfun,Bf)
%
% Calculates range using Brennan 2013 technique modified to use
% differential signal.
% note: Bf is the fractional bandwidth to achieve
%
% (1) Split each chirp into to by subsampling - odds (a) and evens (b)
%   and crop the end of (a) and the start of (b) so that their length is the
%   same but centre frequency different by fd. Redefine Fs and T.
% (2) Window, rotate and FFT a and b.
% (3) Create c = (conj(a).*b). c now has amplitude prod and phase difference
% of a and b (phase diff in direction a-> b). Effective wavelength = 2pi
% ci/delta_w
% (4) Subtract bin centre phase (multiply by reference vector), note this
% is not the same as the ref in Brennan 2013
% (5) Calculate coarse and fine range as before
%
% Craig Stewart
% 2014/8/21

if nargin < 2
    p = 1; % padding not required for this method
end
if nargin < 3
    maxrange = 2000; %m (range to crop output to)
end
if nargin < 4
    winfun = @blackman; % default to blackman window
end
if nargin < 5
    Bf = 4; % default to blackman window
end
%Bf = 4; % effective fractional bandwidth

%% Config
% (number of range bins over which we want a 2pi phase change) i.e. 400%
[nchirps,N] = size(vdat.vif);
if mod(N,2)
    error('N should be even');
end
n = round((N/2)*(1-(1/(Bf+1)))); % number of samples in subchirp
if ~mod(n,2)
    n = n+1; % make sure odd number of samples in each subchirp
    % (for phase centering)
end
ia = 1:2:2*n; % odds from start (length n)
ib = N-2*(n-1):2:N; % evens to end (length n), (N must be even)
%ia = 1:N-1; % test only
%ib = 1:N-1; % test only

%% (1) Split chirps
[a,b] = deal(vdat);

% a (odds)
a.vif = a.vif(:,ia); % subsample and crop data
a.f = a.f(ia); % frequency of Tx
a.f0 = a.f(1); % start frequency
a.fs = vdat.fs/2; % sampling frequency
a = fmcw_derive_parameters(a);
a.t = a.dt*(0:n-1); % sampling times (rel to first)
a.processing = [a.processing, {'Chirp split - first half'}];

% b (evens)
b.vif = b.vif(:,ib); % subsample and crop data
b.f = b.f(ib); % frequency of Tx
b.f0 = b.f(1); % start frequency
b.fs = vdat.fs/2; % sampling frequency
b = fmcw_derive_parameters(b);
b.t = b.dt*(0:n-1); % sampling times (rel to first)
b.processing = [b.processing, {'Chirp split - second half'}];

fd_target = vdat.B/(Bf+1); % centre frequency difference
% between sub chirps a and b that we're aiming for
fd = b.fc-a.fc; % frequency difference achieved
wd = 2*pi*fd; % angular velocity of above

%% (2) Window and FFT
fftmode = 'internal';  % 'internal' 'test'
% The reason for these two modes is that here we require part of the
% processing normally performed by fmcw_range - but only through to the fft
% - we don't require phase correction and calculating the fine range -
% unless we are testing the effect of changing the centre frequency
% wavelength on the determined range
switch fftmode
    case 'internal'
        % This only calculates the raw spectrum - saving doing the phase
        % correction and range calculation performed in fmcw_range
        [a.f,a.specRaw] = winfft(a,p,winfun);
        [b.f,b.specRaw] = winfft(b,p,winfun);
        if a.f~=b.f
            error('a.f~=b.f')
        end
        
    case 'test'
         % Call fmcw_range so that in addition to specRaw we also get
         % rangeFine to check bed depth. 
        [a.rangeCoarse,a.rangeFine,a.specCor,a.specRaw] = fmcw_range(a,p,maxrange,winfun);
        [b.rangeCoarse,b.rangeFine,b.specCor,b.specRaw] = fmcw_range(b,p,maxrange,winfun);
        
        a.range = repmat(a.rangeCoarse,nchirps,1) + a.rangeFine;
        [~,a.bedInd] = max(mean(abs(a.specRaw),1));
        a.bedRange = a.range(:,a.bedInd);
        b.range = repmat(b.rangeCoarse,nchirps,1) + b.rangeFine;
        [~,b.bedInd] = max(mean(abs(b.specRaw),1));
        b.bedRange = b.range(:,b.bedInd);
        dBedRange = b.bedRange-a.bedRange;
        dRange = mean(b.range,1)-mean(a.range,1);
        
        disp(' ')
        disp(['sub-chirp A bed range'])
        disp(['Coarse: ' num2str(a.rangeCoarse(a.bedInd),'%15.10f') ' m'])
        disp(['Fine  : ' num2str(a.rangeFine(a.bedInd),'%15.10f') ' m'])
        disp(['Total : ' num2str(mean(a.bedRange),'%15.10f') ' +- ' num2str(std(a.bedRange),'%15.10f') ' m'])
        disp(['Std er: ' num2str(std(a.bedRange)./sqrt(nchirps),'%15.10f') ' m'])
        disp(' ')
        disp(['sub-chirp B bed range'])
        disp(['Coarse: ' num2str(b.rangeCoarse(b.bedInd),'%15.10f') ' m'])
        disp(['Fine  : ' num2str(b.rangeFine(b.bedInd),'%15.10f') ' m'])
        disp(['Total : ' num2str(mean(b.bedRange),'%15.10f') ' +- ' num2str(std(b.bedRange),'%15.10f') ' m'])
        disp(['Std er: ' num2str(std(b.bedRange)./sqrt(nchirps),'%15.10f') ' m'])
        disp(' ')
        disp(['Range difference B-A: ' num2str(mean(dRange)) ' +- ' num2str(std(dRange)) ' m'])
        disp(['Bed range difference B-A: ' num2str(mean(dBedRange)) ' +- ' num2str(std(dBedRange)) ' m'])
        
        figure
        ax(1) = subplottight(4,1,1); % amp
        ha = plot(a.rangeCoarse,dB(abs(a.specRaw)),'b');
        hold on
        hb = plot(b.rangeCoarse,dB(abs(b.specRaw)),'r');
        legend([ha(1) hb(1)],{'a','b'})
        ylabel('amplitude (dB)')
        
        ax(2) = subplottight(4,1,2); % phase
        plot(a.rangeCoarse,angle(a.specRaw),'b');
        hold on
        plot(b.rangeCoarse,angle(b.specRaw),'r');
        ylabel('phase')
        
        ax(3) = subplottight(4,1,3); % fine range 
        plot(a.rangeCoarse,a.rangeFine,'b');
        hold on
        plot(b.rangeCoarse,b.rangeFine,'r');
        ylabel('fine range (m)')
        
        ax(4) = subplottight(4,1,4); % range diff
        plot(a.rangeCoarse,dRange,'.')
        xlabel('coarse range (m)')
        ylabel('Range difference B-A (m)')
        linkaxes(ax,'x')
        
        % The problem with this is that due to frequency dependent attenuation
        % in the system and bed roughness we have a phase gradient across the 
        % bed reflector. No this shouldn't be a problem as the bin
        % locations don't depend on centre frequency - only on bandwidth
        % I don't know whats going on here...
        keyboard
end

%% (3) Create c
c.specRaw = conj(a.specRaw).*b.specRaw;
c.lambdac = 2*pi*a.ci/wd; % effective wavelength of phase differenced spectrum
c.B = a.B;
c.fc = fd; % effective centre frequency is frequency difference

%% (4) Subtract bin centre phase
% Calculate phase of each range bin centre for correction
n = (0:size(c.specRaw,2)-1);
phiRef = 2*pi*fd*n./(c.B*p); % phase for each range bin centre (measured at t=T/2), given that tau = n/B
phiRef = repmat(phiRef,nchirps,1); % 1 chirp per row
compRef = exp(-1i*phiRef); % unit phasor with conjugate of above phase
c.specCor = compRef.*c.specRaw; % positive frequency half of spectrum with ref phase subtracted

%% (5) Calculate coarse and fine range as before
c.rangeCoarse = a.f*a.ci*a.T/(2*a.B); % Range at the centre of each range bin: eq 14 (rearranged) (p is accounted for in f)
%c.rangeFine = lambdac*angle(c.specCor)/(4*pi); % Distance from centre of range bin to effective reflector: eq 15
%c.rangeFine = angle(c.specCor)./((4*pi/lambdac) - (4*c.rangeCoarse*K/ci^2)); % this is the full equation including the term generated by the last term in (13)
c.rangeFine = fmcw_phase2range(angle(c.specCor),c.lambdac);
%R = c.rangeCoarse + c.rangeFine;

%% Crop output variables to useful depth range only
n = find(c.rangeCoarse<=maxrange,1,'last');
c.rangeCoarse = c.rangeCoarse(1:n);
%c.rangeCoarse = repmat(c.rangeCoarse,nchirps,1); % make output same size as c.rangeFine for consistence
c.rangeFine = c.rangeFine(:,1:n);
c.specRaw = c.specRaw(:,1:n);
c.specCor = c.specCor(:,1:n);

function [f,specRaw] = winfft(v,p,winfun)
% This calculates the windowed fft of the timeseries in p
% this is a cutdown version of fmcw_range (minus the phase correction) but
% it also avoids the loop through all chirps by processing the vdat.vif
% array complete.

% Processing settings
[nchirps,N] = size(v.vif);
xn = 0.5*(N-1); % timeshift amount (steps) prior to fft to get fft measuring phase at t=T/2
if ~mod(N,2)
    disp('Warning: even number of sample in record, not possible to correctly offset to phase centre')
    xn = round(xn);
end

% Measure the sampled IF signal: FFT to measure frequency and phase of IF
deltaf = 1/(v.T*p); % frequency step of FFT
f = 0:deltaf:v.fs/2-deltaf; % frequencies measured by the fft - changed 16 April 2014, was %f = [0:deltaf:fs/2];
nf = length(f); % changed from above 2014/5/22
win = repmat(window(winfun,N),1,nchirps); %chebwin(N);  %rectwin(N); %

vif = v.vif'; % 1 chirp per column
vif = vif-repmat(mean(vif,1),N,1); % de-mean each chirp
vif = win.*vif; % window each chirp
vifpad = zeros(p*N,nchirps);
vifpad(1:length(vif),:) = vif;
%vifpad = circshift(vifpad,-xn,1); % rotate so phase centre at start (this is version for 2014b)
vifpad = circshift(vifpad,[-xn 0]); % rotate so phase centre at start (2012b)
%imagesc(abs(vifpad)), keyboard
fftvif = (sqrt(2*p)/length(vifpad)).*fft(vifpad); % fft and scale for padding
fftvif = fftvif./repmat(rms(win,1),N*p,1); % scale for window
specRaw = fftvif(1:nf,:); % positive frequency half of spectrum up to (nyquist minus deltaf)
specRaw = transpose(specRaw);


