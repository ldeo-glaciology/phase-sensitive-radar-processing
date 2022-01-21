%% Measuring the change in range from 2 time seperated signals
% Define constants
B = 200e6; % Bandwidth (Hz)
c = 3e8; % Speed of light
eps_r = 3.1; % Dielectric constant of ice
T = 1; % Pulse duration
f_c = 300e6; % Center Frequency (Hz)

del_R = c/(2*B*sqrt(eps_r)); % Equation 1b

% Define an arbitrary reflector scenario
% We know del_R is ~43cm, so let's define something at a higher resolution
% Say there are 3 reflectors, one at 10m, one at 20m, and one at 20.2m
R1 = 10;
R2 = 20;
R3 = 21;
tau1 = 2*R1*sqrt(eps_r)/c; % from paragraph between eqns 12 and 13
tau2 = 2*R2*sqrt(eps_r)/c; 
tau3 = 2*R3*sqrt(eps_r)/c; 

% Now we generate our first signal
t = linspace(0,1,40000); % initialize t array (40kHz sampling)
K = 2*pi*B/T; % eqn 11
f1 = K*tau1;
f2 = K*tau2;
f3 = K*tau3;
s1 = exp(1j*f1*t)+exp(1j*f2*t)+exp(1j*f3*t);

% shift the reflectors by a bit and generate second signal
R1 = 10.1;
R2 = 20.2;
R3 = 22;
tau1 = 2*R1*sqrt(eps_r)/c; % from paragraph between eqns 12 and 13
tau2 = 2*R2*sqrt(eps_r)/c; 
tau3 = 2*R3*sqrt(eps_r)/c; 

% Now we generate our second signal
t = linspace(0,1,40000); % initialize t array (40kHz sampling)
K = 2*pi*B/T; % eqn 11
f1 = K*tau1;
f2 = K*tau2;
f3 = K*tau3;
s2 = exp(1j*f1*t)+exp(1j*f2*t)+exp(1j*f3*t);

[rc1,rf1,a1]=get_range(s1);
[rc2,rf2,a2]=get_range(s2);
r1 = rc1+rf1;
r2 = rc2+rf2;

figure()
hold on;
plot(r1,abs(a1),'DisplayName','s1 range');
plot(r2,abs(a2),'DisplayName','s2 range');
xlim([8,24]);
xlabel("Range (m)"); ylabel("Amplitude")
legend;
%% Get range difference
% Try with coarse
[pks1,locs1]=findpeaks(abs(a1),rc1,'MinPeakProminence',1e-3);
[pks2,locs2]=findpeaks(abs(a2),rc2,'MinPeakProminence',1e-3);

shifts = locs2-locs1;
disp("Using coarse range, reflectors shifted by: "+shifts+" m");

% Try with coarse
[pks1,locs1]=findpeaks(abs(a1),r1,'MinPeakProminence',1e-3);
[pks2,locs2]=findpeaks(abs(a2),r2,'MinPeakProminence',1e-3);

shifts = locs2-locs1;
disp("Using full range, reflectors shifted by: "+shifts+" m");

function [Rcoarse,Rfine,amplitude] = get_range(s)
    B = 200e6; % Bandwidth (Hz)
    c = 3e8; % Speed of light
    T = 1; % Pulse duration
    eps_r = 3.1; % Dielectric constant of ice
    f_c = 300e6; % Center Frequency (Hz)
    K = 2*pi*B/T; % eqn 11
    N = size(s,2);
    xn = round(0.5*(N));
    [nchirps,N] = size(s); % only 1 chirp
    win = blackman(N); % define window
    p = 10; % define pad factor
    
    nf = round((p*N)/2 - 0.5); % number of frequencies to recover 
    % Calculate phase of each range bin centre for correction - not sure either
    n = (0:nf - 1)';
    Rcoarse =  n*c/(2*B*sqrt(eps_r)*p); % From figure 6
    phiref = 2*pi*f_c*n./(B.*p) - (K*n.^2)/(2*B.^2*p.^2); % eq 17: phase for each range bin centre (measured at t=T/2), given that tau = n/(B*p)
    [spec,spec_cor] = deal(zeros(nchirps,nf)); % preallocate
    
    vif = s;
    vif = vif-mean(vif); % de-mean
    vif = win.*vif.'; % windowed
    
    % zero padded to length p*N
    vifpad = zeros(p*N,1);
    vifpad(1:length(vif)) = vif;
    vifpad = circshift(vifpad,-xn); % signal time shifted so phase centre at start
    
    fftvif = (sqrt(2*p)/length(vifpad)).*fft(vifpad); % fft and scale for padding 
    fftvif = fftvif./rms(win); % scale for window
    spec = fftvif(1:nf); % positive frequency half of spectrum up to (nyquist minus deltaf)
    comp = exp(-1i*(phiref)); % unit phasor with conjugate of phiref phase
    spec_cor = comp.*fftvif(1:nf); % positive frequency half of spectrum with ref phase subtracted
    lambdac = 0.5608;
    Rfine = angle(spec_cor)./((4*pi/lambdac) - (4*Rcoarse*K/(c/sqrt(eps_r))^2)); 
    % Trim down to relevant sections
    Rfine = Rfine(1:150);
    Rcoarse = Rcoarse(1:150);
    amplitude = spec_cor(1:150);
end