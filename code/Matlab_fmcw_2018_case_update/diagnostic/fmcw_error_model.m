function [R,Pr,Pref,SNR_dB,Rer,phi_rms] = fmcw_error(rho_dB)

% [R,Pr,Pref,SNR_dB,Rer,phi_rms] = fmcw_error(rho_dB)
%
% Calculate theoretical range error as a function of system parameters and range
% using equation (8) of Brennan et al 2013

% Craig Stewart
% 2014/8/14

if nargin == 0
    reflector_type = 'internal'; % 'internal' 'base'
    switch reflector_type
        case 'internal'
            rho_dB = -90; % Internal loss
        case 'base'
            rho_dB = -2; % Bed loss
    end
    n = 100; % number of chirps to average
end
R = [1:10:2000]; % range

%% System parameters
Pt_dBm = 20; % Power transmitted, rel to 1mW at 50 ohms = 10*log10(Pt/Pref)
Pref = 0.001; % watts
Pt = Pref*10^(20/10); % (watts)
G_dB = 10; % Antenna gain (dBi)
G = 10^(G_dB/10); % Antenna gain
dR = 0.43; % m % Depth resolution
T = 1; % 1 second pulse duration
df = 1/T; % bandwidth of frequency bin
fc = 3e8; % centre frequency

%% Environmental parameters
c = 3e8; % speed of light
epsilon = 3.1; % electrical permittivity of ice
lambda = c/(fc*sqrt(epsilon));
L_m_dB = -0.015; % Path Loss dB/m
L_dB = L_m_dB.*2*R; % Total path loss dB
L = 10.^(L_dB./10); % Total Path Loss
rho = 10^(rho_dB/10);

%% Power received
Pr = Pt*G^2*rho*L*lambda^2*dR./(16*pi^2*R.^3); %(watts) % Brennan equation (5)

%% Receiver noise
k = 1.3806488e-23; % Boltzmans constant (used for thermal noise)
%Tc = 20; % Temperature of receiver (deg C)
%Ta = Tc + 273.15; % Absolute temperature of receiver (Kelvin)
Ta = 300; % nominal to match Brennan
Ptn = k*Ta*df; % thermal noise (watts)
NF = 6; % Receiver noise figure (dB) = 10*log10(noise_factor)
F = 10^(NF/10); % Receiver noise factor
Pn = Ptn*F; % Receiver noise power out (watts)

%% Signal to noise ratio
SNR = Pr/Pn;
SNR_dB = 10*log10(SNR);

% Phase noise
phi_rms = 1./sqrt(2*SNR); % or should this have 2pi on top (i.e. is this refering to fractions of a cycle)
%phi_rms = (2*pi)./sqrt(2*SNR);
Rer = lambda*phi_rms/(4*pi);

return

%% Plot
% Echo power level
figure
semilogx(R,10*log10(Pr/Pref))
xlabel('range (m)')
ylabel('Echo power (dBm)')
ylim([-160 20])
xlim([1 1e4])
grid on

% SNR
figure
plot(R,SNR_dB)
xlabel('range (m)')
ylabel('SNR (dB)')
ylim([0 100])
xlim([0 2e3])
grid on

% Range error
figure
h1 = semilogy(R,Rer,'b');
hold on
h2 = plot(R,Rer/sqrt(n),'r')
legend(['Single chirp (power reflection coefficient: ' int2str(round(rho_dB)) ' dB'],[int2str(n) ' chirp average'])
xlabel('range (m)')
ylabel('Measurement rms error (m)')
set(gca,'YMinorGrid','off')
set(gcf,'tag','FMCw_range_error_internal')
%ylim([0 100])
xlim([0 2e3])
grid on

