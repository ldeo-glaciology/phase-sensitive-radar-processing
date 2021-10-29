function fmcw_error_model_plot

% Plot errors calculated by fmcw_error_model

n = 100;
rho_dB = [-2 -75];
phi_lim = (pi/180)*15; % 15 degrees
[R,Pr,Pref,SNR_dB,Rer,phi_rms] = fmcw_error(rho_dB(1));
Rer100 = Rer/sqrt(n);
gi = phi_rms<=phi_lim;
h1 = semilogy(R(gi),Rer(gi),'r');
hold on
gi = (phi_rms/sqrt(n))<=phi_lim;
h2 = plot(R(gi),Rer100(gi),'r','linestyle',':')
[R,Pr,Pref,SNR_dB,Rer,phi_rms] = fmcw_error(rho_dB(2));
Rer100 = Rer/sqrt(n);
gi = phi_rms<=phi_lim;
h3 = plot(R(gi),Rer(gi),'b');
hold on
gi = (phi_rms/sqrt(n))<=phi_lim;
h4 = plot(R(gi),Rer100(gi),'b','linestyle',':')

legend({'Base (n=1)','Base (n=100)','Internal (n=1)','Internal (n=100)'})

%legend( ['Base, n = 1, R\prime = ' int2str(round(rho_dB(1))) ' dB'],...
%        ['n = 100, R\prime = ' int2str(round(rho_dB(1))) ' dB'],...
%        ['n = 1, R\prime = ' int2str(round(rho_dB(2))) ' dB'],...
%        ['n = 100, R\prime = ' int2str(round(rho_dB(2))) ' dB'])

xlabel('range (m)')
ylabel('Measurement rms error (m)')
set(gca,'YMinorGrid','off')
set(gcf,'tag','FMCw_range_error_internal')
%ylim([0 100])
xlim([0 2e3])
ylim([1e-8 0.1])
grid on
    
function [h1,h2] = plot_rer(R,Rer,n,varargin)
% Range error
%figure
hold on
h1 = semilogy(R,Rer,'-',varargin{:});
hold on
h2 = plot(R,Rer/sqrt(n),':',varargin{:})


return
% 
% function plot_epl(R,Pr,Pref)
% %% Plot
% % Echo power level
% figure
% semilogx(R,10*log10(Pr/Pref))
% xlabel('range (m)')
% ylabel('Echo power (dBm)')
% ylim([-160 20])
% xlim([1 1e4])
% grid on
% 
% function plot_snr(R,SNR_dB)
% % SNR
% figure
% plot(R,SNR_dB)
% xlabel('range (m)')
% ylabel('SNR (dB)')
% ylim([0 100])
% xlim([0 2e3])
% grid on
