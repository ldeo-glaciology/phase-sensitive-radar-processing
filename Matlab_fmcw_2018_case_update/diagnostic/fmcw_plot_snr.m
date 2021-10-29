function fmcw_plot_snr(filename,burst)

% fmcw_plot_SNR(vdat) or fmcw_plot_snr(filename,burst)
%
% Plot amplitude, SNR and range error for a burst
%
% SNR and range error are plotted for a chirp and for the burst mean
% (i.e. using standard error of the chirp ensemble). Note the burst is
% split into common attenuator settings before processing.
%
% args:
% filename: filename (or else can be vdat data structure)
% burst: burst number in file

% Craig Stewart
% 2014/10/24

if nargin == 0
    [filename,~] = uigetfile({'*.dat','*.DAT'});
end
if isa(filename,'char')
    if nargin == 1
        burst = 1; % default burst if unspecified
        disp(['Using default burst number: ' int2str(burst)])
    end
    vdat = fmcw_load(filename,burst); % first burst
elseif isa(filename,'struct')
    % we've been passed data not a filename
    vdat = filename;
else
    error('argument 1 should be a filename or the vdat data structure')
end

% Loop through attenuator settings (1 plot per attenuator setting)
vdats = fmcw_burst_split_by_att(vdat);
for ii = 1:length(vdats)
    v = vdats(ii); % Just use the first settting for simplicity (for now)
    att = unique(v.chirpAtt);
    
    c = fmcw_range2(v,2,2000);
    S = fmcw_burst_stats(c.specRaw);
    % Range errors
    cre = v.lambdac*S.stdPhase/(4*pi); % chirp range error
    bre = v.lambdac*S.stePhase/(4*pi); % burst range error
    
    % Plot
    legendTxt = {'chirp','burst average'};
    
    figure
    ax(1) = subplottight(3,1,1);
    hc = plot(c.rangeCoarse,dB(abs(c.specRaw)),'b');
    hc = hc(1); % keep just one handle for legend
    hold on
    hb = plot(c.rangeCoarse,dB(abs(S.mean)),'r');
    set(gca,'XTickLabel',[])
    grid on
    ylabel('Ampltude (dB)')
    legend([hc hb],legendTxt)
    title(['File: ' v.filename ' burst: ' v.Burst ' Att: ' int2str(real(att)) '+' int2str(imag(att)) ' dB'])
    
    ax(2) = subplottight(3,1,2);
    plot(c.rangeCoarse,10*log10(S.stdSNR),'b')
    hold on
    plot(c.rangeCoarse,10*log10(S.steSNR),'r')
    set(gca,'XTickLabel',[])
    xlabel('Range (m)')
    grid on
    %set(gca,'XTickLabel',[],'YAxisLocation','right')
    ylabel('SNR')
    legend(legendTxt)
    
    ax(3) = subplottight(3,1,3);
    plot(c.rangeCoarse,cre,'b')
    hold on
    plot(c.rangeCoarse,bre,'r')
    set(ax(3),'YScale','log')
    grid on
    ylabel('range error (m)')
    legend(legendTxt)
    
    linkaxes(ax,'x')
end
