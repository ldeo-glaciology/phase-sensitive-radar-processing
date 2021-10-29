function bed = fmcw_bed_info(vdat,depthRange,doPlot)

% bpg = fmcw_bed_info(vdat,depthRange,doPlot)
%
% Calculate properties of basal return, including range, amplitude and phase gradient
%
% arg in:
% vdat: radar data structure 
% depthRange: search range for bed (m) (e.g. [240 280])
% doPlot: plot flag (0/1)
%
% arg out:
% returns structure bed with fields:
% .range: fine range (m)
% .amp: amplitude
% .phase: (rad)
% .phaseGrad: phase gradient
% .phaseDevVar: variance phase deviation from linear
 
% Craig Stewart
% 2014/10/21

if nargin == 0
    % Load file
    vdat = fmcw_load;
end
if nargin < 2
    depthRange= [0 2000];
end
if nargin < 3
    doPlot = 1;
end

% Phase processing config
p = 4; % pad val (interpolation)
maxRange = 2000;
win = @blackman;

% Split into seperate attenuator settings
vdats = fmcw_burst_split_by_att(vdat);
vdat = vdats(1); % keeping only the first attenuator setting (so we don't try to average across multiple att sets

%% Individual chirp properties
c = fmcw_range2(vdat,p,maxRange,win);


%% Burst mean properties
S = fmcw_burst_stats(c.specRaw);

% Average all chirps in burst
%vdat = fmcw_burst_mean(vdat);
% Phase process for range
% [rangeCoarse,rangeFine,~,specRaw] = fmcw_range(vdat,p,maxRange,win);

% Find bed
bedmethod = 'maxAmp';
bn = fmcw_findbed(c.rangeCoarse,abs(S.mean),depthRange,bedmethod); % index of bed

% Bed phase evolution throughout burst
ba = abs(c.specRaw(:,bn)); % bed amp (each chirp)
bp = angle(c.specRaw(:,bn)); % phase at bed (each chirp)
br = vdat.lambdac*(bp-bp(1))/(4*pi); % above converted to range
cn = transpose(1:length(ba));

% Calculate best fit to phase evolution throughout burst
fitType = 'exp'; %'exp','asym','poly5'
switch fitType
    case 'exp'
        ft=fittype('exp1');
        f = fit(cn,bp,ft);
        bp_pred = f.a*exp(f.b*cn); % predicted values
    case 'asym' % asymptotic
        ft = fittype('a+(b-a)*exp(c*x)');
        f = fit(cn,bp,ft);
        bp_pred = f.a*exp(f.b*cn); % predicted values
    case 'poly5'
        [b,stats] = robustfit([cn cn.^2 cn.^3 cn.^4],bp); % robust quadratic fit to power curve
        bp_pred = b(1) + b(2).*cn + b(3).*cn.^2 + b(4).*cn.^3 + b(5).*cn.^4; % prediction from fit
        %noisey = abs(stats.resid) > median(p)*noisePowerLimit; % noisey from residual to fit
end
bp_resid = bp-bp_pred;
br_pred = vdat.lambdac*(bp_pred-bp(1))/(4*pi); % above converted to range
br_resid = vdat.lambdac*(bp_resid)/(4*pi); % above converted to range
disp(['bed range std         : ' num2str(std(br)) ' m'])
disp(['bed range residual std: ' num2str(std(br_resid)) ' m'])

% Basic bed characteristics
bed.range = c.rangeCoarse(bn) + c.rangeFine(bn);
bed.amp = abs(S.mean(bn));
bed.ampdB = dB(abs(S.mean(bn)));
bed.phase = angle(S.mean(bn));

% Calculate phase grad
bedInds = [bn-1*p:bn+0*p]; % -1 range bin, i.e. from - 0.41m to bed peak
omega = angle(S.mean(bedInds(end)).*conj(S.mean(bedInds(1))))./(bedInds(end)-bedInds(1)); % phase gradient per bin
dr = diff(c.rangeCoarse(1:2));
bed.phaseGrad = omega./dr; % phase grad per metre (rad/m)

% Calculate phase deviation variance (from linear)
% remove phase gradient by multiplying by negative phasor of gradient
specRawFlat = S.mean.* exp(-1i*omega*(1:length(S.mean))); 
phaseDev = angle(specRawFlat.*conj(specRawFlat(bn))); % Relative to bed phase (to avoid wrapping)
bedInds = [bn-1*p:bn+1*p]; % +- 1 range bin, i.e. +- 0.41m (inside blackman main lobe)
bed.phaseDevVar = var(phaseDev(bedInds));

% Plot
if doPlot
    % evolution throughout burst
    figure
    ax(1) = subplot(4,1,1);
    plot(cn,ba)
    xlabel('chirp no')
    ylabel('bed amp')
    title(vdat.filename,'interpreter','none')
    
    ax(2) = subplot(4,1,2);
    plot(cn,bp_pred,'r');
    hold on
    plot(cn,bp,'b.');
    legend([fitType ' fit'],'data')
    xlabel('chirp no')
    ylabel('bed phase')
    
    ax(3) = subplot(4,1,3);
    plot(cn,1000*br_pred,'r');
    hold on
    plot(cn,1000*br,'b.')
    legend([fitType ' fit'],'data')
    xlabel('chirp no')
    ylabel('bed fine range (mm)')
    
    ax(4) = subplot(4,1,4);
    plot(cn,1000*br_resid,'r');
    %hold on
    %plot(cn,1000*br)    
    xlabel('chirp no')
    ylabel('bed fine range residual (mm)')
    
    linkaxes(ax,'x')
    clear ax
    
    % Mean of burst: bed phase grad
    figure
    ax(1) = subplot(2,1,1);
    plot(c.rangeCoarse,dB(abs(S.mean)),'b');
    hold on
    plot(c.rangeCoarse(bedInds),dB(abs(S.mean(bedInds))),'b.')
    plot(c.rangeCoarse(bn),dB(abs(S.mean(bn))),'bo')
    plot(c.rangeCoarse,dB(abs(specRawFlat)),'r')
    plot(c.rangeCoarse(bedInds),dB(abs(specRawFlat(bedInds))),'r.')
    
    ax(2) = subplot(2,1,2);
    plot(c.rangeCoarse,angle(S.mean),'b');
    hold on
    plot(c.rangeCoarse(bedInds),angle(S.mean(bedInds)),'b.');
    plot(c.rangeCoarse(bn),angle(S.mean(bn)),'bo');
    plot(c.rangeCoarse,angle(specRawFlat),'r')
    plot(c.rangeCoarse(bedInds),angle(specRawFlat(bedInds)),'r.')
    text(c.rangeCoarse(bedInds(end)),angle(S.mean(bedInds(end))),[num2str(bed.phaseGrad,2) ' '],'color','b','HorizontalAl','right','VerticalAl','bot')
    
    linkaxes(ax,'x')
    set(ax,'xlim',[bed.range-5 bed.range+5])
end

%bed
%keyboard
