function [vdat,n] = fmcw_cull_bad(vdat,noisePowerLimit,doPlot)

% Cull grossly contaminated chirps from a FMCW radar burst using chirp
% power level relative to cleanest chirp in burst.
%
% args: vdat = fmcw data structure
% noisePowerLimit = fractional power difference from trend
%
% Craig Stewart
% 2014/5/30

if nargin < 2
    noisePowerLimit = 0.005;
end
if nargin < 3
    doPlot = 0;
end
doDebug = 0;

% Check if there are more than one attenuator settings
if length(unique(vdat.chirpAtt)) > 1
    error('Multiple attenuator settings present in burst - subset into attenuator settings before culling noisy chirps')
end

% Check that there are enough chirps to try...
minChirps = 5;
nChirps = size(vdat.vif,1);
if nChirps<minChirps
    disp([mfilename ' cancelled - Only ' int2str(nChirps) ' chirps in burst.'])
    vdat.processing = [vdat.processing {[mfilename ': cancelled - to few chirps to cull']}];
    n = 0;
    return
end

% Find bad shots
v = vdat.vif;
mvc = repmat(mean(v,2),1,size(v,2));
p = rms(v-mvc,2); % rms power of each chirp
cn = [1:size(v,1)]';

% Two methods
% min: detects chirps more than noisePowerLimit above minimum power chirp
% fit: detects only anomalies more than noisePowerLimit from the trend
% New method using robustfit
fitType = 'min'; %'min','lin','quad','cub','exp';
switch fitType
    case 'min'
        noisey = p > min(p)*(1+noisePowerLimit); % this techniques fails for case of 1 quiet chirp
        bpred = repmat(min(p),size(p)); % for consistency with case:fit
        resid = p - bpred;
    case 'quad'
        [b,stats] = robustfit([cn cn.^2],p); % robust quadratic fit to power curve
        bpred = b(1) + b(2).*cn + b(3).*cn.^2; % prediction from fit
        resid = stats.resid;
    case 'cub'
        [b,stats] = robustfit([cn cn.^2 cn.^3],p); % robust quadratic fit to power curve
        bpred = b(1) + b(2).*cn + b(3).*cn.^2 + b(4).*cn.^3; % prediction from fit
        resid = stats.resid;
    case 'exp' % doesn't work...
        ft = fittype('exp1');
        f = fit(cn,p,ft);
        bpred = f.a*exp(f.b*cn); % predicted values
        resid = p - bpred;
        %plot(ef,x,y)
end
noisey = abs(resid) > median(p)*noisePowerLimit; % noisey from residual to fit
nChirps = size(v,2);
n = sum(noisey);


if n > 0
    % Warn if we're loosing lots of data
    if  n >= nChirps/2
        disp(['Warning: found ' int2str(n) ' bad chirps out of ' int2str(nChirps)])
        doPlot = 1;
        doDebug = 1; % after plotting
    end
    % Record processing then crop
    vdat.processing = [vdat.processing {[mfilename ': removed ' int2str(n) ' contaminated chirps: ' mat2str(vdat.chirpNum(noisey))]}];
    vdat = fmcw_burst_subset(vdat,find(~noisey));
else
    vdat.processing = [vdat.processing {[mfilename ': no bad chirps found']}];
end

if doPlot
    figure
    subplot(1,2,1)
    plot(cn,p,'b.')
    hold on
    plot(cn,bpred,'g')
    plot(cn(noisey),p(noisey),'ro')
    xlabel('chirp')
    ylabel('volts rms')
    title('chirp power')
    
    subplot(1,2,2)
    hist(resid)
    xlabel(['power residual from ' fitType ' fit'])
    ylabel('num occurance')
end
if doDebug
    keyboard
end