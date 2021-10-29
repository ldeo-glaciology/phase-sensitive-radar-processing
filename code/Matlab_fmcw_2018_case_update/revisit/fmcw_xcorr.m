function [iw,ic,cc,lags,pe,pse] = fmcw_xcorr(f,g,fi,maxlag,fe,ge,p)
    
% [rw,ic,cc] = fmcw_xcorr(f,g,r,fii,maxOffset)
%
% Cross-correlate portions of two complex vectors.
% This returns correlation coefficients of 1 for identical circshifted
% vectors irrespective of shift and window size (as long as the matching 
% segment hasn't been shifted off the end of g!).
%
% args:
% f = signal 1 (complex)
% g = signal 2 (complex)
% r = range vector for f and g
% fi = indices of f segment to match 
% maxlag = maximum bin offset between shots to search for. This can be a 
%          scalar (applied as +-) or a 2 element vector [min max] of offsets.
%
% optional args:
% fe = error on f
% ge = error on g
%
% output vectors:
% iw = effective range index of correlation (amplitude centre of mass)
% ic = incoherent correlation (scaled 0-1)
% cc = coherent correlation (complex) (scaled 0-1)
% lags = offsets used (units bin lags - not range units);
% pe = phase error of combined phase difference
% pse = phase standard error of combined phase difference
% 
% note: -angle(cc) is the average phase difference over depthRange
%          abs(cc) is the same as ac if the f-g phase difference is constant over depthRange

% Craig Stewart
% 2014/6/6

if nargin == 0
    % Create test data
    n = 500;
    f = rand(n,1) + 1i*rand(n,1);
    offset = 75;
    g = circshift(f,[offset 0]);
    phaseShift = 0.2345; % radians
    g = g.*exp(1i*phaseShift);
    fi = [250:270];
    maxlag = offset + 10;
    
    % Correlate
    [iw,ic,cc,lags] = fmcw_xcorr(f,g,fi,maxlag);
    [mic,ici] = max(ic); % get index (mci) of best amplitude correlation
    [mcc,cci] = max(cc); % get index (mci) of best complex correlation
    
    % Diplay results
    disp('Using test data:')
    disp(['offset = ' int2str(offset)])
    disp(['phase shift = ' num2str(phaseShift)])
    disp(' ')
    disp('Results:')
    disp('Incoherent:')
    disp(['Offset      = ' int2str(lags(ici))])
    disp(['Correlation = ' num2str(mic)])
    disp('Coherent:')
    disp(['Offset      = ' int2str(lags(cci))])
    disp(['Correlation = ' num2str(mcc)])
    disp(['Amplitude   = ' num2str(abs(mcc))])
    disp(['Phase shift = ' num2str(-angle(mcc))])
    %keyboard
    return
end

% Check inputs
f = f(:);
g = g(:);
fi = fi(:);

% Convert maxlag to lag vector lags
if numel(maxlag)==1 % then this is a maximum lag magnitude
    lags = -abs(maxlag):abs(maxlag); 
elseif numel(maxlag)==2 % then these are absoltue lag limits
    lags = maxlag(1):maxlag(2);
elseif numel(maxlag)>2
    error('maxlag should be a 1 or 2 element vector')
end

% Pad input vectors to cope with end effects (only necessary if we're
% trying to shift off the end of the vector)
n = max(abs(lags));
fi = fi+n;
zp = zeros(n,1);
f = [zp; f; zp];
g = [zp; g; zp];

if nargin>4
    doError = 1;
    fe = fe(:);
    fe = [zp; fe; zp];
    ge = ge(:);
    ge = [zp; ge; zp];
    fec = fe(fi);
else
    doError = 0;
end
fc = f(fi); % crop profile 1 to range
cff0 = sum(abs(fc).^2); % Compute autocorrelation at zero lag
[iw,ic,cc,pe,pse] = deal(zeros(size(lags)));
for ii = 1:length(lags) % loop through offsets
    gc = g(fi+lags(ii)); % shift and crop profile 2
    fg = fc.*conj(gc);
    cgg0 = sum(abs(gc).^2); % Compute autocorrelation at zero lag
    scale = sqrt(cff0*cgg0);
    cc(ii) = sum(fg)./scale; % complex correlation
    ic(ii) = sum(abs(fg))./scale; % amplitude only correlation
    iw(ii) = mean(abs(fg).*fi)/mean(abs(fg))-n; % effective bin of correlation (correcting for zero pad)
    
    % error
    if doError
        % Calculate phase error on this segment of f.g'*
        gec = ge(fi+lags(ii));
        pder = sqrt(fec.^2 + gec.^2); % phase difference standard deviation (fractional error of product)
        ab = abs(fc).*abs(gc); % length of each vector in the sum (at best lag) i.e. weights for errors
        sum_er = sqrt(sum((pder.*ab).^2)); % total combined (weighted) error over the xcorr
        sum_mag = abs(sum(fc.*conj(gc))); % magnitude of the xcorr
        pe(ii) = sum_er./sum_mag; % combined phase error over this window (radians)
        % Calculate the phase standard error taking into account padding
        % and weighting
        ne = sum(abs(fg))./max(abs(fg)); % effective number of samples with weighting
        ne = ne/p; % account for padfactor
        if ne<0
            keyboard
        end
        pse(ii) = pe(ii)./sqrt(ne); % don't use ne-1, or we get complex when ne<1 !
        
%         figure
%         plot(fc,'r')
%         hold on
%         plot(gc,'b')
%         
%         figure
%         plot(abs(fc),'r')
%         hold on
%         plot(abs(gc),'b')
%         keyboard
    end
end

if ~doError
    [pe,pse] = deal(nan*cc);
end
