function fmcw_ts_range3(cfgfile)

% fmcw_ts_range3(cfgfile)
% 
% modified from fmcw_ts_range2 - 
% to start with differencing each chirp internal with the bed - to get rid
% of any effects from chirp to chirp absolute phase differences
%
% Process sequence of FMCW shots to calculate timeseries of bed-internal phase difference. 
%
% intput arg: cfgfile = text file of settings and input data filenames
% this includes manually selected internal and bed depth as well and file
% names and burst numbers to process. If this is absent the user is
% prompted to make one. 
% Output is saved in a matfile with the same name as 
% the cgffile. Data is plotted using fmcw_timeseries_plot.
%
% Processing method consisists of 
% (1) define a cfg file which lists the files to process, the chosen
% internal and bed depths, vertical strain rate for the site and metadata. 
% (2) sequentially load and phase process each burst
% (3) calculate mean phase of internals (angle of vector mean of internals)
% (4) calculate phase of bed (angle of bed)
% (5) calculate internal bed phase diff
% (6) calculate change in (internal bed phase diff) over adjacent bursts
% (i.e. d(phi)/d(time)

% (4) Data saved in .mat file.
%
% Data should then be cleaned with fmcw_ts_clean, this is kept modular to
% avoid having to run this (processor heavy) script more than once.

%
% Craig Stewart


% Plot options
doPlotTest = 0; % plotting intermediate data for debugging

% Select cfg file if not entered
if nargin == 0
    cfgfile = uigetfile('*.cfg','Select a cfg file - or cancel to create new');
end
if cfgfile==0
    % Create cfg file
    cfgfile = fmcw_batchlist_write2;
end

% Read cfg file
[cfg,file] = fmcw_batchlist_read2(cfgfile);

% Processing options
cfg.bedSeachMethod = 'maxAmp';
cfg.intMode = 'peaks'; % 'all','peaks','uninterp'
cfg.intMinPeakSep = 5; % minimmum separation between internal peaks
% Phase processing settings
cfg.padfactor = 1;

% Count bursts so we can preallocate matricies before processing
nbursts = zeros(1,length(file));
disp('Checking for number of bursts to process')
for fn = 1:length(file)
    filename = char(file(fn).name);
    disp(['file ' int2str(fn) '/' int2str(length(file)) ':' filename])
    if isnan(file(fn).burst) % use all bursts from file
        nbursts(fn) = fmcw_nbursts(filename);
    end
end
N = sum(nbursts);
[time] = deal(ones(N,1));
[I,Ie,Ic,Imd,Ipv,bedInd] = deal(ones(N-1,1)); % differences

%% Loop through all files and bursts averaging all chirps within the burst
bi = 0; % burst index
ci = 0; % comparison index (always bi-1)
for fn = 1:length(file)
    filename = char(file(fn).name);
    for burst = 1:nbursts(fn)
        bi = bi+1; % master index of burst
        if bi == 1
            % Load burst 1
            G = getSpec(filename,burst,cfg); % load spectrum
            time(bi) = G.time;
            
            % Find internals
            switch cfg.intMode
                case 'all'
                    isInt = G.rangeCoarse>=min(cfg.intRange) & G.rangeCoarse<=max(cfg.intRange); % is in internal range
                case 'uninterp'
                    isInt = G.rangeCoarse>=min(cfg.intRange) & G.rangeCoarse<=max(cfg.intRange)...
                        & ~mod([1:length(G.rangeCoarse)],cfg.padfactor); % is in internal range and is not an interpolated point
                case 'peaks'
                    % Now only select peaks to get only best SNR
                    dr = mean(diff(G.rangeCoarse));
                    mpdb = cfg.intMinPeakSep/dr; % min peak separation distance (in bins)
                    [~,pkind] = findpeaks(abs(G.mean),'MinPeakDistance',mpdb); % amplitude peaks separated by 4m
                    isPeak = zeros(size(G.mean));
                    isPeak(pkind) = 1;
                    isInt = G.rangeCoarse>=min(cfg.intRange) & G.rangeCoarse<=max(cfg.intRange) & isPeak; % is in internal range
            end
            %intDepth = G.rangeCoarse(isInt);
            
            % Preallocate large matricies
            [FG,FGSEP] = deal(ones(N-1,numel(G.rangeCoarse))); % differences
            [FGI] = deal(ones(N-1,sum(isInt))); % differences internals only
            
        else
            % Compare shots
            ci = ci + 1; % comparison index;
            
            % Find index of bed (at each time)
            bedInd(ci) = fmcw_findbed(G.rangeCoarse,abs(G.mean),cfg.bedSearchRange,cfg.bedSeachMethod); % bed ind
            if bedInd(ci) == G.bedInd
                % We already loaded the file last time with the correct bed
                F = G; 
            else
                % Bed has moved between bursts so reload burst using new bed location
                F = getSpec(G.filename,G.burst,cfg,bedInd(ci)); % load spectrum
            end
            G = getSpec(filename,burst,cfg,bedInd(ci)); % load spectrum
            time(bi) = G.time;
            disp([int2str(ci) ') F: ' F.filename ' b:' int2str(F.burst) ' t:' datestr(F.time)]); % and ' G.filename '(burst:' int2str(G.burst) ')'])
            
            % Compare bursts
            fg = conj(F.mean).*G.mean;
            fgsep = sqrt(F.stePhase.^2 + G.stePhase.^2); % standard error (phase) of fg (assuming uncorrelated noise)
            fgse = sqrt(2)*fgsep.*abs(fg); % standard error of fg
            % Keep a few figures for all internals
            FG(ci,:) = fg;
            FGI(ci,:) = fg(isInt);
            FGSEP(ci,:) = fgsep;
            %FGSEPI(ci,:) = fgsep(isInt);
            
            % Internals
            I(ci) = sum(fg(isInt)); % sum fg over internals
            Ic(ci) = I(ci)./(sqrt(sum(abs(F.mean(isInt)).^2)).*sqrt(sum(abs(G.mean(isInt)).^2))); % coherence
            % Calculate phase variance from coherence
            gamma = Ic(ci);
            In = sum(isInt); % number of internals
            Ipv(ci) = (1/sqrt(2*In))*sqrt(1-abs(gamma)^2)./abs(gamma);% http://www.esa.int/esapub/tm/tm19/TM-19_ptC.pdf eq. 1.18
            Ie(ci) = norm(fgse(isInt)); % total error of I
            Imd(ci) = mean(F.rangeCoarse.*abs(fg))./mean(abs(fg)); % effective mean internal depth (amplitude weighted)
            
            % Bed (for now just a single bin)
            %B(ci) = fg(bedInd(ci)); % fg at bed
            %Bpe(ci) = fgsep(bedInd(ci)); % fg phase standard error at bed
            %B(ci) = sum(fg(isBed))
            %Be(ci) = norm(fgse(isBed)); % total error of I -if we had a
            %bed range rather than a point
            
            % Calculate noise floor
            Fnm(ci) = mean(F.ste(isInt)); % mean noise
            %disp(['F noise mean: ' num2str(dB(Fnm(ci)))])
            Fpnp = Fnm(ci)./(sqrt(2)*abs(F.mean(isInt))); % phase noise prediction assuming white noise in RF
            Gnm(ci) = mean(G.ste(isInt)); % mean noise
            Gpnp = Gnm(ci)./(sqrt(2)*abs(G.mean(isInt))); % phase noise prediction assuming white noise in RF
            % note there appears to be a minimum phase noise of 1e-3 rad...
            % wny - this should eb much lower?
                        
            % Stop to investigate noise
            if doPlotTest || ci==N-1
                title_txt = [int2str(ci) ') F: ' F.filename ' b:' int2str(F.burst) ' t:' datestr(F.time)];
                %disp(['F.time: ' datestr(F.time)])
                %disp(['G.time: ' datestr(G.time)])
                
                dr = mean(diff(F.rangeCoarse));
                mpd = 4; % minimum peak distance from other peak (m)
                mpdb = mpd/dr; % (in bins)
                [pk,pkind] = findpeaks(sqrt(abs(fg)),'MinPeakDistance',mpdb); % amplitude peaks separated by 4m
                
                % Amplitude profile and phase diff at internals
                figure
                ax(1) = subplottight(2,1,1);
                plot(F.rangeCoarse,dB(abs(F.mean)),'r')
                hold on
                plot(G.rangeCoarse,dB(abs(G.mean)),'b')
                %plot(F.rangeCoarse,dB(abs(fg)),'g')
                plot(F.rangeCoarse(pkind),dB(pk),'ro')
                legend('f','g')
                set(gca,'XTickLabel',[])
                title(title_txt)
                
                ax(2) = subplottight(2,1,2);
                %plot(F.rangeCoarse(isInt),angle(fg(isInt)),'k.')
                hold on
                %[h1,h2] = erbar(F.rangeCoarse(pkind),angle(fg(pkind)),fgsep(pkind),-fgsep(pkind),'r','b');
                [h1,h2] = erbar(F.rangeCoarse(isInt),angle(fg(isInt)),fgsep(isInt),-fgsep(isInt),'r','b');
                ylim([-0.3 0.3])
                ylabel('phase difference')
                xlabel('range (m)')
                box on
                
                linkaxes(ax,'x')
                
                % scatter
                figure
                ax(1) = subplottight(2,1,1);
                %scatter(dB(abs(fg(pkind))),angle(fg(pkind)),20,F.rangeCoarse(pkind),'filled')
                %scatter(dB(abs(fg(isInt))),fgsep(isInt),15,angle(fg(isInt)),'filled')
                plot(dB(abs(F.mean(isInt))),F.stePhase(isInt),'b.','markersize',15);
                hold on
                plot(dB(abs(G.mean(isInt))),G.stePhase(isInt),'r.','markersize',15);
                %plot(dB(abs(fg(pkind))),fgsep(pkind),'ko') % show peaks
                colormap jet
                grid on
                %hold
                ylabel('phase standard error')
                %ch = colorbar('East');
                %ylabel(ch,'phase difference (rad)')
                set(gca,'XTickLabel',[])
                %ch.Label.string = 'amplitude (dB)';
                
                
                
                % overlay prediction of noise from this noise floor
                hold on
                plot(dB(fliplr(sort(abs(F.mean(isInt))))),sort(Fpnp),'b'); % predicted phase noise for F
                plot(dB(fliplr(sort(abs(G.mean(isInt))))),sort(Gpnp),'r'); % predicted phase noise for F
                %plot(dB(abs(G.mean(isInt))),Gpnp,'r.'); % predicted phase noise for G
                xlabel('amp')
                title(title_txt)
                
                % Plot phase difference
                ax(2) = subplottight(2,1,2);
                erbar(dB(sqrt(abs(fg(isInt)))),angle(fg(isInt)),fgsep(isInt),-fgsep(isInt),'r','b');
                xlabel('amp')
                ylabel('phase difference')
                linkaxes(ax,'x')
                
%                 figure
%                 hist(dB(noise))
%                 xlabel('noise level')
%                 ylabel('num. bins')
%                 title(title_txt)
                
                %keyboard
            end
        end
    end
end

daysPerYear = 365.25;

% Strain
Br = transpose(F.rangeCoarse(bedInd));
IBr = Br-Imd;
drs = IBr.*diff(time)*(cfg.vsr/daysPerYear);

% Range change
%IB = conj(I).*B; % bed-(mean)internal difference
dr = F.lambdac*angle(I)./(4*pi);
mr = -daysPerYear*(dr-drs)./diff(time); % melr rate (corrected for strain)

% Error (of range change)
%Ipe = Ie./(sqrt(2)*abs(I)); % phase error of I
Ipe = Ipv; % use observed variation in phase from reflectors - rather than 
% combination of phase errors at reflectors to estimate actual phase error
% done through the ESA equation for phase variance based on coherence.

%Bpe = Be./(sqrt(2)*abs(B)); % phase error of B
%IBpe = sqrt(Ipe.^2 + Bpe.^2); % bed-(mean)internal phase standard error
dre = F.lambdac*Ipe./(4*pi);
mre = -365.25*dre./diff(time);

timeC = (time(1:end-1) + time(2:end))./2; % time at centre between measurements
r = cumsum(dr); % range

% % Smoothing...
% % Using singular spectrum analysis
% mrdt = detrend(mr);
% mrt = mr-mrdt; % trend
% n = round(length(mr)/2);
% mrssa = ssa(mrdt,n,'noplot','silent');
% n2 = round(0.8*n);
% mrsdt = ssarecon(mrssa,n2:n,'noplot','silent');
% mrs_ssa = mrsdt + mrt; % retrended to reconstruct smootghed meltrate time series

% Running mean
mrs3 = runningmean(mr,3);
mrs5 = runningmean(mr,5);


%% Plot

% Internal and bed phase
figure
ax(1) = subplottight(3,1,1);
h1 = plot(timeC,angle(I),'col','b');
hold on
%[hb,hp(1)] = erbar(timeC,angle(I),-Ipe,Ipe,'b','k');
%h2 = plot(timeC,angle(B),'col','r');
%[hb,hp(2)] = erbar(timeC,angle(B),-Bpe,Bpe,'r','k');
h = plot(timeC,angle(I),'col','k');
[hb,hp] = erbar(timeC,angle(I),-Ipe,Ipe,'k','k');
grid on
ylabel('inter-shot phase difference (rad)')
legend(h,'bed-internal')
%legend('internal','bed','bed-internal')
xtimelab('tg')
% ax(3) = subplot(4,1,3);
% plot(timeC,dr);
% ylabel('range change (m)')
% xtimelab
ax(2) = subplottight(3,1,2);
h1 = plot(timeC,mr,'col',[0.6 0.6 0.6]);
hold on
erbar(timeC,mr,mre,-mre,'r','b'); % [hb,hp] = 
% Plot smoothed melt rates
hold on
%plot(timeC,mrs_ssa,'k')
%h2 = plot(timeC,mrs3,'g');
%h3 = plot(timeC,mrs5,'c');
%legend([h1 h2 h3],'raw','3-point mean','5-point mean')
grid on
ylabel('meltrate (m/year)')
xtimelab('tg')

ax(3) = subplottight(3,1,3);
h1 = plot(timeC,mr,'.','col',[0.6 0.6 0.6]);
hold on
h2 = plot(timeC,mrs3,'g');
h3 = plot(timeC,mrs5,'c');
legend([h1 h2 h3],'raw','3-point mean','5-point mean')
grid on
ylabel('meltrate (m/year)')
xtimelab

linkaxes(ax,'x')



% 
% % Errors
% % Time series of phase errors and bed index
% figure
% ax(1) = subplot(2,1,1);
% plot(timeC,Ipe,'b')
% %legend('internals','bed')
% xtimelab('tg')
% ylabel('phase error')
% 
% ax(2) = subplot(2,1,2);
% plot(timeC,bedInd,'r')
% ylabel('bed index')
% xtimelab
% linkaxes(ax,'x')

% % Time series of range and range error
% figure
% ax(1) = subplot(2,1,1);
% plot(timeC,r,'-.')
% ylabel('Cumulative range change (m)')
% 
% ax(2) = subplot(2,1,2);
% plot(timeC,dre,'r')
% ylabel('range error')
% linkaxes(ax,'x')

% Timeseries of coherence
figure
plot(timeC,abs(Ic),'.')
ylabel('coherence')
xtimelab


% % Timeseries of phase-difference for each internal
% figure
% % add errorbars
% numInt = sum(isInt);
% for ii = 1:numInt
%     hold on
%     h(ii) = plot(timeC,angle(FGI(:,ii)));
%     %col = get(h(ii),'col');
%     legtxt(ii) = {[num2str(intDepth(ii)) ' m']};
%     %erbar(timeC,angle(FGI(:,ii)),-FGSEPI(:,ii),FGSEPI(:,ii),'col',col);
% end
% xtimelab
% ylabel('phase difference')
% legend(h,legtxt)

% Noise floor timeseries
figure
plot(time(1:end-1),dB(Fnm),'b')
hold on
%plot(timeC,dB(Gnm),'r')
xtimelab
ylabel('noise floor')

keyboard

% cfg.lambdac = vdat.lambdac;

% %% Save results to file
% [~,cfgfilename,~] = fileparts(cfgfile);
% outfilename = ['fmcw_timeseries_' cfgfilename];
% disp(['Saving data to: ' outfilename])
% save(outfilename,'time','temp1','temp2','SR','rangec','range','rangeSD','rangeSE','intIndFx','bedIndFx','intIndPk','bedIndPk','cfg')

% Calculate melt and plot
%fmcw_ts_clean(outfilename)


function X = getSpec(filename,burst,cfg,bedInd)

[vdat] = fmcw_load(filename,burst); % load data
% Now keep only chirps with the first attenuator setting
vdats = fmcw_burst_split_by_att(vdat);
vdat = vdats(1); % Just use the first settting for simplicity (for now)

% Cleaning
doClean = 1;
if doClean
    %frange = [205e6 395e6];
    %vdat = fmcw_cull_freq(vdat,frange);
    
    % drop first 5 chirps
    %nCull = 5; 
    %vdat = fmcw_burst_subset(vdat,[nCull+1:vdat.ChirpsInBurst]);
    
    % Kill remaining outliers
    noisePowerLimit = 0.01;
    doPlot = 0;
    [vdat,~] = fmcw_cull_bad(vdat,noisePowerLimit,doPlot);
    nNoisey = 1;
    vdat = fmcw_cull_noisey(vdat,nNoisey);
end

% Process all chirps to get SNR etc
c = fmcw_range2(vdat,cfg.padfactor,cfg.bedSearchRange(2));
if nargin < 4
    % Find bed
    bedInd = fmcw_findbed(c.rangeCoarse,abs(mean(c.specRaw,1)),cfg.bedSearchRange,cfg.bedSeachMethod); % bed ind
end
% Subtract the bed phase from all depths (for each chirp)
c.specRaw = conj(c.specRaw).*repmat(c.specRaw(:,bedInd),1,size(c.specRaw,2))./abs(mean(c.specRaw(:,bedInd)));
% Take mean and get stats
S = fmcw_burst_stats(c.specRaw);

% Mean of chirps in burst
%vdatm = fmcw_burst_mean(vdat);
% SNR improvement from averaging = SNRgain = 20*log10(sqrt(vdat.ChirpsInBurst))

% Phase processing
%[X.rangeCoarse,~,~,X.specRaw] = fmcw_range(vdatm,cfg.padfactor,cfg.maxRange); % phase process data (from mean shot)
X.rangeCoarse = c.rangeCoarse;
X.mean = S.mean;
X.ste = S.ste;
X.stdPhase = S.stdPhase;
X.stePhase = S.stePhase;
X.time = vdat.TimeStamp;
%X.temp1 = vdatm.Temperature_1;
%X.temp2 = vdatm.Temperature_2;
X.lambdac = vdat.lambdac;
X.filename = filename;
X.burst = burst;
X.bedInd = bedInd;
X.burst = burst;