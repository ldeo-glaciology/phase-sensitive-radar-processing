function fmcw_synth_test(test)

% fmcw_synth_test(test)
%
% Investigate various scenarios by synthesizing different radar profiles
% and plotting or processsing
%
% test options:
% 'range'
% 'close reflectors'
% 'layering'
% 'noise'
% 'vsr'
% 'melt'
% 'bed roughness'
% 'system gain'
%
% Craig Stewart
% 2014/6/1

if nargin ==0
    testOptions = {'range','close reflectors','layering','noise','clipping','vsr','depth shift','deconv','melt','bed roughness','system gain'};
    disp(' ')
    disp('Test scenarios')
    for ii = 1:length(testOptions)
        disp([int2str(ii) ': ' testOptions{ii}])
    end
    ans = input('Enter a test scenario number: ');
    test = testOptions{ans};
end

switch test
    case 'range'
        % Basic test of range detection on single peak
        
        % Test generate profile with one reflector
        %R1 = 3000.035;
        %R1 = 3000.266;
        %R1 = 8000.266; % max range 8600m
        %R1 = 30.05;
        %R1 = 29.965;
        R1 = 300;
        %R1 = 1800.12;
        %R1 = 1800;
        %R1 = 1800.032 - 0.21098;
        A1 = -20;
        N = nan; % no added noise
        S = 0.8; % scale to 0.8 of dynamic range
        %adcBitDepth = inf; % no digitisation
        %adcBitDepth = 12; %
        adcBitDepth = 16; % standard adc
        vdat = fmcw_synth(R1,A1,N,S,adcBitDepth);
        save synth_chirp_onePoint.mat vdat
        
        datamode = 'synth'; % 'field' 'synth'
        switch datamode
            case 'synth'
                vdat1 = fmcw_load('synth_chirp_onePoint.mat'); % this addes necessary metadata
            case 'field'
                %vdat1 = fmcw_load('a01_2013-01-16-0524.DAT'); %- usereal data
                [filename, ~] = uigetfile(['*.dat;*.DAT;*.000;*.mat'],'Choose radar file to plot');
                if isa(filename,'double') % no files chosen
                    return
                end
                disp(['Loading: ' filename])
                vdat1 = fmcw_load(filename); %- usereal data
                vdat1 = fmcw_burst_mean(vdat1);
        end
        p = 2;
        maxrange = 400;
        
        %% process normally
        
        %[rc1nopad,rf1nopad,sr1nopad,s1nopad] = fmcw_range(vdat1,1,maxrange);
        [rc1,rf1,sr1,s1] = fmcw_range(vdat1,p,maxrange);
        r1 = rc1 + rf1; % total range
        
        % Plot
        %fmcw_plot('tap','fl',{'synth_chirp_onePoint.mat'},'mr',maxrange)
        %overlay_reflectors(R1,A1,'b.')
        
        % Check location of bed return
        [~,ii] = max(dB(abs(s1))); % locate maximum strength reflector
        %[~,iinopad] = max(dB(abs(s1nopad))); % locate maximum strength reflector
        disp('Reflector range estimate: ')
        disp(['coarse: ' num2str(rc1(ii),'%12.10f') ' m'])
        disp(['fine  : ' num2str(rf1(ii),'%12.10f') ' m'])
        disp(['total : ' num2str(r1(ii),'%12.10f') ' m'])
        disp(['error : ' num2str(r1(ii)-R1,'%12.10f') ' m'])
        %keyboard
        
        % show some nearby reflectors
        figure
        BF = [1 2 4 8];
        ax(1) = subplottight(length(BF)+1,1,1);
        %range = 1.8;
        range = 4;
        gi = rc1>=(rc1(ii)-range) & rc1<=(rc1(ii)+range);
        %ginopad = find(rc1nopad>=(rc1nopad(iinopad)-range) & rc1nopad<=(rc1nopad(iinopad)+range));
        plot(rc1(gi),dB(abs(s1(gi))),'b')
        hold on
        rch = plot(rc1(gi),dB(abs(s1(gi))),'b.','markersize',12);
        rfh = quiver(rc1(gi),dB(abs(s1(gi))),rf1(gi),zeros(size(rf1(gi))),0,'ShowArrowHead','off','color','m');
        rtoth = plot(rc1(gi)+rf1(gi),dB(abs(s1(gi))),'r.','markersize',12);
        %plot(rc1nopad(ginopad),dB(abs(s1nopad(ginopad))),'ro')
        rchosen = plot(rc1(ii)+rf1(ii),dB(abs(s1(ii))),'ro','markersize',10);
        %         for ii = 1:length(s1nopad(ginopad))
        %             bch = plot(repmat(rc1nopad(ginopad(ii)),1,2),[-120 0],'k');
        %         end
        % Set grid to range bins
        %set(gca,'xtick',rc1(gi))
        halfwaves = R1 + (vdat1.lambdac/2)*[min(round((rc1(gi)-R1)./(vdat1.lambdac/2))):max(round((rc1(gi)-R1)./(vdat1.lambdac/2)))];
        %set(gca,'xtick',halfwaves(1:2:end))
        %         for ii = 1:length(halfwaves)
        %             waveh = plot(repmat(halfwaves(ii),1,2),[-120 0],'col',[0.6 0.6 0.6]);
        %         end
        if strcmp(datamode,'synth')
            rtrueh = plot([R1 R1],[-120 0],'g');
        end
        %legend([rtrueh rch rfh rtoth bch waveh rchosen],{'true range','coarse range (range at bin centres)','fine range','total range estimate','unpadded bin centres','integer half-wavelengths from reflector','chosen range'})
        legend([rch rfh rtoth rchosen],{'coarse range','fine range','total range estimate','chosen range'})
        set(gcf,'tag','fmcw_range_schematic')
        ylabel('Amplitude (dB)')
        xlabel('range (m)')
        ylims = [-100 0];
        set(gca,'xticklabel',[])
        ylim(ylims)
        
        % subplanel label
        x1 = 296.1;
        y1 = -15;
        fonts = 12;
        text(x1,y1,'(a)','fontsize',fonts)
        
        %% Process differentially
        %cnp = fmcw_drange(vdat1,10,maxrange); % no padding
        %rc1nopad = cnp.rangeCoarse;
        %rf1nopad = cnp.rangeFine;
        %sr1nopad = cnp.specCor;
        %s1nopad  = cnp.specRaw;
        
        splab = {'(b)','(c)','(d)','(e)'};
        for jj = 1:length(BF)
            Bf = BF(jj);
            c = fmcw_drange(vdat1,p,maxrange,@blackman,Bf); % padded
            rc1 = c.rangeCoarse;
            rf1 = c.rangeFine;
            sr1 = c.specCor;
            s1  = c.specRaw;
            r1 = rc1 + rf1; % total range
            
            % Plot
            %fmcw_plot('tap','fl',{'synth_chirp_onePoint.mat'},'mr',maxrange)
            %overlay_reflectors(R1,A1,'b.')
            
            % Check location of bed return
            [~,ii] = max(dB(abs(s1))); % locate maximum strength reflector
            %[~,iinopad] = max(dB(abs(s1nopad))); % locate maximum strength reflector
            disp('Reflector range estimate: ')
            disp(['coarse: ' num2str(rc1(ii),'%12.10f') ' m'])
            disp(['fine  : ' num2str(rf1(ii),'%12.10f') ' m'])
            disp(['total : ' num2str(r1(ii),'%12.10f') ' m'])
            disp(['error : ' num2str(r1(ii)-R1,'%12.10f') ' m'])
            %keyboard
            
            % plot
            ax(2) = subplottight(length(BF)+1,1,jj+1);
            gi = rc1>=(rc1(ii)-range) & rc1<=(rc1(ii)+range);
            %ginopad = find(rc1nopad>=(rc1nopad(iinopad)-range) & rc1nopad<=(rc1nopad(iinopad)+range));
            plot(rc1(gi),dB(sqrt(abs(s1(gi)))),'b')
            hold on
            rch = plot(rc1(gi),dB(sqrt(abs(s1(gi)))),'b.','markersize',12);
            rfh = quiver(rc1(gi),dB(sqrt(abs(s1(gi)))),rf1(gi),zeros(size(rf1(gi))),0,'ShowArrowHead','off','color','m');
            rtoth = plot(rc1(gi)+rf1(gi),dB(sqrt(abs(s1(gi)))),'r.','markersize',12);
            %plot(rc1nopad(ginopad),dB(abs(s1nopad(ginopad))),'ro')
            rchosen = plot(rc1(ii)+rf1(ii),dB(sqrt(abs(s1(ii)))),'ro','markersize',10);
            %         for ii = 1:length(s1nopad(ginopad))
            %             bch = plot(repmat(rc1nopad(ginopad(ii)),1,2),[-120 0],'k');
            %         end
            % Set grid to range bins
            %set(gca,'xtick',rc1(gi))
            %halfwaves = R1 + (c.lambdac/2)*[min(round((rc1(gi)-R1)./(c.lambdac/2))):max(round((rc1(gi)-R1)./(c.lambdac/2)))];
            %set(gca,'xtick',halfwaves(1:2:end))
            %         for ii = 1:length(halfwaves)
            %             waveh = plot(repmat(halfwaves(ii),1,2),[-120 0],'col',[0.6 0.6 0.6]);
            %         end
            if strcmp(datamode,'synth')
                rtrueh = plot([R1 R1],[-120 0],'g');
            end
            %legend([rtrueh rch rfh rtoth bch waveh rchosen],{'true range','coarse range (range at bin centres)','fine range','total range estimate','unpadded bin centres','integer half-wavelengths from reflector','chosen range'})
            %legend([rch rfh rtoth rchosen],{'coarse range','fine range','total range estimate','chosen range'})
            set(gcf,'tag','fmcw_range_schematic')
            ylabel('Amplitude (dB)')
            xlabel('range (m)')
            ylim(ylims)
            
            linkaxes(ax,'xy')
            text(x1,y1,splab{jj},'fontsize',fonts)
        end
        %grid on
        %keyboard
        
        %         % Now true range vs estimated
        %         figure
        %         %gi = rc1>=(rc1(ii)-5) & rc1<=(rc1(ii)+5);
        %         %ginopad = rc1nopad>=(rc1nopad(iinopad)-5) & rc1nopad<=(rc1nopad(iinopad)+5);
        %         plot(rc1(gi),rc1(gi),'b.')
        %         hold on
        %         plot(rc1(gi),rc1(gi)+rf1(gi),'r.')
        %         quiver(rc1(gi),rc1(gi),zeros(size(rf1(gi))),rf1(gi),0)
        %         hold on
        %         plot(rc1(gi),repmat(R1,size(rc1(gi))),'g')
        %         plot(rc1nopad(ginopad),rc1nopad(ginopad)+rf1nopad(ginopad),'ro')
        %         set(gca,'xtick',rc1(gi))
        %         set(gca,'ytick',halfwaves)
        %         grid on
        
        
    case 'close reflectors'
        % Test the effect of two closely spaced reflectors on phase
        
        % Test generate some profiles
        r1 = 100;
        r2 = 100.05;
        R1 = [10 r1 300 1000];
        A1 = [-30 -40 -80 -20];
        vdat = fmcw_synth(R1,A1);
        save synth_chirp_100.mat
        R2 = [10 r2 300 1000];
        A2 = [-30 -40 -80 -20];
        vdat = fmcw_synth(R2,A2);
        save synth_chirp_100.1.mat
        %fmcw_plot('rap','fl',{'synth_chirp_1.mat','synth_chirp_2.mat'})
        
        % Combo
        R3 = [10 r1 r2 300 1000];
        A3 = [-30 -40 -40 -80 -20];
        vdat = fmcw_synth(R3,A3);
        save synth_chirp_combo vdat
        fmcw_plot('ap','fl',{'synth_chirp_100.mat','synth_chirp_100.1.mat','synth_chirp_combo.mat'})
        %overlay_reflectors(R1,A1,'r.')
        %overlay_reflectors(R2,A2,'g.')
        overlay_reflectors(R3,A3,'b.')
        
        
        % Zoom in
        set(findobj('tag','aax'),'xlim',[95 105])
        
    case 'layering'
        % Very closely spaced identical layers
        
        R = [10:0.1:15, 20:0.2:25 30:0.5:40 50:60];
        A = -80*ones(size(R));
        vdat = fmcw_synth(R,A);
        save synth_chirp_close_layers vdat
        fmcw_plot('tap','fl',{'synth_chirp_close_layers.mat'})
        overlay_reflectors(R,A,'r.')
        % Zoom in
        set(findobj('tag','aax'),'xlim',[0 80])
        
    case 'noise'
        % Noise floor
        R = [100 200 350];
        A = [-40 -60 -20];
        vdat = fmcw_synth(R,A,-80);
        save synth_chirp_noise80 vdat
        vdat = fmcw_synth(R,A,-100);
        save synth_chirp_noise100 vdat
        vdat = fmcw_synth(R,A,-120);
        save synth_chirp_noise120 vdat
        vdat = fmcw_synth(R,A);
        save synth_chirp_noiseNone vdat
        fmcw_plot('tap','fl',{'synth_chirp_noise80.mat','synth_chirp_noise100.mat','synth_chirp_noise120.mat','synth_chirp_noiseNone.mat'})
        overlay_reflectors(R,A,'r.')
        
        % Zoom in
        %set(findobj('tag','aax'),'xlim',[196 204])
        
    case 'clipping'
        % Show effects of insufficient attenuation
        
        % Create realistic profile
        bedRange = 300;
        bedamp = -20;
        npoints = 1000;
        [R,A] = makeprof(bedRange,bedamp,npoints);
        N = -100;
        S = 0.8; % fraction of ADC range to use (no clipping)
        vdat = fmcw_synth(R,A,N,S);
        save no_clip vdat
        S = 1.1; % mild clipping
        vdat = fmcw_synth(R,A,N,S);
        save clip_10percent vdat
        S = 2; % severe clipping
        vdat = fmcw_synth(R,A,N,S);
        save clip_100percent vdat
        fmcw_plot('tap','fl',{'clip_100percent.mat','clip_10percent.mat','no_clip.mat'})
        overlay_reflectors(R,A,'r.')
        
    case 'vsr'
        % Test determination of vertical strain using fmcw_vsr
        doMakeProf = 1;
        if doMakeProf
            % Create realistic profile 1
            bedRange = 300;
            bedamp = -20;
            npoints = 3*bedRange;
            %npoints = 5*bedRange; % assuming annual layering with 20cm water/year
            %npoints = 20;
            [R1,A1] = makeprof(bedRange,bedamp,npoints);
            N = -110;
            %N = nan;
            S = 0; % fraction of ADC range to use (0 = do not scale)
            vdat = fmcw_synth(R1,A1,N,S);
            t1 = vdat.TimeStamp;
            vdat.synth.A = A1;
            vdat.synth.R = R1;
            save synth_chirp_vsr1 vdat
            
            % Make offset and strained version of shot 1
            strainRate = -0.001678; % per year
            dt = 365;
            offset = pi; % accumulation or cables
            disp(['Creating synthetic profile with strain = ' num2str((1+strainRate)*(dt/365),'%6.4f')])
            R2 = R1*(1+strainRate) + offset;
            A2 = A1;
            vdat = fmcw_synth(R2,A2,N,S); % same reflectors in a different location
            vdat.TimeStamp = t1+dt;
            save synth_chirp_vsr2 vdat
            save stuff R1 A1 R2 A2 strainRate offset N dt
            
        end
        
        % Load files
        vdat1 = fmcw_load('synth_chirp_vsr1.mat');
        vdat2 = fmcw_load('synth_chirp_vsr2.mat');
        load stuff
        disp(['Using synthetic profile with strain = ' num2str((1+strainRate)*(dt/365),'%6.4f')])
        
        % Plot
        %fmcw_plot('tap','fl',{'synth_chirp_vsr1.mat','synth_chirp_vsr2.mat'})
        %overlay_reflectors(R1,A1,'b.')
        %overlay_reflectors(R2,A2,'r.')
        
        % Plot rescaled version for viewing differences caused by strain
        % Phase process
        p = 1;
        maxrange = 400;
        [rc1,rf1,sr1,s1] = fmcw_range(vdat1,p,maxrange);
        r1 = rc1+rf1; % total range
        [rc2,rf2,sr2,s2] = fmcw_range(vdat2,p,maxrange);
        r2 = rc2+rf2; % total range
        
        figure
        % Amplitude
        ax(1) = subplot(2,1,1);
        plot(R1,A1,'b.')
        hold on
        plot(rc1,dB(abs(s1)),'b')
        % now overlay shot - shifted so that we can compare visually with original
        %R2 = R1*(1+strainRate)*(dt/365) + offset;
        rc2c = (rc2 - offset)/(1+strainRate) ; % rescaled so we can visually compare
        %plot(R2,A2,'r.')
        plot(rc2,dB(abs(s2)),'r')
        plot(rc2c,dB(abs(s2)),'g')
        %legend('reflectors 1','profile 1','rescaled profile 2')
        %legend('reflectors 1','profile 1','reflectors 2','profile 2','rescaled profile 2')
        shift = rc2*strainRate + offset;
        s2s = fmcw_shift(rc2,rf2,s2,-shift,vdat2,p);
        plot(rc2,dB(abs(s2s)),'k')
        legend('reflectors','profile','strained','rescaled (no phase shift)','unshifted')
        title(['strain = ' num2str((1+strainRate)*(dt/365),'%6.4f') ' + offset ' num2str(offset) 'm'])
        
        % Phase
        ax(2) = subplot(2,1,2);
        plot(rc1,angle(s1),'b')
        hold on
        plot(rc2,angle(s2),'r')
        plot(rc2c,angle(s2),'g')
        plot(rc2,angle(s2s),'k')
        legend('profile','strained','rescaled (no phase shift)','unshifted')
        %legend('profile 1','profile 2','rescaled profile 2')
        linkaxes(ax,'x')
        
        % Calculate strain rate
        %cfg = fmcw_process_config_default; % load default settings
        cfg.doClean = 0; % Turn off data cleaning
        cfg.maxStrain = abs(strainRate)*1.1;
        cfg.coarseChunkWidth = 15;
        cfg.minCohereCoarse = 0.75; % minimum correlation on lag offset (i.e. confidence we have the right lag)
        cfg.minCohereFine = 0.75; % minimum correlation to use in strain estimate
        cfg.minCoherePhase = 0.3; % minimum correlation to use in strain estimate
        vsr = fmcw_vsr('synth_chirp_vsr1.mat','synth_chirp_vsr2.mat',cfg);
        disp(' ')
        disp(['Actual strain rate   : ' num2str(1000*strainRate,'%6.4f') ' milli strain'])
        disp(['Estimated strain rate: ' num2str(1000*vsr.vsr,'%6.4f') ' milli strain'])
        
        
        
        %keyboard
        
    case 'depth shift'
        % investigate the effect of a half range bin shift on reflectors
        %
        % result: no effect
        
        doMakeProf = 0;
        if doMakeProf
            
            % Create realistic profile 1
            bedRange = 300;
            bedamp = -20;
            npoints = 300;
            [R1,A] = makeprof(bedRange,bedamp,npoints);
            %N = -110;
            N = nan;
            S = 0; % fraction of ADC range to use
            vdat = fmcw_synth(R1,A,N,S);
            save synth_chirp_raw vdat
            
            % Create depth shifted version of above
            % make depths shift = half of the raw (unpadded) range resolution
            % Radar parameters to define coarse range resolution
            er = 3.1;
            c = 3e8; % velocity
            B = 2e8; % bandwidth
            dr = c/(2*B*sqrt(er)); % Brennan eq 14 rearranged
            %depthShift = dr/2;
            depthShift = pi;
            
            R2 = R1 + depthShift;
            vdat = fmcw_synth(R2,A,N,S);
            vdat.synth.depthShift = depthShift;
            save synth_chirp_offset vdat
        end
        
        % Load data
        vdat1 = fmcw_load('synth_chirp_raw.mat');
        vdat2 = fmcw_load('synth_chirp_offset.mat');
        depthShift = depthShift;
        % Phase process
        p = 10;
        maxrange = 400;
        [rc1,~,sr1,s1] = fmcw_range(vdat1,p,maxrange);
        [rc2,~,sr2,s2] = fmcw_range(vdat2,p,maxrange);
        
        % Plot
        figure
        % Amplitude
        ax(1) = subplot(2,1,1);
        plot(rc1,dB(abs(s1)),'b')
        hold on
        % now overlay shot - shifted so that we can compare visually with original
        plot(rc2-depthShift,dB(abs(s2)),'r')
        s2s = fmcw_shift(rc2,s2,-depthShift,vdat1.fc,vdat1.K,vdat1.ci);
        plot(rc2,dB(abs(s2s)),'g')
        legend('raw',['resampled with ' num2str(depthShift) 'm shift'])
        title(['Effect of depthshift'])
        
        % Phase
        ax(2) = subplot(2,1,2);
        plot(rc1,angle(s1),'b')
        hold on
        plot(rc2-depthShift,angle(s2),'r')
        plot(rc2,angle(s2s),'g')
        legend('raw',['resampled with ' num2str(depthShift) 'm shift'])
        
        linkaxes(ax,'x')
        
    case 'deconv'
        % Investigate deconvoluting the range profile by the blackman
        % window to identify individual reflectors
        
        % can't make this work!
        
        % Create realistic profile 1
        bedRange = 300;
        bedamp = -20;
        npoints = 10;
        [R1,A] = makeprof(bedRange,bedamp,npoints);
        %N = -110;
        N = nan;
        S = 0; % fraction of ADC range to use
        vdat = fmcw_synth(R1,A,N,S);
        save synth_chirp_raw vdat
        
        vdat = fmcw_load('synth_chirp_raw.mat');
        % Phase process
        p = 2;
        maxrange = 400;
        winfunh = @blackman;
        [rc1,~,sr,s] = fmcw_range(vdat,p,maxrange,winfunh);
        
        % Deconvolute
        % Make window (odd length)
        n = 5;
        w  = window(@blackman,n);
        wpad = [w; zeros(n*(p-1),1)];
        xn = 0.5*(n-1);
        wpads = circshift(wpad,[-xn 0]); % phase centre
        a = fft(wpads);
        a = fftshift(a);
        win = a(2:end); % remove mean
        % now pad to length of signal
        %winpad = zeros(size(rc1));
        %winpad(1:length(win)) = win;
        %win = winpad;
        
        %figure
        %plot(win)
        
        %[da,dr] = deconv(s.*transpose(blackman(length(s))),win);
        [da,dr] = deconv(s,win);
        nanpad = nan*ones(1,(numel(win)-1)/2);
        dap = [nanpad da nanpad];
        drp = [nanpad dr nanpad];
        
        % Plot
        figure
        % Amplitude
        ax(1) = subplot(2,1,1);
        plot(R1,A,'r.')
        hold on
        plot(rc1,dB(abs(s)),'b')
        plot(rc1,dB(abs(dap)),'g.')
        
        legend('raw','processed','deconvoluted')
        title(['Effect of depthshift'])
        
        % Phase
        ax(2) = subplot(2,1,2);
        plot(rc1,angle(s),'b')
        hold on
        plot(rc1,angle(dap),'g.')
        legend('processed','deconvoluted')
        
        linkaxes(ax,'x')
        
        figure
        subplot(2,1,1)
        plot(abs(win))
        subplot(2,1,2)
        plot(angle(win))
        
        figure
        %plot(dB(win))
        
        keyboard
        
    case 'bed shift'
        % Test fmcw_melt without strain
        doMakeProf = 1;
        if doMakeProf
            % Create realistic profile 1
            bedRange = 300;
            bedamp = -20;
            npoints = 3*bedRange;
            %npoints = 5*bedRange; % assuming annual layering with 20cm water/year
            %npoints = 20;
            [R1,A1] = makeprof(bedRange,bedamp,npoints);
            N = -110;
            %N = nan;
            S = 0; % fraction of ADC range to use (0 = do not scale)
            vdat = fmcw_synth(R1,A1,N,S);
            t1 = vdat.TimeStamp;
            vdat.synth.A = A1;
            vdat.synth.R = R1;
            save synth_chirp_vsr1 vdat
            
            % Make offset and strained version of shot 1
            strainRate = -0.01678; % per year
            dt = 365;
            offset = pi; % accumulation or cables
            disp(['Creating synthetic profile with strain = ' num2str((1+strainRate)*(dt/365),'%6.4f')])
            R2 = R1*(1+strainRate) + offset;
            A2 = A1;
            vdat = fmcw_synth(R2,A2,N,S); % same reflectors in a different location
            vdat.TimeStamp = t1+dt;
            save synth_chirp_vsr2 vdat
            save stuff R1 A1 R2 A2 strainRate offset N dt
            
        end
        
        % Load files
        vdat1 = fmcw_load('synth_chirp_vsr1.mat');
        vdat2 = fmcw_load('synth_chirp_vsr2.mat');
        load stuff
        disp(['Using synthetic profile with strain = ' num2str((1+strainRate)*(dt/365),'%6.4f')])
        
        
        
    case 'melt'
        % Test fmcw_melt without strain
        % Test determination of vertical strain using fmcw_vsr
        doMakeProf = 1;
        if doMakeProf
            % Create realistic profile 1
            bedRange = 260;
            bedamp = -20;
            npoints = 3*bedRange;
            %npoints = 5*bedRange; % assuming annual layering with 20cm water/year
            %npoints = 20;
            [R1,A1] = makeprof(bedRange,bedamp,npoints);
            %N = -90;
            N = nan;
            S = 0; % fraction of ADC range to use (0 = do not scale)
            vdat = fmcw_synth(R1,A1,N,S);
            t1 = vdat.TimeStamp;
            vdat.synth.A = A1;
            vdat.synth.R = R1;
            save synth_chirp_melt1 vdat
            
            % Make offset and bedshifted version of shot 1
            strainRate = 0.001; % 0.01234567; % per year
            daysPerYear = 365.25;
            dt = 365.25;
            strain = strainRate*dt/daysPerYear;
            offset = 0; % accumulation or cables
            bedinds = find(R1>=bedRange);
            meltRate = 1;
            melt = meltRate*dt/daysPerYear;
            % whole profile strain then shift
            R2 = R1*(1+strain) + offset; % depth shift the whole thing (accumulation)
            % bed - melt
            R2(bedinds) = R2(bedinds) - melt; % shift the bed and everything below
            %R2(bedinds(2:end)) = R2(bedinds(2:end)) + 0.3*randn(size(R2(bedinds(2:end)))); % add some noise
            A2 = A1;
            %A2(bedinds) = A2(bedinds) + randn(size(A2(bedinds))); % a little noise on bed amp
            vdat = fmcw_synth(R2,A2,N,S); % same reflectors in a different location
            vdat.TimeStamp = t1+dt;
            save synth_chirp_melt2 vdat
            save meltstuff R1 A1 R2 A2 strainRate meltRate offset melt N dt
            
        end
        
        % Load files
        vdat1 = fmcw_load('synth_chirp_melt1.mat');
        vdat2 = fmcw_load('synth_chirp_melt2.mat');
        load meltstuff
        disp(['Using synthetic profile with strain = ' num2str(strain,'%6.4f')])
        disp(['Using synthetic profile with melt = ' num2str(melt,'%6.4f')])
        
        % Process for melt
        cfg.doClean = 0;
        cfg.errorMethod = 'assumedNoiseFloor'; % 'emperical' 'assumedNoiseFloor'
        cfg.chunkWidth = 4;
        cfg.bedSearchRange = [100 inf];
        cfg.ampThreshdB = -60;
        cfg.doBulkAllignment = 1; %
        cfg.minDepth = 10;
        cfg.firnDepth = 60;
        
        cfg.bedShiftMethod = 'rangeDiff'; % 'xcorr' 'rangeDiff' (xcorr more reliable...)
        ard = fmcw_melt('synth_chirp_melt1.mat','synth_chirp_melt2.mat',cfg)
        
        %cfg.bedShiftMethod = 'xcorr'; % 'xcorr' 'rangeDiff' (xcorr more reliable...)
        %axc = fmcw_melt('synth_chirp_melt1.mat','synth_chirp_melt2.mat',cfg)
        
        keyboard
        
    case 'bed roughness'
        % Test the effect of multiple closely spaced reflectors on phase
        % (i.e. the effect of a rough bed)
        
        % Test generate some profiles
        Anom = -80;
        astd = 5;
        n = 100;
        Rnom = 300; % nominal depth
        Rstep = 5; % separation between peaks so we can see them
        rstd = [0.01 0.02 0.04 0.06];
        
        for ii = 1:length(rstd)
            % Roughness 0
            R = Rnom + (ii-1)*Rstep + rstd(ii)*randn(1,n); % reflector ranges
            A = Anom + astd*randn(1,n); % reflector amplitudes
            vdat = fmcw_synth(R,A); % make synthetic return
            save synth_chirp_roughbed.mat
            fmcw_plot('ap','fl','synth_chirp_roughbed.mat','nl')
            h = overlay_reflectors(R,A,'k.');
            set(h,'HandleVisibility','off')
            lt(ii) = {[num2str(rstd(ii)) ' m']};
        end
                
        % Zoom in
        set(findobj('tag','aax'),'xlim',[Rnom-Rstep Rnom+ii*Rstep])
        axes(findobj('tag','aax'))
        legend(lt)
        
    case 'system gain'
        % Test the effect of frequency dependent system gain
        % e.g. high frequency attenuation often seen in field data
        
        % Test generate a profile
        A = -20;
        Rnom = 300; % nominal depth
        Rstep = 5; % separation between peaks so we can see them
        sg400 = [1 0.8 0.6 0.4 0.2]; % system gain at 400MHz
        
        for ii = 1:length(sg400)
            % Roughness 0
            R = Rnom + (ii-1)*Rstep; % reflector ranges
            vdat = fmcw_synth(R,A); % make synthetic return
            save synth_chirp_sys_at.mat
            % Reload so that we turn vdat.v - into vdat.vif
            vdat = fmcw_load('synth_chirp_sys_at.mat',1);
            calArray = [1 sg400(ii)]; % system gain
            vdat = fmcw_cal(vdat,calArray);
            vdat.v = vdat.vif; % have to overwrite v with modified signal as this is loaded again in fmcw_plot
            save synth_chirp_sys_at.mat
            
            fmcw_plot('tap','fl','synth_chirp_sys_at.mat','nl')
            lt(ii) = {['sg@400MHz: ' num2str(sg400(ii)) ]};
        end
                
        % Zoom in
        set(findobj('tag','aax'),'xlim',[Rnom-Rstep Rnom+ii*Rstep])
        axes(findobj('tag','aax'))
        legend(lt)
end

function h = overlay_reflectors(R,A,marker)
% Plot specified reflectors on aplitude profile
aax = findobj('tag','aax');
axes(aax);
hold on
h = plot(R,A,marker)

function [R,A] = makeprof(bedRange,bedamp,n)
% Make "realistic" looking profile
%
% n = number of reflectors

R = (bedRange+100) * rand(1,n); % more layers than we can resolve
% Kill those close to bed
R(abs(R-bedRange)<1) = [];
R = [R bedRange]; % include specified bed range
range = [0    20  65  90 bedRange-0.1 bedRange bedRange+20 bedRange+100];
amp =   [-50 -30 -40  -70   -120        bedamp  bedamp-70 bedamp-100] -dB(n.^(1/4)); % trying to scale for num samples...
A = interp1(range,amp,R);
ampstd = 5;
A = A + ampstd*randn(size(A)); % random variation in reflector amp

