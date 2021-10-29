function fmcw_plot_chirp_power(file,burst,noisePowerLimit)

% fmcw_plot_chirp_power(file,burst,noisePowerLimit)
% 
% Plot chirp power throughout a burst, identifying outliers

% Craig Stewart
% 2014/5/30

%% Select and load data
if nargin == 0
    [file,~] = uigetfile({'*.dat;*.DAT;*.000;*.mat','Radar files: .dat, .DAT, .000 and .mat'},'Choose radar files to plot power','multiselect','on');
    if isa(file,'double') % no files chosen
        return
    end
end
if isa(file,'char')
    file = {file};
end

if nargin <2
    disp('Showing burst 1 only')
    burst = 1; %you can plot more than 1 if you want
end
if nargin < 3
    noisePowerLimit = 0.001; % fraction of median noise power
end


for ii = 1:length(file)
    for bi = 1:length(burst)
        vdat = fmcw_load(file{ii},burst(bi));
        v = vdat.vif;
        mv = repmat(mean(v),size(v,1),1); % mean of each sample (across all chirps)
        ap = rms(v-mean(v(:)),2); % absolute rms voltage
        np = rms(v-mv,2); % anomaly from burst mean rms voltage
        cn = [1:size(v,1)]'; % chirps number
        
        % Highlight noisey chirps
        % Olde method
        %     noisePowerLimit = 0.01;
        %     noisey = p > min(p)*(1+noisePowerLimit);
        %     plot(cn(noisey),p(noisey),'r.')
        
        % New method using robustfit
        %[b,stats] = robustfit([cn cn.^2],p); % robust quadratic fit to power curve
        %bpred = b(1) + b(2).*cn + b(3).*cn.^2; % prediction from fit
        
        %noisey = abs(stats.resid) > median(p)*noisePowerLimit; % noisey from residual to fit
        
        % New method using robustfit
        fitType = 'exp'; %'lin','quad','cub','exp';
        switch fitType
            case 'quad'
                [b,stats] = robustfit([cn cn.^2],np); % robust quadratic fit to power curve
                bpred = b(1) + b(2).*cn + b(3).*cn.^2; % prediction from fit
                resid = stats.resid;
            case 'exp'
                ft=fittype('exp1');
                f = fit(cn,np,ft);
                bpred = f.a*exp(f.b*cn); % predicted values
                resid = np - bpred;
        end
        
        noisey = abs(resid) > median(np)*noisePowerLimit; % noisey from residual to fit
        
        % Plot noise power per chirp
        figure
        ax(1) = subplot(3,1,1);
        plot(cn,ap,'r') % absolute power
        xlabel('chirp')
        ylabel('chirp power (rms volts)')
        
        ax(2) = subplot(3,1,2);
        hold on
        plot(cn,np,'b') % noise power
        plot(cn,bpred,'g') % fit
        plot(cn(noisey),np(noisey),'go') % noisey
        xlabel('chirp')
        ylabel('noise power (rms volts)')
        title([file{ii} ': burst: ' int2str(burst(bi)) ' chirp power'],'interpreter','none')
        
        frq = 1:1:length(mv(1,:));
        for i=1:length(mv(1,:))
        figure(3)
        plot(frq(i),mv(i,1),'o') % absolute power
        hold on;
        xlabel('chirp')
        ylabel('mean power (rms volts)')
        end
%         
%         figure
%         [n,x] = hist(resid);
%         hold on
%         %bar(x,n,'FaceColor','none')
%         plot(x,n,'.')
%         xlabel(['noise power residual from ' fitType ' fit'])
%         ylabel('num occurance')
    end
end
