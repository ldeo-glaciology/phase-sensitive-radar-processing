function PD = fmcw_kwnprofile(varargin)
%{
 ProfileDescriptor structure
 PD         Time-ordered structure array with fields:
    Fnam    Filename
    Vt      VM2 Time of first burst
    Burst   Burst number within file
    T1      Temperature 1
    T2      Temperature 2
    Batt    Battery voltage
    Lat     Latitude
    Lon     Longitude
    Gt      GPS time
%}

Ncartoon = 1000; % Number of burst averages to plot in the cartoon
maxRange = 1000;
win = 0;

if nargin == 1
    if isempty(varargin{1})
        [PDfile, PDpath] = uigetfile('*.mat','Select profile descriptor file');
        load([PDpath, PDfile]);
    else
        load(varargin{1});
    end
    
    % Check for valid descriptor file
    if ~isfield(PD,'Vt')
        
        % If not valid, message user and exit
        fprintf('Input file not a valid profile descriptor - exiting\n');
        return
    end
else
    
    % Obtain set of filenames
    [filenames, path] = uigetfile('*.dat','Enter data filenames','Multiselect','on');
    if ~iscell(filenames)
        tmp = filenames;
        filenames = cell(1);
        filenames{1} = tmp;
    end
    
    PD = struct();
    % Analyse files to create descriptor structure
    Nfiles = length(filenames);
    PD.Ntotal = 0;  % Running total of Bursts
    for File = 1:Nfiles
        for Burst = 1:1e6
            vdat = fmcw_load([path, filenames{File}],Burst);
            if vdat.Code == -4
                break
            elseif vdat.Code == 0
                PD.Ntotal = PD.Ntotal + 1;
                % Populate PD
                PD.Fnam{PD.Ntotal} = vdat.filename;
                PD.Vt(PD.Ntotal) = vdat.TimeStamp;
                PD.Burst(PD.Ntotal) = vdat.Burst;
                PD.T1(PD.Ntotal) = vdat.Temperature_1;
                PD.T2(PD.Ntotal) = vdat.Temperature_2;
                PD.Batt(PD.Ntotal) = vdat.BatteryVoltage;
                PD.NChirps(PD.Ntotal) = vdat.ChirpsInBurst;
                
                %                 PD(PD.Ntotal).Lat(PD.Ntotal) = vdat.latitude;
                %                 PD(PD.Ntotal).Lon(PD.Ntotal) = vdat.longitude;
                %                 PD(PD.Ntotal).Gt(PD.Ntotal) = vdat.GPStime
            else
                fprintf('Error reading Burst %d from file %s\n',Burst,filenames{File});
                break
            end
        end
    end
    
    % Sort on time
    [PD.Vt, I] = sort(PD.Vt);
    PD.Fnam = PD.Fnam(I);
    PD.Burst = PD.Burst(I);
    PD.T1 = PD.T1(I);
    PD.T2 = PD.T2(I);
    PD.Batt = PD.Batt(I);
    PD.NChirps = PD.NChirps(I);
    %     PD.Lat = PD.Lat(I);
    %     PD.Lon = PD.Lon(I);
    %     PD.Gt = PD.Gt(I);
    clear I;
    
    [p,f,e] = fileparts(PD.Fnam{1});
    eval(['save ' [p, '\PD_',f, '.mat'], ' PD '])
    
end


%% Plot cartoon

% Obtain chirp characteristics
vdat = fmcw_load(PD.Fnam{1},1);
bin2range = vdat.ci/2/vdat.B; range2bin = 1/bin2range;
range = [1:size(vdat.vif,2)] * bin2range;
if PD.Ntotal > 1
    PD.Dt = mean(diff(PD.Vt)./PD.NChirps(1:end-1));
else
    PD.Dt = vdat.T/24/3600;
end

maxBin = round(maxRange*range2bin);
Spc = zeros(maxBin,min(PD.Ntotal,Ncartoon)); Ang = Spc;
tim = zeros(1,min(PD.Ntotal,Ncartoon));

selBurst = 0;
dBurst = max(floor(PD.Ntotal/Ncartoon),1);
for Burst = 1:dBurst:PD.Ntotal
    selBurst = selBurst + 1;
    vdat = fmcw_load(PD.Fnam{Burst},PD.Burst(Burst));
    tmp = fft(mean(vdat.vif,1).*blackman(vdat.SamplesPerChirp)');
    Spc(:,selBurst) = 20*log10(abs(tmp(1:maxBin)));
    Ang(:,selBurst) = angle(tmp(1:maxBin));
    tim(selBurst) = PD.Vt(Burst);
end
CartoonHan = figure;
imagesc(1:selBurst,range(1:maxBin)',Spc);
colormap jet
uAng = detrend(unwrap(Ang,[],2),'constant');
PhaseHan = figure;
imagesc(1:selBurst,range(1:maxBin)',uAng);
colormap jet

%datetick('x','HH:MM:SS','keeplimits')

%%
while 1
    if exist('linehan1'), delete(linehan1); end
    if exist('linehan2'), delete(linehan2); end
    % Request bracket
    figure(CartoonHan);
    fprintf('Click limits for detailed plot, right click to exit\n')
    [B(1), tmp, button] = ginput(1);
    if button ~= 1, break; end
    yl = ylim; xl = xlim;
    linehan1 = line([B(1),B(1)],yl,'color',[1,0,0]);
    [B(2), tmp, button] = ginput(1);
    if button ~= 1, break; end
    linehan2 = line([B(2),B(2)],yl,'color',[1,0,0]);
    B = sort(B);
    B = B * dBurst;
    
    % Plot full profile figure over requested bracket
    %B(1) = find(PD.Vt >= B(1),1);
    %B(2) = find(PD.Vt >= B(2),1);
    B = B - 0.5;
    burst1 = ceil(B(1));
    burst2 = ceil(B(2));
    
    frac1 = B(1) - floor(B(1));
    frac2 = B(2) - floor(B(2));
    
    fc = round(frac1 * PD.NChirps(burst1));
    FirstChirp = fc;
    if burst1 > 1
        FirstChirp = sum(PD.NChirps(1:burst1-1)) + fc;
    end
    lc = round(frac2 * PD.NChirps(burst2));
    LastChirp = lc;
    if burst2 > 1
        LastChirp = sum(PD.NChirps(1:burst2-1)) + lc;
    end
    
    Spc = zeros(maxBin, LastChirp - FirstChirp+1);
    a = zeros(maxBin, LastChirp - FirstChirp+1);
%    vif = zeros(LastChirp - FirstChirp+1,vdat.SamplesPerChirp);
    t = zeros(1, LastChirp - FirstChirp+1);
    
    chirpCount = 1;
    
    vdat = fmcw_load(PD.Fnam{burst1},PD.Burst(burst1));
    if burst2 == burst1
%        vif(chirpCount:(LastChirp-FirstChirp+1),:) = vdat.vif(FirstChirp:LastChirp,:);
        for chirp = fc:lc
            t(chirpCount) = PD.Vt(burst1) + PD.Dt * (chirp-1);
            spectr = fft(vdat.vif(chirp,:).*blackman(vdat.SamplesPerChirp)');
            tmp = 20*log10(abs(spectr));
            Spc(:,chirpCount) = tmp(1:maxBin);
            a(:,chirpCount) = angle(spectr(1:maxBin));
            chirpCount = chirpCount + 1;
        end
    else
%        vif(FirstChirp:PD.NChirps(burst1),:) = vdat.vif(FirstChirp:PD.NChirps(burst1),:);
        for chirp = fc:PD.NChirps(burst1)
            t(chirpCount) = PD.Vt(burst1) + PD.Dt * (chirp-1);
            spectr = fft(vdat.vif(chirp,:).*blackman(vdat.SamplesPerChirp)');
            tmp = 20*log10(abs(spectr));
            Spc(:,chirpCount) = tmp(1:maxBin);
            a(:,chirpCount) = angle(spectr(1:maxBin));
            chirpCount = chirpCount + 1;
        end
        for burst = burst1+1:burst2-1
            
            vdat = fmcw_load(PD.Fnam{burst},PD.Burst(burst));
 %           vif(chirpCount:(chirpCount + PD.NChirps(burst)-1),:) = vdat.vif(1:PD.NChir ps(burst),:);
            for chirp = 1:PD.NChirps(burst)
                t(chirpCount) = PD.Vt(burst) + PD.Dt * (chirp-1);
                spectr = fft(vdat.vif(chirp,:).*blackman(vdat.SamplesPerChirp)');
                tmp = 20*log10(abs(spectr));
                Spc(:,chirpCount) = tmp(1:maxBin);
                a(:,chirpCount) = angle(spectr(1:maxBin));
                chirpCount = chirpCount + 1;
            end
        end
%        vif(chirpCount:(chirpCount + LastChirp-1),:) = vdat.vif(1:LastChirp,:);
        vdat = fmcw_load(PD.Fnam{burst2},PD.Burst(burst2));
        for chirp = 1:lc
            t(chirpCount) = PD.Vt(burst2) + PD.Dt * (chirp-1);
            spectr = fft(vdat.vif(chirp,:).*blackman(vdat.SamplesPerChirp)');
            tmp = 20*log10(abs(spectr));
            Spc(:,chirpCount) = tmp(1:maxBin);
            a(:,chirpCount) = angle(spectr(1:maxBin));
            chirpCount = chirpCount + 1;
        end
    end
    
%     AvgSpc = zeros(maxBin,chirpCount - win);
%     for chirp = 1:chirpCount-win-1
%         tmp = 20*log10(abs(fft(mean(vif(chirp:chirp+win,:),1).*blackman(vdat.SamplesPerChirp)')));
%         AvgSpc(:,chirp) = tmp(1:maxBin);
%     end
    
%%
%     figure; imagesc(t(1:chirpCount-win-1),range(1:maxBin)',AvgSpc)
%     datetick('x','HH:MM:SS','keeplimits')
    
    figure('units','normalized','position',[0.04,0.1,0.92,0.6]); imagesc(t,range(1:maxBin)',Spc)
    colormap jet
    if t(end)-t(1) < 0.5
        datetick('x','HH:MM','keeplimits')
    elseif t(end)-t(1) < 2
        datetick('x','HH','keeplimits')
    elseif t(end)-t(1) < 10
        datetick('x','dd/mm','keeplimits')
    else
        datetick('x','dd/mm','keeplimits')
    end
    
%%

%     sz = size(a);
% 
%     % Unwrap along rows (time) dimension
%     una = unwrap(a,[],2);
% 
%     %remove mean
%     muna = una - repmat(mean(una,2),1,sz(2));
% 
%     % Interpolate to 5-minute intervals
%     %t = Info.Time(1):1/12/24:Info.Time(end);
%     %imuna = interp1(Info.Time,muna',t);
%     %ib = interp1(Info.Time,b',t);
% 
%     % Assuming really noisy unwrapping means low signal strength, use variance
%     % with a threshold to select for decent records
%     tol = 0.02;
%     vmuna = zeros(size(muna,1),1);
%     for i = 1:size(muna,1);
%         %vimuna = var(detrend(imuna','linear'),0,2);
%         vmuna(i) = var(cheb(muna(i,:),0.1,0.99));
%     end
%     ind = find(vmuna > tol);
% 
%     % Replace noisy phase records with nan
%     seluna = muna;
%     clear muna
%     seluna(ind,:) = nan;
%     %selib = ib;
%     %selib(:,ind) = nan;
% 
%     % get depth estimate to go with it (assuming 3.18, 200e6, 300e6 and pad=2)
%     Z = (0:sz(2)-1) * 0.2106;
% 
%     figure;imagesc(seluna); colormap jet
%     caxis([-1.8,0.5]);
%%
    fprintf('Any key to return to cartoon\n')
    pause
end

