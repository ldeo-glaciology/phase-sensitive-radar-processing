function fmcw_profile(filelist)

% fmcw_profile(filelist)
%
% Plots fmcw data in profile mode
% Craig Stewart
% 2013-1-7 (Scott Base)

% Settings
padfactor = 2;
maxrange = 2000;
nthBurst = 1;
nthChirp = 10;
samplingList = {'allChirps','nthChirp','meanBurst'};
sampling = samplingList{1};
doPlotTemp = 1;

if nargin == 0
    [filelist, pathname] = uigetfile(['*.dat;*.DAT;*.000'],'Choose radar files to plot','multiselect','on');
    if isa(filelist,'double') % no files chosen
        return
    end
    if ischar(filelist) % cover the case of a single file selected (where uigetfile returns a char
        filelist = {filelist};
    end
    for ii = 1:length(filelist) % Convert all filenames to path/file for consistency
        filelist(ii) = {[pathname filelist{ii}]};
    end
end
if ischar(filelist) % cover the case of a single file selected (where uigetfile returns a char
    filelist = {filelist};
end

% 
% % Calculate how many profiles we;re keeping
% for FileNo = 1:length(filelist)
%     filename = filelist{FileNo};
%     nBursts(FileNo) = fmcw_nbursts(filename);
%     numBurstsUsed(FileNo) = length(1:nthBurst:nBursts(FileNo));
% end
% N = sum(numBurstsUsed);


A = [];
BurstTime = [];
ChirpTime = [];
Tem1 = [];
Tem2 = [];
for FileNo = 1:length(filelist)
    filename = filelist{FileNo};
    nBursts = fmcw_nbursts(filename);
    for burstNo = 1:nthBurst:nBursts
        vdat = fmcw_load(filename,burstNo); % load data
        [~,name,~] = fileparts(filename);
        disp(['Loaded file ' name ', burst ' int2str(burstNo)]);
       
        % Subsample
        switch sampling
            case 'nthChirp'
                chirpList = [1:nthChirp:vdat.ChirpsInBurst];
                vdat = fmcw_burst_subset(vdat,chirpList);
            case 'meanBurst'
                vdat = fmcw_burst_mean(vdat);
        end

        % phase process data
        [rc,rf,sr,su] = fmcw_range(vdat,padfactor,maxrange); %ELIZ: this function only seems to have three outputs, but requesting 4

        a = 20*log10(abs(sr));
        
        % Concatenate with other bursts
        ChirpTime = [ChirpTime; vdat.chirpTime];
        BurstTime = [BurstTime; vdat.TimeStamp];
        A = [A a'];
        Tem1 = [Tem1 vdat.Temperature_1];
        Tem2 = [Tem2 vdat.Temperature_2];
    end
end

% Check that there is more than one chirp
if size(A,2)==1
    disp('Only one chirp on profile: plotting cancelled')
    return
end

% Plot profile
%chirpnum = 1:size(A,2);
figure
if doPlotTemp
    subplot(2,1,1)
end
%sh = surf((ChirpTime-ChirpTime(1))*24*60,rc,A,'edgecolor','none');

sh = surf([1:size(A,2)],rc,A,'edgecolor','none');
hold on
view(0,90)
colorbar('East')
set(gca,'ydir','rev')
axis tight
%xlabel('ChirpTime from start (minutes)')
xlabel('chirp')
ylabel('range (m)')
%save profile ChirpTime rc A
[PATHSTR,NAME,EXT] = fileparts(filelist{1});
if length(filelist)==1
    title(NAME,'interpreter','none')
end
set(gcf,'tag',[NAME '_ZScope'])

if doPlotTemp
% Plot instrument temperature

%figure
subplot(2,1,2)
plot(BurstTime,Tem1,'r',BurstTime,Tem2,'b')
ylabel('Temperature')
legend('Temp1','Temp2')
axis tight
xtimelab
end