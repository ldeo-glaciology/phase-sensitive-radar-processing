function fmcw_meanchirp(Nchirps,filelist)

% fmcw_meanchirp_oneatt(Nchirps,filelist)
% 
% Calculate mean chirp from chirps spread over multiple bursts and files
% saves mean chirp in new .mat file.
%
% Uses chirp stacking (to burst level) during read to conserve memory
%
% args in:
% Nchirps: number of chirps to average
% filelist: list of files to use chirps from
%
% note: if filelist contains fewer chirps than requested this will use all
% the chirps available.
%
% Craig Stewart
% 2013-10-02
% 2014/12/3: added optional cleaning to remove noisey shots from each burst

if nargin == 0
    %error('Not enough args: must specify number of chirps to average')
    Nchirps = input('Enter the number of chirps to average (may be vector): ');
end
if nargin < 2
%    datafid = fopen('fmcw_data_dir.txt','r');
%    if datafid == -1
        initpath = pwd;
%    else
%        initpath = fgets(datafid); fclose(datafid);
%    end
    %[filelist, pathname] = uigetfile([initpath,'\*.dat;*.DAT;*.000;*.mat'],'Choose radar file to use chirps from','multiselect','on');
    [filelist, pathname] = uigetfile('*.dat;*.DAT;*.000;*.mat','Choose radar file to use chirps from','multiselect','on');
%    datafid = fopen('fmcw_data_dir.txt','w'); fprintf(datafid,'%s',pathname); fclose(datafid);
end
if isa(filelist,'char')
    filelist = {filelist};
end

SamplesPerChirp = 40000;

% Get filename root
% outputfilenameroot = input('Enter name of average set (default=''test''): ','s');
% if isempty(outputfilenameroot)
%     outputfilenameroot = 'test';
% end
outputfilenameroot = 'meanchirp'; % just use standard name for output files

% Loop through number to average
for nci = 1:length(Nchirps)
    disp(' ')
    disp(['Averaging chirps for N = ' int2str(Nchirps(nci))])
    clear file burstno v n
    
    % Loop through files
    nchirpstoget = Nchirps(nci);
    ii = 0; % master burst number
    for fn = 1:length(filelist)
        Filename = [pathname,filelist{fn}];
        fmt = fmcw_file_format(Filename);
        Burst = 0;
        
        vdatt.Code = 0;
        while nchirpstoget > 0 && vdatt.Code == 0 % i.e. still need chirps and last read was good
            Burst = Burst + 1;
            % Loop through bursts
            doClean = false; % slower but weeds out bad chirps on the fly
            if doClean
                vdatt = fmcw_load(Filename,Burst);
                if vdatt.Code == 0
                    vdatt = fmcw_keep_clean_chirps(vdatt,round(size(vdatt.vif,1)*0.7),0); % best 70% of chirps
                    if size(vdatt.vif,1)>nchirpstoget
                        vdatt.vif = vdatt.vif(1:nchirpstoget,:);
                    end
                    % Now burst stack for consistency with LongBurst case below
                    vdatt.v = sum(vdatt.vif);
                    vdatt.NChirp = size(vdatt.vif,1);
                end
            else % memory efficient loading of long bursts using burst stacking (cleaning not supported)
                if fmt==3
                    vdatt = LongBurstRMB3(Filename, Burst, SamplesPerChirp, nchirpstoget);
                elseif fmt ==4
                    vdatt = LongBurstRMB4(Filename, Burst, SamplesPerChirp, nchirpstoget);
                elseif fmt ==5
                    vdatt = LongBurstRMB5(Filename, Burst, SamplesPerChirp, nchirpstoget);
                else
                    error('Unsupported file format')
                end
            end
            if vdatt.Code == 0
                vdat = vdatt; % Keep temporary variable if we've had a good read
                
                ii = ii + 1;
                file(ii) = {Filename};
                burstno(ii) = Burst;
                v(ii,:) = vdat.v; % stacked NChirp chirps from this burst
                n(ii) = vdat.NChirp; % number of chirps used
                nchirpstoget = nchirpstoget - vdat.NChirp;
                time(ii) = vdat.TimeStamp;
                disp(['Loaded ' int2str(vdat.NChirp) ' chirps from file ' filelist{fn} ' burst ' int2str(Burst)])
            end
        end
    end
    nchirpsread = Nchirps(nci) - nchirpstoget;
    
    % Now find average of all chirps
    vdat.v = sum(v,1)/sum(n);
    vdat.ChirpsInBurst = 1;
    vdat.TimeStamp = mean(time);
    vdat.processing = {['Multi-chirp mean produced by ' mfilename ' using ' int2str(sum(n)) ' chirps.']};
    %vdat.Code = -10; % matlab file
    %vdat.Attenuator_1 = mean(vdat.Attenuator_1);
    %vdat.Attenuator_2 = mean(vdat.Attenuator_2);
    outputfilelist(nci) = {[outputfilenameroot '_n' sprintf('%06u',nchirpsread) '.mat']};
    eval(['save ' outputfilelist{nci} ' vdat']) % note we're just using vdat from the last burst read, but this shouldn't be used anyway??
    
    if nchirpstoget > 0
        disp(['Warning: only ' int2str(nchirpsread) ' chirps in files selected'])
    end
end
% hold on
% plot(vm)
% keyboard

fmcw_plot(outputfilelist,'maxrange',2000)