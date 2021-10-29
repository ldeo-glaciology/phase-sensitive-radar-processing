function data = fmcw_plot(plotop,varargin)

% Plot fmcw shot
%
% Plot FMCW radar data
%
% args:
%
% plotop = string of plot options: containing code letters of:
%   type of plots -
%     r = raw (voltage vs sample number from start of burst whole burst)
%     t = time (voltage vs time from start of chirp for whole burst)
%     a = amplitude (processed data amplitude in range domain) (default)
%     p = phase (processed data phase in range domain)
%   display options - (override display defaults)
%     no = no overlay - don't overlay data on same axis
%     nl = no legend
%     c = plot individual chirps (rather than burst mean)
%
% vararg = parameter name-value pair list with the following possible parameters:
%
% 'filelist' (or 'fl'), followed by  = single filename (char), or file list
%   (cell array) of files to plot, or can use 'last', which will plot the
%   same files as last time
%
% 'burstlist' (or 'bl'), followed by burst numbers to plot. This will plot
%   the same bursts from each file - or can be 'all').  List of burst numbers
%   is in standard Matlab format eg [3:7] or [1:4,8].
%
% 'chirplist' (or 'cl'), followed by chirp numbers to plot. This will plot
%   the same chirps from each burst - or can be 'all').  List of chirp numbers
%   is in standard Matlab format eg [3:7] or [1:4,8].
%
% 'maxrange' (or 'mr'), followed by maximum depth to plot to, in m (default = 1500).
%
% 'maxbursts' (or 'mb'), followed by maximum number of bursts to plot, (default = 10).
%
% 'maxchirps' (or 'mc'), followed by maximum number of chirps to plot (default = 100).
%
% 'win' followed by window function handle (default @blackman)

% 'attn' followed by desired attenuation setting (1st or 2nd) (default =
% all)

% 'attnset' followed by desired attenuation setting e.g. 18 or [18,28]
%
% 'sn' followed by desired bounds for snip signal [istart, iend] 

% 'plot_on' turns plotting on and off (1, 0)
%
% Craig Stewart
% 2013-Sep-30
% 2014/5/6 - major re-factor

% Updated by Elizabeth Case
% 2019

%%
%close all

%% Settings
maxbursts = 1; % per file
maxchirps = 100; % per burst
maxrangedefault = 1000; % default range to crop data to (m) (was 800, changed to 1500 - jk, 07/24/18
%maxrangedefault=2000;
p = 2; % padfactor;
defaultplotop = 'a';

%% Check args
if nargin == 0 % plot options not specified
    % Default plot options
    plotop = defaultplotop;
% else
%     plotop = [defaultplotop plotop];
%     plotop = strrep(plotop,'aa','a');
end
%paramnames = {'filelist','burstlist','chirplist'};

if nargin > 1 % burst not specified
    for n = 1:2:length(varargin)-1
        switch varargin{n}
            case 'filelist'
                filelist = varargin{n+1};
            case 'fl'
                filelist = varargin{n+1};
            case 'burstlist'
                burstlist = varargin{n+1};
            case 'bl'
                burstlist = varargin{n+1};
            case 'chirplist'
                chirplist = varargin{n+1};
            case 'cl'
                chirplist = varargin{n+1};
            case 'maxrange'
                maxrange = varargin{n+1};
            case 'mr'
                maxrange = varargin{n+1};
            case 'maxbursts'
                maxbursts = varargin{n+1};
            case 'mb'
                maxbursts = varargin{n+1};
            case 'maxchirps'
                maxchirps = varargin{n+1};
            case 'mc'
                maxchirps = varargin{n+1};
            case 'win'
                win = varargin{n+1};
            case 'attn'
                whichattn = varargin{n+1};
            case 'attnset'
                attnset = varargin{n+1};
            case 'sn'
                sn = varargin{n+1};
            case 'plot_on'
                plot_on = varargin{n+1};
            otherwise
                disp(['Warning: param ' varargin{n} ' not recognised'])
        end
    end
end

% Filelist not defined
if ~exist('filelist','var')
    %     datafid = fopen('fmcw_data_dir.txt','r');
    %     if datafid == -1
    %         initpath = pwd;
    %     else
    %         initpath = fgets(datafid); fclose(datafid);
    %     end
    %    [filelist, pathname] = uigetfile([initpath,'*.dat;*.DAT;*.000;*.mat'],'Choose radar file to plot','multiselect','on');
    [filelist, pathname] = uigetfile(['*.dat;*.DAT;*.000;*.mat'],'Choose radar file to plot','multiselect','on');

    if isa(filelist,'double') % no files chosen
        return
    end
    datafid = fopen('fmcw_data_dir.txt','w'); fprintf(datafid,'%s',pathname); fclose(datafid);
    if ischar(filelist) % cover the case of a single file selected (where uigetfile returns a char
        filelist = {filelist};
    end
    for ii = 1:length(filelist) % Convert all filenames to path/file for consistency
        filelist(ii) = {[pathname filelist{ii}]};
    end
    % Save filelist for future use
    filelistfid = fopen('lastlist.txt','w');
    for ii = 1:length(filelist), fprintf(filelistfid,'%s\n',filelist{ii}); end
    fclose(filelistfid);
end

% Read filelist
if isa(filelist,'char') % if file list is a string this could be a datafile name or the filename of a datafilename list
    if strcmp(filelist,'last')
        filelist = 'lastlist.txt'; % special case - load last used files
    end
    [~,~,ext] = fileparts(filelist);
    if strcmp(ext,'.txt')
        [fid,msg]  = fopen(filelist,'rt');
        clear filelist
        if fid==-1
            %disp(['Trouble reading filelist file: ' filelist])
            %error(msg)
            disp('plot history not available in this dir - manually select file')
            return
        else
            ii = 0;
            while ~feof(fid)
                ii = ii+1;
                filelist(ii) = {fgetl(fid)};
            end
            fclose(fid);
        end
    else
        filelist = {filelist}; % This is just a data filename (
    end
end

% Burstlist not defined
if ~exist('burstlist','var')
    burstlist = 'all';
    %disp(['Burst number not specified: defaulting to first ' int2str(maxbursts)])
end
% chirplist not defined
if ~exist('chirplist','var')
    chirplist = 'all';
end
% maxrange not defined
if ~exist('maxrange','var')
    maxrange = maxrangedefault;
end
% win not defined
if ~exist('win','var')
    win = @blackman;
end

if ~exist('plot_on','var')
    plot_on = 1;
end

% Display comand line for this plot sequence
if 1 %nargin <2
    disp('To repeat this plot command use this command:')
    fileliststr = 'last';
    if isa(burstlist,'char')
        burstliststr = ['''' burstlist ''''];
    else
        burstliststr = mat2str(burstlist);
    end
    if isa(chirplist,'char')
        chirpliststr = ['''' chirplist ''''];
    else
        chirpliststr = mat2str(chirplist);
    end
    
    disp(['fmcw_plot(''' plotop ''',''filelist'',''' fileliststr ''',''burstlist'',' burstliststr ',''chirplist'',' chirpliststr ')'])
end

% Replace character options for burst and chirp lists with numeric
if isa(burstlist,'char')
    if strcmp(burstlist,'all')
        burstlist = 1:maxbursts; %
    else
        error(['Unrecognised burstlist option ' burstlist])
    end
end
if isa(chirplist,'char')
    if strcmp(chirplist,'all')
        chirplist = [1:maxchirps];
    end
end

% Q: Faster if all defaults just set ahead of time?
%% Parse plot settings
if isempty(strfind(plotop,'no')) % stop plot overlay
    overlay = 1; % default = overlay on
else
    overlay = 0;
end
if isempty(strfind(plotop,'nl')) % stop legend
    showlegend = 1; % default = legend on
else
    showlegend = 0;
end
if isempty(strfind(plotop,'r')) % plot raw data (multi chirp timeseries)
    plot_raw = 0;
else
    plot_raw = 1;
end
if isempty(strfind(plotop,'t')) % time domain voltage data (chirps overlayed)
    plot_time = 0;
else
    plot_time = 1;
end
if isempty(strfind(plotop,'a')) % processed amplitude data
    plot_amp = 0;
else
    plot_amp = 1;
end
if isempty(strfind(plotop,'p')) % processed phase data
    plot_phase = 0;
else
    plot_phase = 1;
end
if isempty(strfind(plotop,'c')) % plot individual chirps (rather than average all chirps in burst)
    DoAverageBurst = 1;
else
    DoAverageBurst = 0;
end

%% Loop through loading/plotting data
for FileNo = 1:length(filelist)
    [path,name,ext] = fileparts(filelist{FileNo});
    filename = [name ext];
    if strcmp(ext,'.mat')
        burstlist = 1;
        %disp('Loading single burst from mat file')
    end
    
    getBurst = 1;
    BurstNo = 0;
    while BurstNo<length(burstlist) && getBurst
        BurstNo = BurstNo + 1;
        thisburst = burstlist(BurstNo); 
        
        % Load and process data
        
        vdat = fmcw_load(filelist{FileNo},thisburst); % load data
      
        
        if vdat.Code == -4 % burst not found in file
            %disp(['Only ' int2str(thisburst-1) ' bursts in file ' filename])
            %return
            getBurst = 0;
        elseif vdat.Code == -5
            disp(['No chirp starts found in file ' filename ' burst ' int2str(thisburst) ': - Corrupted data?'])
            %getBurst = 0;
        else
            % we have a good data file
            
   
            % Crop chirplist to a max of ChirpsInBurst
            BurstChirpList = chirplist(chirplist<=vdat.ChirpsInBurst);
            if ~isempty(BurstChirpList)
                

                % Crop data to these chirps only
                vdat = fmcw_burst_subset(vdat,BurstChirpList);
                
                % Snips data if bounds are given
                if exist('sn','var')
                    vdat = fmcw_burst_snip_signal(vdat, sn(1), sn(2));
                    snip = ['snipped to ' int2str(sn(1)) ':' int2str(sn(2))];
                end
                
                if ~exist('snip','var')
                    snip = [''];
                end
                


                % Split burst into various attenuator settings
                vdats = fmcw_burst_split_by_att(vdat);
                
                % Chooses attenuator setting
                if exist('attnset','var')
                    [col,whichattn] = find([vdats.Attenuator_1] == attnset);
                end
                
                if ~exist('whichattn','var')
                    whichattn = 1:1:vdat.NAttenuators;
                end
                
                vdats = vdats(whichattn);
                
                
                for asn = 1:length(vdats) % attenuator setting number
                    % Generate labels for each chirp to plot
                    shotNamePrefix = [strrep(name,'_','-') ' b:' int2str(thisburst) ' c: '];
                    if DoAverageBurst
                        vdat = fmcw_burst_mean(vdats(asn));
                        chirpname = {[shotNamePrefix ' avg' int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB ' snip]};
                    else
                        vdat = vdats(asn);
                        for ci = 1:size(vdat.vif,1)
                            chirpname(ci) = {[shotNamePrefix int2str(vdat.chirpNum(ci))] snip};
                        end
                    end
                    
                    % phase process data
                    [rc,rf,sr,su] = fmcw_range(vdat,p,maxrange,win);
                    
                    data.vdat = vdats(asn);
                    data.rc = rc;
                    data.amp = 20*log10(abs(su));
                            
                    %% Make plots
                    
                    if plot_on == 1
                        % Figure 1: all data in single timeseries
                        if plot_raw
                            figure; % open a new figure for each burst
                            hold on
                            plot(vdat.v);
                            % Mark the start of each chirp with a red line
                            m = repmat(vdat.Startind',2,1);% m = m'; m = m(:);
                            y = repmat([0; 2.5],1,vdat.ChirpsInBurst);
                            plot(m,y,'r');
                            title(['Raw timeseries: file: ' filename ' burst:' int2str(thisburst)],'interpreter','none')
                            xlabel('sample num')
                            ylabel('voltage (V)')
                        end

                        % Figure 2: time domain of shots stacked on each other
                        if plot_time
                            [tfig,tax] = open_tfig(overlay,vdat,asn);
                            for ci = 1:size(vdat.vif,1)
                                h_sig = plot(vdat.t,vdat.vif(ci,:),'tag',chirpname{ci},'userdata','fmcw'); % signal
                            end
                            updateaxis(tax)
                        end

                        % Figure 3: amplitude profile
                        if plot_amp
                            [afig,aax] = open_afig(overlay,asn);
                            for ci = 1:size(vdat.vif,1)
                                
                                tc = rc./vdat.ci;
                                h_sig = plot(20*log10(abs(su(ci,:))),tc,'tag',chirpname{ci},'userdata','fmcw');
                            end
                            updateaxis(aax)
                        end

                        % Figure 4: phase profile
                        if plot_phase
                            [pfig,pax] = open_pfig(overlay,asn);
                            for ci = 1:size(vdat.vif,1)
                                h_sig = plot(rc,angle(su(ci,:)),'tag',chirpname{ci},'userdata','fmcw');
                            end
                            updateaxis(pax)
                        end
                    end
                end % ends phase profile
            else
                disp('Chirp not found in burst')
            end
        end
    end
end
% Link x axes
if exist('aax','var') && exist('pax','var')
    linkaxes([aax pax],'x')
end


%--------------------------------------------------------------------------
function [tfig,tax] = open_tfig(overlay,vdat,nattn)
% Find or open figure

%new figure for each attn setting for overlay
    tagf = sprintf('tfig%d',nattn);
    taga = sprintf('tax%d',nattn);
    
if overlay
    % Find previous tfig figures

    tfig = findobj('tag',(tagf));
    tax = findobj('tag',(taga));
    if isempty(tfig)
        clear tfig tax
    end
end
if ~exist('tfig','var')
   
    tfig = figure('units','normalized','tag',(tagf),'userdata','fmcw'); % ,'position',[0.05,0.1,0.43,0.6]
    tax = axes('tag',(taga));
    title('FMCW time domain stacked chirps')
    hold on
    box on
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    ylim([-0.25 2.75])
    h_sat = plot(vdat.t,repmat(0,size(vdat.t)),'r','HandleVisibility','off'); % ADC saturation level
    h_sat = plot(vdat.t,repmat(2.5,size(vdat.t)),'r','HandleVisibility','off'); % ADC saturation level
end
figure(tfig(1))
hold on
%xlim([0,maxrange]);

%--------------------------------------------------------------------------
function [afig,aax] = open_afig(overlay,nattn)
% Find or open figure

    tagf = sprintf('afig%d',nattn);
    taga = sprintf('aax%d',nattn);
    
if overlay

    % Find previous tfig figures
    afig = findobj('tag',(tagf));
    aax = findobj('tag',(taga));
    if isempty(afig)
        clear afig aax
    end
end
if ~exist('afig','var')
      
    afig = figure('units','normalized','position',[0.1,0.5,0.5,0.4],'tag',(tagf),'userdata','fmcw'); %
    aax = axes('tag',(taga));
    title('FMCW range domain amplitude')
    box on
    xlabel('Range (m)');
    ylabel('amplitude (dB Vrms)')
end
figure(afig(1));
hold on
%xlim([0,maxrange]);

%--------------------------------------------------------------------------
function [pfig,pax] = open_pfig(overlay,nattn)
% Find or open figure

    tagf = sprintf('pfig%d',nattn);
    taga = sprintf('pax%d',nattn);
    
if overlay
    % Find previous tfig figures
    pfig = findobj('tag',(tagf));
    pax = findobj('tag',(taga));
    if isempty(pfig)
        clear pfig pax
    end
end
if ~exist('pfig','var')
    pfig = figure('units','normalized','position',[0.1,0,0.5,0.4],'tag',(tagf),'userdata','fmcw'); % ,'position',,[0.1,0.5,0.8,0.5]
    pax = axes('tag',(taga));
    title('FMCW range domain phase')
    box on
    xlabel('Range (m)');
    ylabel('Phase (rad)')
end
figure(pfig(1));
xlim([0 100]);
hold on
%xlim([0,maxrange]);

%--------------------------------------------------------------------------
function updateaxis(ax)
h = findobj('parent',ax,'userdata','fmcw');
if ~isempty(h)
% Update line colors
n = length(h);
if n>1
    for ii = 1:n
        set(h(ii),'col',getcol([0 n+1],ii,jet))
    end
end

% Update legend
shotnames = get(h,'tag');
l = legend(h,shotnames);
set(l, 'Interpreter', 'none')
end

%--------------------------------------------------------------------------
function col = getcol(clims,val,cmap)
% function col = getcol(clims,val,cmap)
%
% Gets a RGB colour vector for value "val" given color limits
% "clims" and colormap "cmap" (default jet)
%
% Craig Stewart 23-Aug-2007

if nargin == 2
    cmap = jet;
end
% reshape val to a column vector
val = val(:);

% Check that clim is a 2 element vector
if numel(clims)~=2
    disp(['numel(clims) = ' int2str(numel(clims))])
    error('clims must have 2 elements')
end

% Clip out of bound values to the colorlimits
val = min(val,repmat(max(clims),size(val)));
val = max(val,repmat(min(clims),size(val)));
% Normalise the data to the colorlimits
norm_val = interp1(clims,[1 size(cmap,1)],val);

% Interpolate to get the color
col = interp2([1:size(cmap,2)],[1:size(cmap,1)],cmap,...
    repmat([1:size(cmap,2)],numel(val),1),repmat(norm_val,1,size(cmap,2)));
%--------------------------------------------------------------------------






















%**************************************************************************
% JUNK here


% plot_profile
%                 if plot_profile
%                     clear pfig
%                     if overlay
%                         % Find previous tfig figures
%                         profig = findobj('tag','profig');
%                         if isempty(profig)
%                             clear profig
%                         end
%                     end
%                     if ~exist('profig','var')
%                         pfig = figure('units','normalized','position',[0.1,0,0.8,0.4],'tag','profig','userdata'); % ,'position',,[0.1,0.5,0.8,0.5]
%                         title('FMCW profile')
%                         box on
%                         set(gca,'ydir','rev')
%                         xlabel('Chirp #');
%                         ylabel('Range (m)')
%                         hold on
%                     end
%                     figure(profig(1));
%                     ylim([0,maxrange]);
%
%                     chirpno = 1:;
%                 end