function fmcw_plot(filelist,varargin)

% Plot fmcw shot
%
% Plot FMCW radar data
%
% args:
%
% filelist: list of files to plot. Single filename (char), or file list 
% (cell array) or 'last', which will plot the same files as last time, or
% '[]', which will allow file selection from browser.
%
% vararg: = parameter name-value pair list with the following possible parameters:
%
% 'plotop' (or 'po') followed by string of plot options: containing code letters:
%   type of plots -
%     r = raw (voltage vs sample number from start of burst whole burst)
%     t = time (voltage vs time from start of chirp for whole burst)
%     a = amplitude (and phase) (processed data amplitude in range domain) (default)
%
%   display options - (override display defaults)
%     no = no overlay - don't overlay data on same axis
%     nl = no legend
%     c = plot individual chirps (rather than burst mean)
%
% 'burstlist' (or 'bl'), followed by burst numbers to plot. This will plot
%   the same bursts from each file - or can be 'all').  List of burst numbers
%   is in standard Matlab format eg [3:7] or [1:4,8].
%
% 'chirplist' (or 'cl'), followed by chirp numbers to plot. This will plot
%   the same chirps from each burst - or can be 'all').  List of chirp numbers
%   is in standard Matlab format eg [3:7] or [1:4,8].
%
% 'maxrange' (or 'mr'), followed by maximum depth to plot to, in m (default = 3500).
%
% 'maxbursts' (or 'mb'), followed by maximum number of bursts to plot, (default = 10).
%
% 'maxchirps' (or 'mc'), followed by maximum number of chirps to plot (default = 10).
%
% 'color' (or 'c'), followed by line color (default = chirps maped by colormap)
%
% 'label' (or 'l'), follwed by text to add to legend labels (after auto text)
%
% 'win' followed by window function handle (default @blackman)
%
% Craig Stewart
% 2013-Sep-30
% 2014/5/6 - major re-factor

%% Settings
maxbursts = 10; % per file
maxchirps = 10; % per burst
maxrangedefault = 3500; % default range to crop data to (m)
p = 2; % padfactor;

%% Check args
if nargin == 0 % plot options not specified
    filelist = [];
end
%paramnames = {'filelist','burstlist','chirplist'};

if nargin > 1 % burst not specified
    for n = 1:2:length(varargin)-1
        switch varargin{n}
            case 'plotop'
                plotop = varargin{n+1};
            case 'po'
                plotop = varargin{n+1};
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
            case 'color'
                color = varargin{n+1};
            case 'c'
                color = varargin{n+1};
            case 'label'
                label = varargin{n+1};
            case 'l'
                label = varargin{n+1};
            case 'win'
                win = varargin{n+1};
            otherwise
                disp(['Warning: param ' varargin{n} ' not recognised'])
        end
    end
end

% Filelist not defined
if isempty(filelist)
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

%% Set default behaviour
% Default plot options
if ~exist('plotop','var')
    plotop = 'a';
    %disp(['Burst number not specified: defaulting to first ' int2str(maxbursts)])
end
% Check there is a plot style specified
if isempty([strfind(plotop,'r') strfind(plotop,'t') strfind(plotop,'a')])
    plotop = ['a' plotop];
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
if exist('color','var')
    colorMode = 'manual';
else
    colorMode = 'auto';
    color = 'k';
end
% line label not defined
if ~exist('label','var')
    label = '';
end
% win not defined
if ~exist('win','var')
    win = @blackman;
end

% Display command line for this plot sequence
if 1; %nargin <2
    fileliststr = '{';
    for ii = 1:length(filelist)
        fileliststr = [fileliststr '''' filelist{ii} ''',']; %#ok<AGROW>
    end
    fileliststr = [fileliststr(1:end-1) '}'];
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
    cmd = ['fmcw_plot(' fileliststr ',''plotop'',''' plotop ''',''burstlist'',' burstliststr ',''chirplist'',' chirpliststr ',''maxrange'',' mat2str(maxrange)];
    if strcmp(colorMode,'manual')
        cmd = [cmd ',''color'',' color];
    end
    cmd = [cmd ')'];
    
    clipboard('copy',cmd)
    disp('To repeat this plot command use this command (paste from clipboard):')
    disp(cmd)
    disp('')
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
    plot_amp_phase = 0;
else
    plot_amp_phase = 1;
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
        fileburstlist = 1;
    else
        fileburstlist = burstlist;
        %disp('Loading single burst from mat file')
    end
    
    getBurst = 1;
    BurstNo = 0;
    while BurstNo<length(fileburstlist) && getBurst
        BurstNo = BurstNo + 1;
        thisburst = fileburstlist(BurstNo);
        
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
            
            %% Plot Figure 1: all data in single timeseries
            if plot_raw
                figure; % open a new figure for each burst
                rax = axes;
                sn = 1:length(vdat.v);
                plot(vdat.v,'k');
                % Overlay chirps as used
                hold on
                for ci = 1:length(vdat.Startind)
                    %h = plot(sn(vdat.Startind(ci):vdat.Endind(ci)),vdat.v(vdat.Startind(ci):vdat.Endind(ci)),'.','color',color,'userData','fmcw');
                    h = plot(sn(vdat.Startind(ci):vdat.Endind(ci)),vdat.v(vdat.Startind(ci):vdat.Endind(ci)),'userData','fmcw');
                    %setappdata(h,'colorMode',colorMode)
                end
                
                % Mark the start of each chirp with a red line
                m = [vdat.Startind(:)'; vdat.Startind(:)'];% m = m'; m = m(:);
                y = repmat([0; 2.5],1,length(vdat.Startind));
                plot(m,y,'r');
                title(['Raw timeseries: file: ' filename ' burst:' int2str(thisburst)],'interpreter','none')
                xlabel('sample num')
                ylabel('voltage (V)')
            end
            
            % Crop chirplist to a max of ChirpsInBurst
            BurstChirpList = chirplist(chirplist<=vdat.ChirpsInBurst);
            if ~isempty(BurstChirpList)
                
                % Crop data to these chirps only
                vdat = fmcw_burst_subset(vdat,BurstChirpList);
                %vdat = vdats;
                
                % Split burst into various attenuator settings
                vdats = fmcw_burst_split_by_att(vdat);
                for AttSetNum = 1:length(vdats) % attenuator setting number
                    % Generate labels for each chirp to plot
                    shotNamePrefix = [strrep(name,'_','-') ' b:' int2str(thisburst) ' c:'];
                    if DoAverageBurst
                        vdat = fmcw_burst_mean(vdats(AttSetNum));
                        chirpname = {[shotNamePrefix ' avg' int2str(real(vdat.chirpAtt)) '+' int2str(imag(vdat.chirpAtt)) 'dB ' label]};
                    else
                        vdat = vdats(AttSetNum);
                        for ci = 1:size(vdat.vif,1)
                            chirpname(ci) = {[shotNamePrefix sprintf('%03d',vdat.chirpNum(ci)) ' ' label]};
                        end
                    end
                    
                    % phase process data
                    [rc,~,~,su] = fmcw_range(vdat,p,maxrange,win);
                    
                    %% Make plots
                    % Figure 2: time domain of shots stacked on each other
                    if plot_time
                        if BurstNo==1 && AttSetNum ==1
                            [tax] = open_tax(overlay,vdat); % Open new figure only once
                        else
                            [tax] = open_tax(1,vdat); % overlay on
                        end
                        axes(tax), hold on
                        for ci = 1:size(vdat.vif,1)
                            h = plot(vdat.t,vdat.vif(ci,:),'tag',chirpname{ci},'color',color,'userData','fmcw'); % signal
                            setappdata(h,'colorMode',colorMode)
                        end
                    end
                    
                    % Figure 3: amplitude/phase profile
                    if plot_amp_phase
                        if BurstNo==1 && AttSetNum ==1
                            [aax,pax] = open_apax(overlay); % Open new figure only once
                        else
                            [aax,pax] = open_apax(1); % overlay on
                        end
                        axes(pax), hold on
                        for ci = 1:size(vdat.vif,1)
                            h = plot(rc,angle(su(ci,:)),'tag',chirpname{ci},'color',color,'userData','fmcw');
                            setappdata(h,'colorMode',colorMode)
                        end
                        axes(aax), hold on
                        for ci = 1:size(vdat.vif,1)
                            h = plot(rc,20*log10(abs(su(ci,:))),'tag',chirpname{ci},'color',color,'userData','fmcw');
                            setappdata(h,'colorMode',colorMode)
                        end
                    end
                    
                end 
            else
                disp('Chirp not found in burst')
            end
        end
    end
end

% Update line colors
if exist('rax','var'), refreshLineCol(rax), end % refreshLegend(rax), note can't do legend for this as we need new tags first
if exist('tax','var'), refreshLineCol(tax), refreshLegend(tax), end
if exist('aax','var'), refreshLineCol(aax), refreshLegend(aax), refreshLineCol(pax), linkaxes([aax pax],'x'), end


%--------------------------------------------------------------------------
function tax = open_tax(overlay,vdat)
% Find or open figure
if overlay
    % Find previous tfig figures
    tfig = findobj('tag','tfig');
    if ~isempty(tfig)
        tax = findobj('parent',tfig,'tag','tax');
    else
        clear tfig
    end
end
if ~exist('tfig','var')
    figure('units','normalized','tag','tfig'); % ,'position',[0.05,0.1,0.43,0.6] ,'userdata','fmcw' % tfig = 
    tax = axes('tag','tax');
    title('FMCW time domain stacked chirps')
    hold on
    box on
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    ylim([-0.25 2.75])
    plot(vdat.t,repmat(0,size(vdat.t)),'r','HandleVisibility','off'); % ADC saturation level
    plot(vdat.t,repmat(2.5,size(vdat.t)),'r','HandleVisibility','off'); % ADC saturation level
end
%xlim([0,maxrange]);

%--------------------------------------------------------------------------
function [aax,pax] = open_apax(overlay)
% Find or open figure
if overlay
    % Find previous apfig figures
    apfig = findobj('tag','apfig');
    if ~isempty(apfig)
        apfig = apfig(1); % defaults to plot into the last open fig if we have many
        aax = findobj('parent',apfig,'tag','aax');
        pax = findobj('parent',apfig,'tag','pax');
    else
        clear apfig
    end
end
if ~exist('apfig','var')
    figure('units','normalized','position',[0.5,0.1,0.5,0.8],'tag','apfig'); % apfig = 
    
    % Amp subplot
    aax = subplottight(2,1,1);
    set(aax,'tag','aax','xticklabel',[])
    %title('FMCW range domain amplitude')
    box on
    %xlabel('Range (m)');
    ylabel('amplitude (dB Vrms)')
    
    % Phase subplot
    pax = subplottight(2,1,2);
    set(pax,'tag','pax');
    box on
    xlabel('Range (m)');
    %ylim([-3.5 3.5])
    %set(pax,'YLimMode','auto') % to allow rescale
    ylabel('Phase (rad)')
end
%xlim([0,maxrange]);

%--------------------------------------------------------------------------
function refreshLineCol(ax)
if ishandle(ax)
    h = findobj('parent',ax,'userData','fmcw');
    % Find only lines with auto color
    colorModeIsAuto = false(size(h));
    for ii = 1:length(h)
        colorModeIsAuto(ii) = strcmp(getappdata(h(ii),'colorMode'),'auto');
    end
    h = h(colorModeIsAuto);
    % Update line colors
    if ~isempty(h)
        n = length(h);
        if n>1
            for ii = 1:n
                set(h(ii),'col',getcol([0 n+1],ii,jet))
            end
        end
    end
end

%--------------------------------------------------------------------------
function refreshLegend(ax)
if ishandle(ax)
    h = findobj('parent',ax,'userData','fmcw');
    
    % Update legend
    if ~isempty(h)
        shotnames = get(h,'tag');
        l = legend(h,shotnames);
        set(l, 'Interpreter', 'none')
    end
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