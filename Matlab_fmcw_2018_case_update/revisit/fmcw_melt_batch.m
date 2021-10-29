function fmcw_melt_batch(batchfile)

% FMCW radar - batch process sites for melt and strain using 
% config file batchfile listing file pairs
%
%
% Craig Stewart
% 2014/6/17

% Settings
doManualCheck = 1;

if nargin == 0
    batchfile = uigetfile(['*.cfg'],'Choose batchlist for files to process');
    if isa(batchfile,'double') % no files chosen
        return
    end
end

% Read batch config file: this includes a config file and list of file pairs
% note all general config settings should come before the file list
[fid,mes] = fopen(batchfile);
if fid==-1
    error(mes)
end
ii = 0;
cfg = [];
while ~feof(fid)
    s = fgetl(fid);
    if ~isempty(s) && ~strcmp(s(1),'%') && length(s)>3
        if strncmp(s,'cfg.',4)
            % We have a config line
            try
                eval(s)
            catch ME
                fclose(fid);
                disp(['Error evaluating config setting in: ' batchfile])
                disp(['line: ' s])
                keyboard
                rethrow(ME)
            end
        else
            % We have a file pair definition
            ii = ii+1;
            if ii == 1
                cfgDefault = cfg; % save master config settings
                clear cfg % as cfg will be reused for local (file pair specific) config
            end
            [f1, remain] = strtok(s);
            file1(ii) = {f1};
            [f2, cfgdef] = strtok(remain);
            file2(ii) = {f2};
            % The rest of the line is for custom config settings:
            % this is a series of executable matlab commands defining the config file
            % e.g. cfg = cfgfilename; cfg.notes = 'demo'
            while numel(cfgdef)>0 && strcmp(cfgdef(1),' ')
                cfgdef = cfgdef(2:end); % remove leading spaces
            end
            if isempty(cfgdef) || strcmp(cfgdef(1),'%')
                CFG(ii) = {cfgDefault};
            else
                try
                    eval(cfgdef)
                catch ME
                    fclose(fid);
                    disp(['Error evaluating config setting in: ' batchfile])
                    disp(['line: ' s])
                    disp(['config: ' cfgdef])
                    keyboard
                    rethrow(ME)
                end
                CFG(ii) = {update_cfg(cfgDefault,cfg)};
                clear cfg
            end
            %CFG(ii)
        end
    end
end
fclose(fid);

%C = textscan(fid,'%s%s');
%[file1,file2] = deal(C{:});

%% Main processing loop
for ii = 1:length(file1)
    % Check whether we have already processed this file (by looking for
    % output)
    [~,name1,~] = fileparts(file1{ii});
    [~,name2,~] = fileparts(file2{ii});
    outfile(ii) = {[name1 name2 '_meltsite.mat']};
    
    if exist(outfile{ii},'file')
        disp(['Skipping files ' file1{ii} ' with ' file2{ii}])
        disp(' - already processed')
        disp(' ')
    else           
        disp(['Processing file: ' file1{ii} ' with ' file2{ii} ' (' int2str(ii) '/' int2str(length(file1)) ')'])
        site = fmcw_melt(file1{ii},file2{ii},CFG{ii});
        
        if doManualCheck
            disp('Results:')
            disp(site)
            disp(' ')
            doSaveResults = input('Save output? (y)/n: ','s');
            if isempty(doSaveResults)
                doSaveResults = 'y';
            end
            if strcmp(doSaveResults,'y')
                % Save data
                save(outfile{ii},'site');
            end
        else
            save(outfile{ii},'site');
        end
        close all
        %clc
    end
end

fmcw_melt_batch_join(batchfile,outfile)


function cfg = update_cfg(cfg,cfgnew)

% Overwrite cfg settings with those present in cfg2
% overwrite defaults
fieldnames = fields(cfgnew);
for ii = 1:length(fieldnames);
    thisfield = fieldnames{ii};
    cfg = setfield(cfg,thisfield,getfield(cfgnew,thisfield));
end
