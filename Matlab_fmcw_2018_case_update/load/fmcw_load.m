function vdat = fmcw_load(filename,burst)

% vdat = fmcw_load(filename,burst)
%
% Load FMCW radar burst and metadata
%
% input: filename and burst number
% output: vdat - structure of metadata and data

% Craig Stewart
% 2013 April 24
% 2013 September 30 - corrected error in vif scaling
% 2014/5/20 time stamp moved here from fmcw_derive_parameters (so that this
% is not overwritted later)
% 2014/5/21 changed how radar chirp is defined (now using chirp gradient as
% fundamental parameter)
% 2014/5/22 fixed bug in chirptime
% 2014/8/21 moved make odd length to external (called from fmcw_range)
% 2014/10/22 KWN - edited to allow for new headers in RMB2 files

% Check args
if nargin == 0
    [filename, ~] = uigetfile(['*.dat;*.DAT;*.000;*.mat'],'Choose radar file to plot','multiselect','on');
    disp(filename)
    if isa(filename,'double') % no files chosen
        return
    end
end
if nargin < 2
    burst = 1;
end

% Settings
SamplesPerChirp = 40000; %SHOULD THIS BE 40001?
doShowFileType = 0;

%% Load data and reshape array
[~,~,ext] = fileparts(filename);
if strcmp(ext,'.mat')
    load(filename); % note when you load a mat file you get whatever burst was stored in this - not the one you selected

%ELIZ: NEED TO ADD 6TH FILE FORMAT -- THIS NEW ONE; HOW TO DIFFERENTIATE?
else
    FileFormat = fmcw_file_format(filename);
    if FileFormat == 6
       vdat = LoadBurstRMB6(filename, burst, SamplesPerChirp);%Data from 2018 Juneau Icefield run (any changes issued within ApRES?)
    elseif FileFormat == 5
        vdat = LoadBurstRMB5(filename, burst, SamplesPerChirp);% Data from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)
    elseif FileFormat == 4
        vdat = LoadBurstRMB4(filename, burst, SamplesPerChirp); % Data from after Oct 2013  (RMBB??)
        %disp('file format4')
    elseif FileFormat == 3
        vdat = LoadBurstRMB3(filename, burst, SamplesPerChirp); % Data from Jan 2013 (RMBA?)
    elseif FileFormat == 2
        vdat = LoadOldBurstCompat(filename, burst, SamplesPerChirp); % Data from Prototype FMCW radar (nov 2012) (RMBA?)
    end
end
if doShowFileType
    disp(['Detected file format: ' int2str(FileFormat)])
end

% Check file was found
switch(vdat.Code)
    case -1
        fprintf('Unable to open file: %s\n',filename);
        return
    case -2
        fprintf('Corrupt header in burst - fatal %d\n',vdat.Burst);
        return
    case -4
        disp(['Burst ' int2str(burst) ' not found in file ' filename]);
        return
end

% Extract just good chirp data from voltage record and rearrange into
% matrix with one chirp per row
% note: you can't just use reshape as we are also cropping the 20K samples
% of sync tone etc which occur after each 40K of chirp.
vdat.vif = zeros(vdat.ChirpsInBurst,SamplesPerChirp); % preallocate array
AttSet = vdat.Attenuator_1 + 1i*vdat.Attenuator_2; % unique code for attenuator setting

% Load each chirp into a row
chirpInterval = 1.6384/(24*3600); % days
for chirp = 1:vdat.ChirpsInBurst
    vdat.vif(chirp,:) = vdat.v(vdat.Startind(chirp):vdat.Endind(chirp));
    vdat.chirpNum(chirp,1) = chirp; % chirp number in burst
    vdat.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); % attenuator setting for chirp
    vdat.chirpTime(chirp,1) = vdat.TimeStamp + chirpInterval*(chirp-1); % time of chirp
end

%% Add metadata to structure
% Sampling parameters
vdat.SamplesPerChirp = SamplesPerChirp;
vdat.fs = 4e4; % sampling frequency
vdat.f0 = 2e8; % start frequency
%vdat.fc = 3e8; % start frequency
vdat.K = 2*pi*2e8; % chirp gradient in rad/s/s (200MHz/s)
%vdat.f0 = vdat.f0 + (vdat.K/(4*pi))/vdat.fs; % start frequency
vdat.processing = {};

% Environment settings
material = 'ice'; % or 'cable'
switch material
    case 'ice'
        vdat.er = 3.1;
    case 'cable'
        disp('****************************************************************')
        disp(' ')
        disp('WARNING USING Er = 1.5 FOR PROCESSING TEST SHOTS WITH DELAY LINE ONLY')
        disp(' ')
        disp('****************************************************************')
        disp(' ')
        vdat.er = 1.5;
end

% Derived constants
vdat = fmcw_derive_parameters(vdat);

% Create time and frequency stamp for samples
vdat.t = vdat.dt*(0:size(vdat.vif,2)-1); % sampling times (rel to first)
vdat.f = vdat.f0 + vdat.t.*vdat.K/(2*pi);

vdat.filename = filename;% Calibrate

%vdat = fmcw_burst_split_by_att(vdat);


ca13 = [1 6]; % 2013
%ca14 = [1 2]; % 2014
%ca = [1 4];
%vdat = fmcw_cal(vdat,ca13);
