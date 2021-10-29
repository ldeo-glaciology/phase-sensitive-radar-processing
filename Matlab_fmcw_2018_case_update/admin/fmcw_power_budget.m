function fmcw_power_budget

% Calculate power and memory budget for ApRES deployment planning

% Craig Stewart
% 2014-1-5

% Deployment
Duration = 365; % days

% Settings from config file
RepSecs = 1800; % seconds
NSubBursts = 20;
nAttenuators = 1;
Average = 1;
N_ADC_SAMPLES = 60000;

% constants
Current = 0.7; %A
ChirpDuration = 1.64; % seconds
Overhead = 5; % seconds
ByteDepth = 2; % 16bit storage

% Derived
ChirpsperBurst = NSubBursts*nAttenuators;
% power
AHperChirp = Current*ChirpDuration/3600;
AHperBurst = AHperChirp*ChirpsperBurst;

% Data
BytesperChirp = ByteDepth*N_ADC_SAMPLES;
if Average
    BytesperBurst = BytesperChirp*nAttenuators;
else
    BytesperBurst = BytesperChirp*ChirpsperBurst;
end

% Totals
nBursts = (24*3600*Duration)/RepSecs;
AH = nBursts*AHperBurst;
GB = nBursts*BytesperBurst/1e9;
disp(['Power  required: ' sprintf('%4.1f',AH) ' AH'])
disp(['Memory required: ' sprintf('%4.1f',GB) ' GB'])
    
