function cfg = fmcw_process_config_juneau(filename2)

% Default configuration file for fmcw_process
% filename1 is baseline (e.g. first measurement)
% accepts filename2 as comparison file; allows for iterative processing
%
% Edit this and rename for custom processing settings.

% Default config
cfg.notes = 'juneau config file';

% Use test data
cfg.useTestFiles = 0;
%cfg.filename1 = '/Users/elizabeth/Google Drive/home_research/projects/juneau_icefield/divide/continuous/DIR2018-07-27-1332/DATA2018-07-28-0328.DAT';
%cfg.filename2 = '/Users/elizabeth/Google Drive/home_research/projects/juneau_icefield/divide/continuous/DIR2018-07-27-1332/DATA2018-07-28-0144.DAT';
%cfg.filename2 = '/Users/elizabeth/Google Drive/home_research/projects/juneau_icefield/divide/continuous/DIR2018-07-27-1332/DATA2018-07-27-2305.DAT';
%cfg.filename1 = '/Users/elizabeth/Google Drive/home_research/projects/juneau_icefield/divide/continuous/DIR2018-07-27-1332/DATA2018-07-27-2011.DAT';
%cfg.filename2 = '/Users/elizabeth/Google Drive/home_research/projects/juneau_icefield/divide/continuous/DIR2018-07-27-1332/DATA2018-07-27-1626.DAT';
%cfg.filename1 = '/Users/elizabeth/Google Drive/home_research/projects/ApRES/data/Ronne_ApRES_KN/Site5cApRESdata 3/DATA2015-01-10-1448.DAT';
%cfg.filename1 = '/Users/elizabeth/Google Drive/home_research/projects/juneau_icefield/2019/apres/Jonny093/DIR2018-07-25-2237/DATA2018-07-25-2311.DAT';
 
%cfg.filename2 = filename2;

% Pre-processing
cfg.fRange = [2e8 4e8]; % Tx frequency range to use - full range
%cfg.fRange = [2.5e8 3.5e8];  % half bandwidth
%cfg.fRange = [2e8 3e8];  % half bandwidth
%cfg.fRange = [3e8 4e8];  % half bandwidth
cfg.doManualChirpSelect = 0;
cfg.fchirplist = [1:100];
cfg.gchirplist = [57:100];
cfg.doClean = 0;
cfg.noisePowerLimit = 0.0015; % allowable percent differrence in power from quadtratic fit to chirpnum-power trend
cfg.nNoisest = 10;

% phase-processing
cfg.p = 10; % padfactor (interpolation factor gi generating profile)
% note: higher pad factor increases the likelyhoog of estimating the
% correct bin  lag.
cfg.maxRange = 600; % maximum bed range
cfg.winFun = @blackman;
%cfg.winFun = @blackmanharris; % less spectral leakage so better to pick up near bed internals

% Range error estimate
cfg.errorMethod = 'emperical'; % 'emperical' 'assumedNoiseFloor'
cfg.noiseFloordB = -100; %-100; % Assumed level of noise

% Find bed
cfg.maxDepthMethod = 'manual'; % manual config auto
cfg.maxDepthConfig = 600; %420
cfg.bedMethod = 'maxAmp'; % 'ampThresh' 'maxAmp' 'xcor' ???
cfg.ampThreshdB = -90; %-50; % dB - minimum bed strength
cfg.bedSearchRange = [200 inf]; % bed search range (m)

% Bulk lag matching (co-registration)
cfg.doBulkAllignment = 1; %
cfg.bulkAlignRange = [40 100];
cfg.maxOffsetM = 10; % 10m recoverable offset near surface

% Chunk lag matching (depth dependant co-registration)
cfg.minDepth = 15; %10; % to avoid breakthrough and cables etc (cables 2m in 2013, 5m in 2013).
cfg.bedBuffer = 10; % m buffer to exclude spectral leakage from bed return
cfg.coarseChunkWidth = 15; % long segments more uniquely define lag - except if there is high strain
cfg.maxStrain = 0.04; % maximum strain magnitude to search for
cfg.minAmpCor = 0.93; % Minimum amplitude correlation to use
cfg.minAmpCorProm = 0.05; % Minimum difference between max correlation and next best local maximum

% Chunk phase difference
cfg.doUseCoarseOffset = 1; % uses coarse offset determined above to specify rough lag for fine offset
cfg.doPolySmoothCoarseOffset = 1;
cfg.polyOrder = 2;
cfg.phaseDiffMethod = 'xcor'; % 'xcor' 'peakDiff'
cfg.chunkWidth = 4; % between 4 to 8 is a good compromise
cfg.doSmartUnwrap = 0; % phase difference tracking to determine bin lag -
% note SmartUnwrap only works with high pad factors - 10 works.
% high strain rates may cause problems...

% Strain fitting
cfg.minCohereCoarse = 0.9; % minimum correlation on lag offset (i.e. confidence we have the right lag)
cfg.minCohereFine = 0.9; % minimum coherence to use in strain estimate
cfg.minCoherePhase = 0.9; % phase coherence over win...
cfg.firnDepth = 25; % default = 70; used to select depth range for vertical strain estimate; ELIZ changed to 25 for juneau
cfg.fitMethod = 'robust'; % 'menke' 'robust' 'regress'
cfg.fit.type = 'linear'; % 'linear' 'quadratic'
cfg.fit.errorMetric = 'standardError'; % 'standardError' 'alpha'
cfg.minPointsToFit = 3; % minimum number of depth offset to use in estimate

% Bed shift
cfg.bedShiftMethod = 'xcorr'; % 'xcorr' 'rangeDiff' (xcorr more reliable... 
cfg.xcorBedWin = [-2 1]; % m bed window limits - m from bed peak - to use in xcorr or rangediff
%but rangeDiff more accurate if it gets the interger amigiuity correct...)
cfg.rangeBedWin = [-1 0]; % m bed window limits - m from bed peak - to use in xcorr or rangediff

cfg.bed.maxWrapFraction = 0.25; % maximum tolerable % of wrapping present in bed window
% note wrapping fraction wont exceed 0.5 as it will just choose the other range...
%cfg.maxBedOffset = 10;

% Melt
cfg.doMeltEstimate = 0;
cfg.useUDStrainRate = 0; % Override strain rate estimate with used defined estimate
% This is useful when calculating melt over short periods where errors in
% the strain estimate dominate total melt errors
cfg.udvsr = 0.001; % user defined strain rate (/year) to use for melt
cfg.udvsre = 0.0001; % user defined strain rate error (/year) to use for melt

% Output
cfg.verbose = 1; % results to screen
cfg.doSaveOutput = 0; % save mat file

% Plots
cfg.doPlotAll = 1; % plot lots of other stuff

% Individual control of all plots
cfg.doPlotAlignBulk = 1;
cfg.doPlotAlignCoarse = 1;
cfg.doPlotAlignFine = 1;
cfg.doPlotFit = 1;
cfg.doPlotBedShift = 1;
cfg.doPlotMelt = 1;
cfg.doPlotResid = 1;

