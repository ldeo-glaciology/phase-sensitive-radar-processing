function cfg = fmcw_process_config_vsr

% Default configuration file for fmcw_process
%
% Edit this and rename for custom processing settings.

% Default config
cfg.notes = 'default config file';

% Use test data
cfg.useTestFiles = 0;
cfg.filename1 = 'c11_2013-01-25-0053.DAT';
cfg.filename2 = 'c11_2014-01-19_223232.dat';

% Pre-processing
cfg.fRange = [2e8 4e8]; % Tx frequency range to use - full range
%cfg.fRange = [2.5e8 3.5e8];  % half bandwidth
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
cfg.maxRange = 2000; % maximum bed range
cfg.winFun = @blackman;
%cfg.winFun = @blackmanharris; % less spectral leakage so better to pick up near bed internals

% Range error estimate
cfg.errorMethod = 'emperical'; % 'emperical' 'assumedNoiseFloor'
cfg.noiseFloordB = -100; % Assumed level of noise

% Find bed (or max range for vsr)
cfg.maxDepthMethod = 'config'; % manual config auto
cfg.maxDepthConfig = 1000;
cfg.bedMethod = 'maxAmp'; % 'ampThresh' 'maxAmp' 'xcor' ???
cfg.ampThreshdB = -50; % dB - minimum bed strength
cfg.bedSearchRange = [100 inf]; % bed search range (m)

% Bulk lag matching (co-registration)
cfg.doBulkAllignment = 1; %
cfg.bulkAlignRange = [60 90];
cfg.maxOffsetM = 10; % 10m recoverable offset near surface

% Chunk lag matching (depth dependant co-registration)
cfg.minDepth = 10; % to avoid breakthrough and cables etc (cables 2m in 2013, 5m in 2013).
cfg.bedBuffer = 10; % m buffer to exclude spectral leakage from bed return
cfg.coarseChunkWidth = 15; % long segments more uniquely define lag - except if there is high strain
cfg.maxStrain = 0.005; % maximum strain magnitude to search for
cfg.minAmpCor = 0.95; % Minimum amplitude correlation to use
cfg.minAmpCorProm = 0.05; % Minimum difference between max correlation and next best local maximum

% Chunk phase difference
cfg.doUseCoarseOffset = 0; % uses coarse offset determined above to specify rough lag for fine offset
cfg.doPolySmoothCoarseOffset = 1;
cfg.polyOrder = 3;
cfg.phaseDiffMethod = 'xcor'; % 'xcor' 'peakDiff'
cfg.chunkWidth = 8; % between 4 to 8 is a good compromise
cfg.doSmartUnwrap = 1; % phase difference tracking to determine bin lag -
% note SmartUnwrap only works with high pad factors - 10 works.
% high strain rates may cause problems...

% Strain fitting
cfg.minCohereCoarse = 0.6; % minimum correlation on lag offset (i.e. confidence we have the right lag)
cfg.minCohereFine = 0.6; % minimum coherence to use in strain estimate
cfg.minCoherePhase = 0.6; % phase coherence over win...
cfg.firnDepth = 70; % used to select depth range for vertical strain estimate
cfg.fitMethod = 'robust'; % 'menke' 'robust' 'regress'
cfg.fit.type = 'linear'; % 'linear' 'quadratic'
cfg.fit.errorMetric = 'standardError'; % 'standardError' 'alpha'
cfg.minPointsToFit = 3; % minimum number of depth offset to use in estimate

% Melt
cfg.doMeltEstimate = 0;
cfg.useUDStrainRate = 0; % Override strain rate estimate with used defined estimate
% This is useful when calculating melt over short periods where errors in
% the strain estimate dominate total melt errors
cfg.udvsr = 0.001; % user defined strain rate (/year) to use for melt
cfg.udvsre = 0.0001; % user defined strain rate error (/year) to use for melt

% Bed shift
cfg.bedShiftMethod = 'xcorr'; % 'xcorr' 'rangeDiff' (xcorr more reliable... 
cfg.xcorBedWin = [-2 1]; % m bed window limits - m from bed peak - to use in xcorr or rangediff
%but rangeDiff more accurate if it gets the interger amigiuity correct...)
cfg.rangeBedWin = [-1 0]; % m bed window limits - m from bed peak - to use in xcorr or rangediff
cfg.bed.maxWrapFraction = 0.25; % maximum tolerable % of wrapping present in bed window
% note wrapping fraction won't exceed 0.5 as it will just choose the other range...
%cfg.maxBedOffset = 10;

% Output
cfg.verbose = 1; % results to screen
cfg.doSaveOutput = 0; % save mat file

% Plots
cfg.doPlotAll = 0; % plot lots of other stuff

% Individual control of all plots
cfg.doPlotAlignBulk = 0;
cfg.doPlotAlignCoarse = 1;
cfg.doPlotAlignFine = 1;
cfg.doPlotFit = 1;
cfg.doPlotBedShift = 1;
cfg.doPlotMelt = 1;
cfg.doPlotResid = 1;

