function BED = fmcw_bed_info_batch

% BED = fmcw_bed_info_batch
%
% Calculate bed info for a batch of files using fmcw_bed_info
%
% Craig Stewart
% 2014/10/21

if 1
    [file,path] = uigetfile({'*.dat','*.DAT'},'multiselect','on');
    if isa(file,'char')
        file = {file};
    end
else
    % All 2013 ross shots
    file = {'a01_2013-01-16-0524.DAT','a02_2013-01-17-2029.DAT','a03_2013-01-17-2111.DAT','a04_2013-01-17-2231.DAT','a05_2013-01-17-2334.DAT','b01_2013-01-24-0312.DAT','b02_2013-01-24-0402.DAT','b03_2013-01-24-0455.DAT','b04_2013-01-24-0536.DAT','b05_2013-01-24-0624.DAT','b06_2013-01-24-2117.DAT','b07_2013-01-24-2207.DAT','b08_2013-01-24-2258.DAT','c01_2013-01-25-2333.DAT','c02_2013-01-25-2249.DAT','c03_2013-01-25-2205.DAT','c04_2013-01-25-2122.DAT','c05_2013-01-25-2039.DAT','c06_2013-01-25-0449.DAT','c07_2013-01-25-0410.DAT','c08_2013-01-25-0332.DAT','c09_2013-01-25-0254.DAT','c10_2013-01-25-0219.DAT','c11_2013-01-25-0053.DAT','d01_2013-01-26-0346.DAT','d02_2013-01-26-0440.DAT','d03_2013-01-26-0523.DAT','d04_2013-01-27-2104.DAT','d05_2013-01-27-2204.DAT','d06_2013-01-27-2246.DAT','d07_2013-01-27-2334.DAT','d08_2013-01-28-0041.DAT','d09_2013-01-28-0116.DAT','d10_2013-01-28-0150.DAT','d11_2013-01-28-0258.DAT','e07_2013-01-29-0130.DAT','e08_2013-01-29-0006.DAT','e09_2013-01-28-2258.DAT','e10_2013-01-28-2205.DAT','e11_2013-01-28-2124.DAT','e12_2013-01-28-0519.DAT'};
end

depthRange = [0 2000]; % depth search range for bed

BED = [];
doPlot = 0;
for ii = 1:length(file)
    vdat = fmcw_load(file{ii});
    bed = fmcw_bed_info(vdat,depthRange,doPlot);
    bed.file = file{ii};
    bed.time = vdat.TimeStamp
    BED = [BED; bed];
    disp([bed.file ': ' num2str(bed.phaseGrad,'%4.2f')])
end
