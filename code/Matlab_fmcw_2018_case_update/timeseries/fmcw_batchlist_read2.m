function [config,file] = fmcw_batchlist_read2(filename)

% Reads batchlist config text file used with fmcw_timeseries
%
% Craig Stewart
% 2013-10-01

fid = fopen(filename,'rt');
hdr1 = fgetl(fid);
% check batchlist format
ver = strtok(hdr1,':');
required_version = 'fmcw_timeseries_config_v2';
if ~strcmp(ver,required_version)
    fclose(fid);
    error(['Incorrect config file format. Expecting file type: ' required_version '. Input file type: ' ver])
end

config.notes = fgetl(fid);
vsrtxt = fgetl(fid);
config.vsr = sscanf(vsrtxt,'%*s%f');
tline = fgetl(fid);
config.intRange = sscanf(tline,'%*s%f%f');
tline = fgetl(fid);
config.bedSearchRange = sscanf(tline,'%*s%f%f');

% Loop reading all filenames and burst numbers
ii = 0;
while ~feof(fid)
    ii = ii+1;
    tline = fgetl(fid);
    [datafile,remain] = strtok(tline);
    file(ii).name = {datafile};
    if isempty(remain)
        file(ii).burst = nan;
    else
        file(ii).burst = str2num(remain);
    end
end