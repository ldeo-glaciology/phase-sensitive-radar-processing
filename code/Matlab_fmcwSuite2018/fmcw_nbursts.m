function nbursts = fmcw_nbursts(filename)

% nbursts = fmcw_nbursts(filename)
% 
% Determine the number of bursts within a file by header keyword count
% - only fast for unix at this stage...
%
% Craig Stewart
% 2014/3/25

% Check file format to determine what to search for in the header
fmt = fmcw_file_format(filename);
if fmt >= 5;
    txt = 'SW_Issue='; % Data from RMB2 after Oct 2014
elseif fmt == 4;
    txt = 'SubBursts'; % Data from after Oct 2013
elseif fmt == 3;
    txt = 'Burst'; % Data from Jan 2013
elseif fmt == 2;
    txt = 'RADAR'; % Data from Prototype FMCW radar (nov 2012)
else
    %fmt = 0; % unknown file format
    error('Unknown file format - check file')
end

% Count instance of txt
if isunix
    % Unix method
    %[status,output] = system(['tr -cs ''A-Za-z'' ''\n'' < ' filename ' | grep -c "' txt '"'],'-echo');
    % If we are given a file only determine the full path from matlab which
    [p,f,x] = fileparts(filename);
    if isempty(p)
        filename = which(filename);
    end
    cmd = ['tr -cs ''A-Za-z'' ''\n'' < ' filename ' | grep -c "' txt '"'];
    [status,output] = system(cmd);
    if status==0
        nbursts = str2num(output);
    else
        error(output)
    end
else
    fid = fopen(filename,'r');
    C = textscan(fid,'%c'); % read whole file into character array
    t = (transpose(C{:}));
    nbursts = numel(findstr(t,txt));
    fclose(fid);
end
