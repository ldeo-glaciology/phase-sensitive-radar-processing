function fmcw_help(pathname) %,mode)

% fmcw_help(mode)
%
% Get help on fmcw processing scripts

if nargin == 0
    pathname = pwd;
end
if nargin <2
    mode = 's'; %short
end

% List m files
d = dir([pathname filesep '*.m']);
for fn = 1:length(d)
    fid = fopen(d(fn).name,'rt');
    if fid ~= -1
        try
            for ii = 1:5
                tline = fgetl(fid);
            end
        catch ME
            fclose(fid);
            rethrow(ME)
        end
        fclose(fid);
        disp([d(fn).name ':' tline(2:end)])
    end
end