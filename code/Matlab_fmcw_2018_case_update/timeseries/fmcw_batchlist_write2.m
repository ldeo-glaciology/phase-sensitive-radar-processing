function filename = fmcw_batchlist_write2

% Interactively create list of files/bursts to process with
% fmwc_ts_range3.
%
% note: this is modified from fmcw_batchlist_write to provide ranges for
% internals and bed rather than single values
%
% Craig Stewart
% 2013-6-27

filename = input('Enter filename for batchlist: ','s');
[~,~,ext] = fileparts(filename);
if ~strcmp(ext,'.cfg')
    filename = [filename '.cfg'];
end
notes = input('Enter header notes for file: ','s');
vsr = input('Enter site strain rate (/year): ');

%% Select files
%getmorefiles = 1;
%while getmorefiles
datafilename = uigetfile({'*.dat;*.DAT;*.000','Radar files: .dat, .DAT, and .000'},'Choose radar files to calculate meltrate','multiselect','on');
if ischar(datafilename)
    datafilename = {datafilename};
end
    
%     if isempty(filename)
%         getmorefiles = 0;
%     else
        %         for ii = 1:length(filename)
        %             burst = input(['Enter burst num for shot ' int2str(ii) ' (file:' shot(ii).filename '): ']);
        %         end
        %     else
%    end
%end
    
%% Plot the shot and choose internal and bed ranges
fmcw_plot('a','filelist',datafilename{1},'burstlist',1)
uiwait(msgbox('Click on the internal to zoom','modal'))
[xi,~] = ginput(1);
% Zoom in
xrange = 75;
set(gca,'xlim',[xi-xrange xi+xrange]);
uiwait(msgbox('Click on the internal range to select','modal'))
[xi,yi] = ginput(2);
axis tight

uiwait(msgbox('Click on the bed to zoom','modal'))
[xb,~] = ginput(1);
% Zoom in
set(gca,'xlim',[xb-xrange xb+xrange]);
uiwait(msgbox('Click on the bed range to select','modal'))
[xb,yb] = ginput(2);
axis tight

hold on
ih = plot(xi,yi,'g+');
bh = plot(xb,yb,'r+');
legend([ih bh],'internal range','bed range')

%% Write the output file
try
    disp(['writing filelist to ' filename])
    [fid,msg] = fopen(filename,'wt');
    if fid == -1
        error(msg)
    end
        
    fprintf(fid,'%s\n','fmcw_timeseries_config_v2: config file for batch processing BAS FMCW radar data with fmcw_timeseries');
    fprintf(fid,'%s\n',notes);
    fprintf(fid,'%s\n',['vsr: ' num2str(vsr)]);
    fprintf(fid,'%s','intdep: ');
    fprintf(fid,'%6.2f ',xi(1)); % internal depth
    fprintf(fid,'%6.2f\n',xi(2)); % internal depth
    fprintf(fid,'%s','beddep: ');
    fprintf(fid,'%6.2f ',xb(1)); % bed depth
    fprintf(fid,'%6.2f\n',xb(2)); % bed depth
    for ii = 1:length(datafilename)-1
        fprintf(fid,'%s\n',datafilename{ii});
    end
    fprintf(fid,'%s',datafilename{end});
    fclose(fid);
    %keyboard
catch ME
    fclose(fid);
    rethrow(ME)
end