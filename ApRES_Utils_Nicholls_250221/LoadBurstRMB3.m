function vdat = LoadBurstRMB3(Filename, Burst, SamplesPerChirp)

% vdat = LoadBurstRMB3(Filename, Burst, SamplesPerChirp)
%
% Read FMCW data file from RMB-A? (Data from Jan 2013)

WperChirpHdr = 0;
MaxHeaderLen = 1024;
burstpointer = 0;
vdat.Code = 0;
fid = fopen(Filename,'r');
if fid >= 0
    fseek(fid,0,'eof');
    filelength = ftell(fid);
    BurstCount = 1;
    while BurstCount <= Burst && burstpointer <= filelength - MaxHeaderLen
        fseek(fid,burstpointer,'bof');
        A = fread(fid,MaxHeaderLen,'*char');
        A = A';
        searchind = strfind(A, 'Samples:');
        if ~isempty(searchind)
            try
                searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                vdat.Nsamples = sscanf(A(searchind(1)+8:searchCR(1)+searchind(1)),'%d');
                WperChirpCycle = vdat.Nsamples + WperChirpHdr;
                searchind = strfind(A, 'SubBursts in burst:');
                searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                vdat.ChirpsInBurst = sscanf(A(searchind(1)+20:searchCR(1)+searchind(1)),'%d');
                
                searchind = strfind(A, '*** End Header ***');
                
                burstpointer = burstpointer + searchind(1) + 20;
            catch
                vdat.Code = -2;
                vdat.Burst = BurstCount;
                return
            end
        end
        WordsPerBurst = vdat.ChirpsInBurst * WperChirpCycle;
        if BurstCount < Burst && burstpointer <= filelength - MaxHeaderLen
            burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle*2;
        end
        BurstCount = BurstCount + 1;
    end
    
    % Extract remaining information from header
    searchind = strfind(A, 'Time stamp:');
    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
    try
        td = sscanf(A(searchind(1)+11:searchCR(1)+searchind(1)),...
            '%d-%d-%d %d:%d:%d');
        vdat.TimeStamp = datenum(td(1),td(2),td(3),td(4),td(5),td(6));
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Temperature 1:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Temperature_1 = sscanf(A(searchind(1)+14:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Temperature 2:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Temperature_2 = sscanf(A(searchind(1)+14:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Battery voltage:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.BatteryVoltage = sscanf(A(searchind(1)+16:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Attenuator 1:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Attenuator_1 = sscanf(A(searchind(1)+13:searchCR(1)+searchind(1)),'%f,',4);
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Attenuator 2:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Attenuator_2 = sscanf(A(searchind(1)+13:searchCR(1)+searchind(1)),'%f,',4);
    catch err
        vdat.Code = 1;
    end
    
    fseek(fid,burstpointer-1,'bof');
    if BurstCount == Burst+1
        [vdat.v count] = fread(fid,WordsPerBurst,'*int16','ieee-le');
        if count < WordsPerBurst
            vdat.Code = 2;
        end
        vdat.v = double(vdat.v);
        vdat.v = vdat.v * 2.5 / 2^16 + 1.25;
        vdat.Startind = ((WperChirpHdr+1):WperChirpCycle:WperChirpCycle*vdat.ChirpsInBurst)';
        vdat.Endind = vdat.Startind + SamplesPerChirp - 1;
        vdat.Burst = Burst;
    else
% Too few bursts in file
        vdat.Burst = BurstCount - 1;
        vdat.Code = -4;
    end
    fclose(fid);
else
% Unknown file
    vdat.Code = -1;
end

% Clean temperature record (wrong data type?)
bti1 = find(vdat.Temperature_1>300); % bad temperature indices
if ~isempty(bti1)
    %disp('Cleaning temperature over 300C')
    vdat.Temperature_1(bti1) = vdat.Temperature_1(bti1)-512;
end
bti2 = find(vdat.Temperature_2>300); % bad temperature indices
vdat.Temperature_2(bti2) = vdat.Temperature_2(bti2)-512;
