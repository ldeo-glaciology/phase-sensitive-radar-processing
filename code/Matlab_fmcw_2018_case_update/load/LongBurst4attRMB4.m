function vdat = LongBurst4attRMB4(Filename, Burst, SamplesPerChirp)

% vdat = LongBurst4attRMB4(Filename, Burst, SamplesPerChirp)
%
% Load data from RMB-B burst with 4 attenuator settings. Burst average by stacking to limit memory use
% Used for averaging very long bursts.

WperChirpHdr = 0;
MaxHeaderLen = 500;
burstpointer = 0;
vdat.Code = 0;
NChirp = 0;
vdat.v = zeros(SamplesPerChirp,4);
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
                vdat.ChirpsInBurst = sscanf(A(searchind(1)+19:searchCR(1)+searchind(1)),'%d');
                
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
        vdat.Attenuator_1 = sscanf(A(searchind(1)+13:searchCR(1)+searchind(1)),'%f',4);
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Attenuator 2:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Attenuator_2 = sscanf(A(searchind(1)+13:searchCR(1)+searchind(1)),'%f',4);
    catch err
        vdat.Code = 1;
    end

    fseek(fid,burstpointer-1,'bof');
    if BurstCount == Burst+1
        while vdat.Code == 0 % && NChirp<2
            [tmp count] = fread(fid,SamplesPerChirp,'*int16','ieee-le');
            if count < SamplesPerChirp
                vdat.Code = 2;
            else
                vdat.v(:,mod(NChirp,4)+1) = vdat.v(:,mod(NChirp,4)+1) + double(tmp);
                NChirp = NChirp + 1;
                fseek(fid,burstpointer-1+vdat.Nsamples * NChirp * 2,'bof');
            end
        end
        fclose(fid);
        vdat.v = vdat.v*4/NChirp;
        vdat.v = vdat.v * 2.5 / 2^16;
        vdat.Startind = 1;
        vdat.Endind = vdat.Startind + SamplesPerChirp - 1;
        vdat.Burst = Burst;
    else
% Too few bursts in file
        vdat.Burst = BurstCount - 1;
        vdat.Code = -4;
    end
else
% Unknown file
    vdat.Code = -1;
end

figure;hold;col = 'brck';
for i=1:4
    [spec,Range] = PerformCFFT(vdat.v(:,i));
    plot(Range,20*log10(abs(spec))-10*log10(50)-10*log10(0.001),col(i))
end
