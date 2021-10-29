function vdat = LoadBurstRMB5(Filename, Burst, SamplingFrequency)

% vdat = LoadBurstRMB5(Filename, Burst, SamplesPerChirp)
%
% Read FMCW data file from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)

% Corrected so that Sampling Frequency has correct use (ie, not used in
% this case)

MaxHeaderLen = 1500;
burstpointer = 0;
vdat.Code = 0;
fid = fopen(Filename,'r');
if fid >= 0
    fseek(fid,0,'eof');
    filelength = ftell(fid);
    BurstCount = 1;
    while BurstCount <= Burst && burstpointer <= filelength - MaxHeaderLen
        
        % Read header
        fseek(fid,burstpointer,'bof');
        A = fread(fid,MaxHeaderLen,'*char');
        A = A';
        
        % Read all aspects of header that define length of burst
        SearchString = 'N_ADC_SAMPLES=';
        searchind = strfind(A,SearchString);
        if ~isempty(searchind)
            try
                searchCR = strfind(A(searchind(1):end),char(10));
                vdat.Nsamples = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d');
                WperChirpCycle = vdat.Nsamples;
                SearchString = 'NSubBursts=';
                searchind = strfind(A,SearchString);
                searchCR = strfind(A(searchind(1):end),char(10));
                vdat.SubBurstsInBurst = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d');
                SearchString = 'Average=';
                searchind = strfind(A, SearchString);
                if isempty(searchind)
                    vdat.Average = 0; %cls 9/jan/14 -average not included in mooring deploy
                else
                    searchCR = strfind(A(searchind(1):end),char(10));
                    vdat.Average = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d');
                end
                
                SearchString = 'nAttenuators=';
                searchind = strfind(A, SearchString);
                searchCR = strfind(A(searchind(1):end),char(10));
                vdat.NAttenuators = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%d',1);
                
                SearchString = 'TxAnt=';
                searchind = strfind(A, SearchString);
                searchCR = strfind(A(searchind(1):end),char(10));
                vdat.TxAnt = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)-1),'%d,',8);
                
                SearchString = 'RxAnt=';
                searchind = strfind(A, SearchString);
                %                searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                searchCR = strfind(A(searchind(1):end),char(10));
                vdat.RxAnt = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)-1),'%d,',8);
                
                % Calculate number of chirps in each burst
                ind = find(vdat.TxAnt~=1);
                vdat.TxAnt(ind) = [];
                ind = find(vdat.RxAnt~=1);
                vdat.RxAnt(ind) = [];
                
                if vdat.Average
                    vdat.ChirpsInBurst = 1;
                else
                    vdat.ChirpsInBurst = vdat.SubBurstsInBurst * length(vdat.TxAnt) * ...
                        length(vdat.RxAnt) * vdat.NAttenuators;
                end
                
                % Move burstpointer to end of header
                SearchString = '*** End Header ***';
                searchind = strfind(A, SearchString);
                burstpointer = burstpointer + searchind(1)-1 + length(SearchString)+2;
            catch
                vdat.Code = -2;
                vdat.Burst = BurstCount;
                return
            end
        end
        
        % (WperChirpCycle is the number of samples in the chirp)
        WordsPerBurst = vdat.ChirpsInBurst * WperChirpCycle;
        if BurstCount < Burst && burstpointer <= filelength - MaxHeaderLen
            if vdat.Average
                burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle*4;
            else
                burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle*2;
            end
        end
        BurstCount = BurstCount + 1;
    end
    
    % Extract remaining information from header for required burst
    try
        
        SearchString = 'Attenuator1=';
        searchind = strfind(A, SearchString);
        searchCR = strfind(A(searchind(1):end),char(10));
        vdat.Attenuator_1 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f,',vdat.NAttenuators);
        
        SearchString = 'AFGain=';
        searchind = strfind(A, SearchString);
        %                searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        searchCR = strfind(A(searchind(1):end),char(10));
        vdat.Attenuator_2 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f,',vdat.NAttenuators);
        
        SearchString = 'Time stamp=';
        searchind = strfind(A, SearchString);
        if isempty(searchind)
            vdat.Code = -4;
            return
        end
        searchCR = strfind(A(searchind(1):end),char(10));
        td = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),...
            '%d-%d-%d %d:%d:%d');
        vdat.TimeStamp = datenum(td(1),td(2),td(3),td(4),td(5),td(6));
    catch err
        vdat.Code = 1;
    end
    
    SearchString = 'Temp1=';
    searchind = strfind(A, SearchString);
    try
        searchCR = strfind(A(searchind(1):end),char(10));
        vdat.Temperature_1 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    SearchString = 'Temp2=';
    searchind = strfind(A, SearchString);
    try
        %        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        searchCR = strfind(A(searchind(1):end),char(10));
        vdat.Temperature_2 = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    SearchString = 'BatteryVoltage=';
    searchind = strfind(A, SearchString);
    try
        searchCR = strfind(A(searchind(1):end),char(10));
        vdat.BatteryVoltage = sscanf(A(searchind(1)+length(SearchString):searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    
    fseek(fid,burstpointer,'bof');
    if BurstCount == Burst + 1
        if vdat.Average == 2
            [vdat.v count] = fread(fid,WordsPerBurst,'*uint32','ieee-le');
        elseif vdat.Average == 1
            fseek(fid,burstpointer,'bof');
            [vdat.v count] = fread(fid,WordsPerBurst,'*real*4','ieee-le');
        else
            [vdat.v count] = fread(fid,WordsPerBurst,'*uint16','ieee-le');
        end
        if count < WordsPerBurst
            vdat.Code = 2;
        end
        vdat.v(vdat.v<0) = vdat.v(vdat.v<0) + 2^16;
        vdat.v = single(vdat.v);
        vdat.v = vdat.v * 2.5 / 2^16;
        if vdat.Average == 2
            vdat.v = vdat.v / (vdat.SubBurstsInBurst * vdat.NAttenuators);
        end
        vdat.Startind = (1:WperChirpCycle:WperChirpCycle*vdat.ChirpsInBurst)';
        vdat.Endind = vdat.Startind + WperChirpCycle - 1;
        vdat.Burst = Burst;
    else
        % Too few bursts in file
        vdat.Burst = BurstCount - 1;
        vdat.Code = -4;
        %keyboard
    end
    fclose(fid);
else
    % Unknown file
    vdat.Code = -1;
end

% Clean temperature record (wrong data type?)
bti1 = find(vdat.Temperature_1>300);
if ~isempty(bti1)
    vdat.Temperature_1(bti1) = vdat.Temperature_1(bti1)-512;
end
bti2 = find(vdat.Temperature_2>300);
vdat.Temperature_2(bti2) = vdat.Temperature_2(bti2)-512;
