function H = fmcw_ParametersRMB1b(filename)

% For this version of the radar, system clock and sampling frequency are
% fixed.
fsysclk = 1e9;
fs = 40e3;

[path,~,~] = fileparts(filename);

fid = fopen(filename,'rt');
A = fread(fid,1000,'*char');
A = A';
fclose(fid);
loc = strfind(A,'Samples: ');
searchCR = strfind(A(loc(1):end),[char(10)]);
H.Nsamples = sscanf(A(loc(1)+length(['Samples: ']):searchCR(1)+loc(1)),'%d\n');


% Set default values, to be overwritten if config file found with relevant
% Reg values 
H.noDwellHigh =1;
H.noDwellLow = 0;
H.startFreq = 2.0000e+08;
H.stopFreq = 4.0000e+08;
H.rampUpStep = 4.0000e+03;
H.rampDownStep = 4.0000e+03;
H.tstepUp = 2.0000e-05;
H.tstepDown = 2.0000e-05;

% Attempt to open the config.ini file
fid = fopen([path,'\Config.ini'],'rt');

if fid ~= -1
    
    %read the whole file. Use the = sign to populate a second cell with values
    reg = textscan(fid,'%s %s','delimiter', '=');
    fclose(fid);
    
    for k = 1:length(reg{1})
        switch(char(reg{1,1}{k,1}))
            case 'Reg01' %Control Function Register 2 (CFR2)—Address 0x01 Four bytes
                %Bit 19 (Digital ramp enable)= 1 = Enables digital ramp generator functionality.
                %Bit 18 (Digital ramp no-dwell high) 1 = enables no-dwell high functionality.
                %Bit 17 (Digital ramp no-dwell low) 1 = enables no-dwell low functionality.
                %With no-dwell high, a positive transition of the DRCTL pin initiates a positive slope ramp, which
                %continues uninterrupted (regardless of any activity on the DRCTL pin) until the upper limit is reached.
                %Setting both no-dwell bits invokes a continuous ramping mode of operation;
                val = char(reg{1,2}{k,1});
                ind = strfind(val,'"');
                if ~isempty(ind), val(ind) = []; end
                val = dec2bin(hex2dec(val)); val = fliplr(val);
                H.noDwellHigh = str2num(val(18+1));
                H.noDwellLow = str2num(val(17+1));
                
                %        case 'Reg08' %Phase offset word Register (POW)—Address 0x08. 2 Bytes dTheta = 360*POW/2^16.
                %         val = char(reg{1,2}(k));
                %         H.phaseOffsetDeg = hex2dec(val(1:4))*360/2^16;
                
            case 'Reg0B' %Digital Ramp Limit Register—Address 0x0B
                %63:32 Digital ramp upper limit 32-bit digital ramp upper limit value.
                %31:0 Digital ramp lower limit 32-bit digital ramp lower limit value.
                val = char(reg{1,2}{k,1});
                ind = strfind(val,'"');
                if ~isempty(ind), val(ind) = []; end
                H.startFreq = hex2dec(val(9:end))*fsysclk/2^32;
                H.stopFreq = hex2dec(val(1:8))*fsysclk/2^32;
                
            case 'Reg0C'  %Digital Ramp Step Size Register—Address 0x0C
                %63:32 Digital ramp decrement step size 32-bit digital ramp decrement step size value.
                %31:0 Digital ramp increment step size 32-bit digital ramp increment step size value.
                val = char(reg{1,2}{k,1});
                ind = strfind(val,'"');
                if ~isempty(ind), val(ind) = []; end
                H.rampUpStep = hex2dec(val(9:end))*fsysclk/2^32;
                H.rampDownStep = hex2dec(val(1:8))*fsysclk/2^32;
                
            case 'Reg0D'  %Digital Ramp Rate Register—Address 0x0D
                %31:16 Digital ramp negative slope rate 16-bit digital ramp negative slope value that defines the time interval between decrement values.
                %15:0 Digital ramp positive slope rate 16-bit digital ramp positive slope value that defines the time interval between increment values.
                val = char(reg{1,2}{k,1});
                ind = strfind(val,'"');
                if ~isempty(ind), val(ind) = []; end
                H.tstepUp = hex2dec(val(5:end))*4/fsysclk;
                H.tstepDown = hex2dec(val(1:4))*4/fsysclk;
        end
    end
end

% Set derived configuration parameters
H.fs = fs;
H.nstepsDDS = round(abs((H.stopFreq - H.startFreq)/H.rampUpStep));%abs as ramp could be down
H.chirpLength = H.nstepsDDS * H.tstepUp;

% If number of ADC samples collected is less than required to collect
% entire chirp, set chirp length to length of series actually collected
if H.chirpLength*H.fs > H.Nsamples
    H.chirpLength = H.Nsamples / H.fs;
end
H.nchirpSamples = round(H.chirpLength*H.fs);
H.K = 2*pi*(H.rampUpStep/H.tstepUp); % chirp gradient (rad/s/s)

if(H.stopFreq > 400e6)
    H.rampDir = 'down';
else
    H.rampDir = 'up';
end

if(H.noDwellHigh && H.noDwellLow)
    H.rampDir = 'upDown';
    H.nchirpsPerPeriod = 1.6384/(H.chirpLength);
end
