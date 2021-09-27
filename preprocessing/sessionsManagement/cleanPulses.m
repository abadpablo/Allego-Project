
function [data] = cleanPulses(ts, varargin)
% [data] = cleanPulses(ts, varargin)
%
% Find square tsses
%
% INPUTS
% ts            - List of artifacts that will be removed from data (in seconds)
% <OPTIONALS>
% filename
% basepath
% correctDC     - Logical variable to indicate if DC is corrected, default
%                   false
% ch            - List of channels to clean pulses, default all
% winArt        - window for artefact removal, in seconds, default 0.0005s
% winDC         - window for DC removal, in seconds, default 0.005s
%
% OUTPUTS
% data         
%
% Manu Valero 2018

% Default parameters
filename = split(pwd,'\'); filename = filename{end};

% Parse options
p = inputParser;
addParameter(p,'filename',filename,@isstr);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'correctDC',false, @islogical);
addParameter(p,'ch','all');
addParameter(p,'winArt',0.0005,@isnumeric);
addParameter(p,'winDC',0.005,@isnumeric);
addParameter(p,'dataIn',[],@isnumeric);

parse(p, varargin{:});

filename = p.Results.filename;
basepath = p.Results.basepath;
correctDC = p.Results.correctDC;
ch = p.Results.ch;
winArt = p.Results.winArt;
winDC = p.Results.winDC;
dataIn = p.Results.dataIn;

try [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    fs = sessionInfo.rates.wideband;
    nChannels = sessionInfo.nChannels;
catch
    warning('SessionInfo file not found.');
end

if ischar('ch') && strcmpi(ch, 'all')
    ch = 1:nChannels;
else
    ch = ch + 1;
end
% fs = 1250;
ts = int32(ts * fs);
winArt = winArt * fs;
winDC = winDC * fs;

disp('Creating .dat back up...');
copyfile(strcat(filename,'.dat'),'copy_bin.dat'); % create .dat back up
m = memmapfile(fullfile(basepath,strcat(filename,'.dat')),'Format','int16','Writable', true);
data=reshape(m.Data,nChannels,[]);

% We need to load .lfp file (1250) instead of .dat file (30000)
% lfp = bz_GetLFP([1]);
% lfp.data = lfp.data';
% lfp.timestamps = lfp.timestamps';

for hh = ch
    fprintf('Channel #%i of %i...\n',hh,length(ch))
    for ii = 1:size(ts,2)-1 % tsses
        % remove dc
        if correctDC
            data_(hh,ts(1,ii):ts(2,ii)) = data(hh,ts(1,ii):ts(2,ii)) -...
                median(data(hh,ts(1,ii):ts(1,ii)+winDC)) + ...
                median(data(hh,ts(1,ii)-winDC:ts(1,ii)));
            
            data(hh,ts(ii):ts(ii+1)) = data(hh,ts(ii):ts(ii+1)) -...
                median(data(hh,ts(ii):ts(ii)+winDC)) + ...
                median(data(hh,ts(ii)-winDC:ts(ii)));
        end
        % remove artifacts
        for jj = 1:size(ts,1)
            art = int32([ts(jj,ii) - winArt ts(jj,ii) + winArt]);
            data(hh,art(1):art(2))=int16(interp1(double(art),...
                double(data(hh,art)),double(art(1):art(2))));
        end
    end
end

d = zeros(2,max([size(ts,2) size(ts,2)]));
d(1,1:size(ts,2)) = ts;
d(2,1:size(ts,2)-1) = ts(2:end);

if d(1,1) > d(2,1)
    d = flip(d,1);
end

if d(2,end) == 0; d(2,end)= nan; end

    
data1 = lfp.data;
data2 = lfp.data;
data = lfp.data;

for hh = ch
    fprintf('Channel #%i of %i...\n',hh,length(ch))
    for ii = 1:size(d,2)-1 % tsses
        % remove dc
        if correctDC
            data(hh,d(1,ii):d(2,ii)) = lfp.data(hh,d(1,ii):d(2,ii)) -...
                median(lfp.data(hh,d(1,ii):d(1,ii)+winDC)) + ...
                median(lfp.data(hh,d(1,ii)-winDC:d(1,ii)));
            
        end
        % remove artifacts
        for jj = 1:size(d,1)
            art = int32([d(jj,ii) - winArt d(jj,ii) + winArt]);
            data(hh,art(1):art(2))=int16(interp1(double(art),...
                double(lfp.data(hh,art)),double(art(1):art(2))));
        end
    end
end

tamano = get(0,'ScreenSize');
figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
plot(data(1,:));

% disp('Cleaning data on disk...');
% m.Data=data(:);

end