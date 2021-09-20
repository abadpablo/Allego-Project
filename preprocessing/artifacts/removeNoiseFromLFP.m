function removeNoiseFromLFP(basepath,varargin)
% Edit dat files with different options.
%
% USAGE
%   editDatFile(basepath,ints,varargin)
% 
% INPUT
% basepath      If not provided, takes pwd
% threshold     Intervals to be edited.
% method        'substractMedian' or 'substratMean' (defaut)
% keepDat       Default, false.
%
% <optional>
% option        'remove' or 'zeroes' (default). 
%
% Manu Valero-BuzsakiLab 2021
%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'ch','all');
addParameter(p,'method','substractMean',@ischar);
addParameter(p,'keepDat',false,@islogical);

warning('Performing median/mean substraction!! Dat file will be compromised!! ');
parse(p,varargin{:});
ch = p.Results.ch;
method = p.Results.method;
basepath = p.Results.basepath;
keepDat = p.Results.keepDat;

% Get elements
prevPath = pwd;
cd(basepath);

xml = LoadParameters;
fileTargetAmplifier = dir('amplifier*.dat');
if isempty(fileTargetAmplifier)
    filename = split(pwd,filesep); filename = filename{end};
    fileTargetAmplifier = dir([filename '*.dat']);
end

if size(fileTargetAmplifier,1) == 0
    error('Dat file not found!!');
end

if ischar(ch) && strcmpi(ch, 'all')
    ch = 1:length(xml.channels);
end
nChannels = xml.nChannels;
session = loadSession();
% duration = 60;
duration = session.general.duration;
frequency = xml.rates.wideband;

fid = fopen(fileTargetAmplifier.name,'r'); 
filename = fileTargetAmplifier.name;
C = strsplit(fileTargetAmplifier.name,'.dat'); 
filenameOut = [C{1} '_temp.dat'];
fidOutput = fopen(filenameOut,'a');

lfp = bz_GetLFP('all');
data = lfp.data';

if exist([basepath filesep strcat(session.general.name,'.MergePoints.events.mat')],'file')
    
	load(strcat(session.general.name,'.MergePoints.events.mat'));
end



while 1    
    data = fread(fid,[nChannels frequency*duration], 'int16');
    if isempty(data)
        break;
    end
    
    if strcmpi('substractMedian',method)
        m_data = median(data(ch,:));
    elseif strcmpi('substractMean',method)
        m_data = mean(data(ch,:));
    end

    for ii = 1:length(ch)
        data_aux(ch(ii),:) = int16(double(data(ch(ii),:)) - double(m_data));
    end
    fwrite(fidOutput,data,'int16');
    
end
fclose(fid);
fclose(fidOutput);

if ~keepDat
    copyfile(filename, [C{1} '_original.dat']);
end
delete(filename);
movefile(filenameOut, filename);

cd(prevPath);

end