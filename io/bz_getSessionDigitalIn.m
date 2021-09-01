function [digitalIn] = bz_getSessionDigitalIn(ch,varargin)
% [pul, val, dur] = getPulses(d,varargin)
%
% Find digital In pulses
%
% INPUTS
% ch            Default all.
% <OPTIONALS>
% fs            Sampling frequency (in Hz), default 30000, or try to
%               recover for rhd 
% offset        Offset subtracted (in seconds), default 0.
% periodLag     How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)    
% filename      File to get pulses from. Default, digitalin.dat file with folder
%               name in current directory
%
%
% OUTPUTS
%               digitalIn - events struct with the following fields
% ints          C x 2  matrix with pulse times in seconds. First column of C 
%               are the beggining of the pulses, second column of C are the end of 
%               the pulses.
% dur           Duration of the pulses. Note that default fs is 30000.
% timestampsOn  Beggining of all ON pulses
% timestampsOff Beggining of all OFF pulses
% intsPeriods   Stimulation periods, as defined by perioLag
% 
% MV-BuzsakiLab 2019
% Based on Process_IntanDigitalChannels by P Petersen

% Parse options
if exist('ch') ~= 1
    ch = 'all';
end

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'fs',30000,@isnumeric)
addParameter(p,'offset',0,@isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'periodLag',1,@isnumeric)


parse(p, varargin{:});
fs = p.Results.fs;
offset = p.Results.offset;
filename = p.Results.filename;
lag = p.Results.periodLag;
basepath = p.Results.basepath;


if ~isempty(dir('*.xml'))
    sess = bz_getSessionInfo(pwd,'noPrompts',true);
end

if ~isempty(dir('*DigitalIn.events.mat'))
    disp('Pulses already detected! Loading file.');
    file = dir('*DigitalIn.events.mat');
    load(file.name);
    return
end

% Find subfolder recordings
cd(basepath);
[sessionInfo] = bz_getSessionInfo(basepath,'noPrompts',true);

if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    count = 1;
    
    for ii=1:size(MergePoints.foldernames,2)
        if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*DigitalIn.events.mat']))
            cd([basepath filesep MergePoints.foldernames{ii}]);
            fprintf('Getting digitalIn in %s folder \n',MergePoints.foldernames{ii});
            tempDigitalIn{count} = bz_getDigitalIn('all','fs',fs);
            digitalInFolder(count) = ii;
            count = count + 1;
        end
    end
    cd(basepath)
else
    error('missing MergePoints, quiting...')
       
end


%% Concatenate and sync timestamps and digitalIn fields

ts = []; subSessions = []; maskSessions = []; maskSessions = [];
timestampsOn = cell(size(tempDigitalIn{1}.timestampsOn)); 
timestampsOff = cell(size(tempDigitalIn{1}.timestampsOff)); 
ints = cell(size(tempDigitalIn{1}.ints)); 
dur = cell(size(tempDigitalIn{1}.dur)); 
intsPeriods = cell(size(tempDigitalIn{1}.intsPeriods));

if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        for ii = 1:length(digitalInFolder)
            if strcmpi(MergePoints.foldernames{digitalInFolder(ii)},tempDigitalIn{ii}.folder)
                
                timestamp_0 = MergePoints.timestamps(digitalInFolder(ii),1);
                
                timestampsOn_val = cellfun(@(x) x+timestamp_0, tempDigitalIn{ii}.timestampsOn,'UniformOutput',false);
                timestampsOn = cellfun(@(x,y) [x y], timestampsOn,timestampsOn_val,'UniformOutput',false);
                
                timestampsOff_val = cellfun(@(x) x+timestamp_0, tempDigitalIn{ii}.timestampsOff,'UniformOutput',false);
                timestampsOff = cellfun(@(x,y) [x y], timestampsOff,timestampsOff_val,'UniformOutput',false);
                
                timestampsInts_val = cellfun(@(x) x+timestamp_0, tempDigitalIn{ii}.ints,'UniformOutput',false);
                ints = cellfun(@(x,y) [x y], ints,timestampsInts_val,'UniformOutput',false);
                
                timestampsDur_val = cellfun(@(x) x+timestamp_0, tempDigitalIn{ii}.dur,'UniformOutput',false);
                dur = cellfun(@(x,y) [x y], dur,timestampsDur_val,'UniformOutput',false);
                
                timestampsIntsPeriods_val = cellfun(@(x) x+timestamp_0, tempDigitalIn{ii}.intsPeriods,'UniformOutput',false);
                intsPeriods = cellfun(@(x,y) [x; y], intsPeriods,timestampsIntsPeriods_val,'UniformOutput',false);
             
                folder{ii} = tempDigitalIn{ii}.folder;
            else
                error('Folders name does not match!!');
            end
        end
    else
        warning('No MergePoints file found. Concatenating timestamps...');
        for ii = 1:length(trackFolder)
            sumTs = max(ts)+ tempTracking{ii}.timestamps;
            subSessions = [subSessions; [sumTs(1) sumTs(end)]];
            ts = [ts; sumTs];
        end
end


digitalIn.timestampsOn = timestampsOn;
digitalIn.timestampsOff = timestampsOff;
digitalIn.ints = ints;
digitalIn.dur = dur;
digitalIn.intsPeriods = intsPeriods;
digitalIn.folders = folder;
%% Saving digitalIn  ...

if exist('digitalIn')==1
    try save([sess.FileName '.DigitalIn.events.mat'],'digitalIn');
    catch
        save('digitalIn.events.mat','digitalIn');
    end
else
    digitalIn = [];
end

end