function [behavior] = getSessionBehaviour_v2(varargin)
% Compute beahviour over all session subfolders (OpenField, YMaze) and
% linearize for linearTrack and TMaze

% USAGE

%   [behavior] = getSessionBehaviour(varargin)

% INPUTS
% basePath                      (default: pwd) basePath for the recording file, 
%                                    in buzcode format:
% forceReload                   Force detection (boolean, default false)
% verbose                       Default false
% saveMat                       Default true
% maze                          OpenField or YMaze. 
%
% OUTPUT
%
%   Pablo Abad 2021
%   Modified by Pablo Abad 23/08 in order to linearize if necessary
%   (linearTrack and TMaze) or getting the behaviour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'verbose',false,@islogical)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'maze','OpenField',@ischar)
parse(p,varargin{:});
forceReload = p.Results.forceReload;
basepath = p.Results.basepath;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
maze = p.Results.maze;

%% Deal with inputs
filename = split(pwd,filesep); filename = filename{end};
if ~isempty(dir([basepath filesep filename '.Behavior.mat'])) || forceReload
    disp('Behavior already detected! Loading file.');
    file =dir([basepath filesep filename '.Behavior.mat']);
    load(file.name);
    return
end

% if isempty(maze)
%     fileTarget = dir('*SessionArmChoice*');
%     if isempty(fileTarget)
%         maze = 'linearMaze';
%     else
%         maze = 'tMaze';
%     end
% end

% I need to automatize this part in order to get the YMaze or the OpenField
% position. (YMaze should have a SessionArmChoice maybe ¿?¿?¿)

%% Find subfolder recordings
cd(basepath);
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
C = strsplit(sessionInfo.session.name,'_');
sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
count = 1;
for ii = 1:size(sess,1)
    if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']))
        cd([basepath filesep sess(ii).name]);
        fprintf('Computing behaviour in %s folder \n',sess(ii).name);
        
        if ~isempty(dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']))
            file = dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']);
            load(file.name);
        end
        
        
        if strcmpi(tracking.apparatus.name,'yMaze') | strcmpi(tracking.apparatus.name,'YMaze Apparatus')
%             behaviorTemp.(sess(ii).name)= behaviorOpenField('verbose',verbose);
            behaviorTemp.(sess(ii).name)= behaviorYMaze('verbose',verbose); %% I need to create this function for the Ymaze paradigm
        elseif strcmpi(tracking.apparatus.name,'Open Field') || strcmpi(tracking.apparatus.name,'OpenField')
            behaviorTemp.(sess(ii).name)= behaviorOpenField('verbose',verbose,'tracking',tracking);
        elseif strcmpi(tracking.apparatus.name,'linearTrack')
            behaviorTemp.(sess(ii).name) = linearizeLinearMaze('verbose',verbose);
        elseif strcmpi(tracking.apparatus.name,'TMaze')
            behaviorTemp.(sess(ii).name) = linearizeArmChoice('verbose',verbose);
        end
        trackFolder(count) = ii; 
        count = count + 1;
    end
end
cd(basepath);

efields = fieldnames(behaviorTemp);
tracking = getSessionTracking;
x = []; y = []; timestamps = []; 
x_head = []; y_head = [];
x_tail = []; y_tail = [];
lin = []; armMask = []; trialMask = []; recMask = [];
startPoint = []; rReward = []; lReward = []; startDelay = []; endDelay = []; intersection = [];
startPointTrials = []; endDelayTrials = []; visitedArm = []; choice = []; expectedArm = [];
recordings = []; recordingsTrial = []; description = [];
direction = []; trialsDirection = [];

for ii = 1:size(efields,1)
    if ~isfield(behaviorTemp.(efields{ii}).masks, 'direction')
        behaviorTemp.(efields{ii}).masks.direction = behaviorTemp.(efields{ii}).masks.arm;
        behaviorTemp.(efields{ii}).masks.trialsDirection = behaviorTemp.(efields{ii}).masks.trials;
    end    
end

% SubSessionsMask
for ii=1:size(efields,1)
    behaviorTemp.(efields{ii}).events.subSessionMask = tracking.events.subSessionsMask == ii;
    behaviorTemp.(efields{ii}).events.subSessionMask(behaviorTemp.(efields{ii}).events.subSessionMask == 0) = [];
end

if size(tracking.events.subSessions,1) == size(efields,1)
    disp('Correctiong timestamps for session recording...');
    for ii = 1:size(efields,1)
        preRec = tracking.events.subSessions(ii,1);
        timestamps = [timestamps; behaviorTemp.(efields{ii}).timestamps + preRec];
%         timestamps{ii} = behaviorTemp.(efields{ii}).timestamps + preRec;
        x = [x; behaviorTemp.(efields{ii}).position.x];
        y = [y; behaviorTemp.(efields{ii}).position.y];
%         x{ii} = behaviorTemp.(efields{ii}).position.x;
%         y{ii} = behaviorTemp.(efields{ii}.position.y);
        if isfield(behaviorTemp.(efields{ii}),'headposition')
            x_head = [x_head; behaviorTemp.(efields{ii}).headposition.x];
            y_head = [y_head; behaviorTemp.(efields{ii}).headposition.y];
%             x_head{ii} = behaviorTemp.(efields{ii}).headposition.x];
%             y_head{ii} = behaviorTemp.(efields{ii}).headposition.y];
        end
        if isfield(behaviorTemp.(efields{ii}),'tailposition')
            x_tail = [x_tail; behaviorTemp.(efields{ii}).tailposition.x];
            y_tail = [y_tail; behaviorTemp.(efields{ii}).tailposition.y];
        end
        if isfield(behaviorTemp.(efields{ii}),'zone')
            zone = behaviorTemp.(efields{ii}).zone;
        end
            
        lin = [lin; behaviorTemp.(efields{ii}).position.lin];
        armMask = [armMask; behaviorTemp.(efields{ii}).masks.arm];
        trialMask = [trialMask; behaviorTemp.(efields{ii}).masks.trials];
        recMask = [recMask; ii * ones(size(behaviorTemp.(efields{ii}).masks.trials))];
        startPoint = [startPoint; behaviorTemp.(efields{ii}).events.startPoint + preRec];
        rReward = [rReward; behaviorTemp.(efields{ii}).events.rReward + preRec];
        lReward = [lReward; behaviorTemp.(efields{ii}).events.lReward + preRec];
        startDelay = [startDelay; behaviorTemp.(efields{ii}).events.startDelay + preRec];
        endDelay = [endDelay; behaviorTemp.(efields{ii}).events.endDelay + preRec];
        intersection = [intersection; behaviorTemp.(efields{ii}).events.intersection + preRec];
        recordings = [recordings; preRec];
        startPointTrials = [startPointTrials; behaviorTemp.(efields{ii}).trials.startPoint + preRec];
        endDelayTrials = [endDelayTrials; behaviorTemp.(efields{ii}).trials.endDelay + preRec];
        visitedArm = [visitedArm; behaviorTemp.(efields{ii}).trials.visitedArm];
        choice = [choice; behaviorTemp.(efields{ii}).trials.choice];
        expectedArm = [expectedArm; behaviorTemp.(efields{ii}).trials.expectedArm];
        direction = [direction; behaviorTemp.(efields{ii}).masks.direction];
        trialsDirection = [trialsDirection; behaviorTemp.(efields{ii}).masks.trialsDirection'];
        recordingsTrial = [recordingsTrial; ii*ones(size(behaviorTemp.(efields{ii}).trials.startPoint,1),1)];
        description{ii} = behaviorTemp.(efields{ii}).description;
    end
else
    warning('Number of behavioral recordings do not match!')
end


% generate maps, one for each arm and for each recording
maps = [];
% directionList = unique(direction);
directionList = 1;
count = 1;
for ii = 1:length(efields)
    for jj = 1:length(directionList)
%         maps{count}(:,1) = timestamps(direction==directionList(jj) & recMask==ii);
%         maps{count}(:,2) = lin(direction==directionList(jj) & recMask==ii);
        maps{count}(:,1) = timestamps;
        maps{count}(:,2) = x;
        maps{count}(:,3) = y;
        count = count + 1;
    end
end

% populate behavior
behavior.timestamps = timestamps;

behavior.position.lin = lin;
behavior.position.x = x;
behavior.position.y = y;

if exist('x_head','var') && ~isempty(x_head)
    behavior.headposition.x = x_head;
    behavior.headposition.y = y_head;
end

if exist('x_tail','var') && ~isempty(x_tail)
    behavior.tailposition.x = x_tail;
    behavior.tailposition.y = y_tail;
end

if exist('zone','var') && ~isempty(zone)
    behavior.zone = zone;
end

behavior.masks.arm = armMask;
behavior.masks.trials = trialMask;
behavior.masks.direction = direction;
behavior.masks.trialsDirection = trialsDirection;
behavior.masks.recording = recMask;

behavior.maps = maps;

behavior.description = description;

behavior.events.startPoint = startPoint;
behavior.events.rReward = rReward;
behavior.events.lReward = lReward;
behavior.events.startDelay = startDelay;
behavior.events.endDelay = endDelay;
behavior.events.intersection = intersection;
behavior.events.recordings = recordings;

behavior.trials.startPoint = startPointTrials;
behavior.trials.endDelay = endDelayTrials;
behavior.trials.visitedArm = visitedArm;
behavior.trials.choice = choice;
behavior.trials.expectedArm = expectedArm;
behavior.trials.recordings = recordingsTrial;

if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.Behavior.mat'], 'behavior');
end




end
