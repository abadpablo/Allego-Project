
function [tracking] = getSessionTracking(varargin)
%
% Gets position trackign for each sub-session and concatenate all of them so they are 
% aligned with LFP and spikes. Default is recording with Basler, and requiere avi videos 
% and at least one tracking LED. There is an alternative in case OptiTrack was used. 
% Needs to be run in main session folder. 
%
% USAGE
%
%   [tracking] = getSessionTracking(varargin)
%
% INPUTS
%   basePath       -(default: pwd) basePath for the recording file, in buzcode format:
%   roiTracking    - 2 x R, where 1C is x and 2C is y. By default it
%                   considers the whole video. With the option 'manual' allows to draw
%                   a ROI.
%   roiLED         - 2 x R. 'manual' for drawing the ROI.
%   roisPath       - provide a path with ROI mat files ('roiTRacking.mat'
%                   and 'roiLED.mat'). By default try to find it in
%                   basePath or upper folder.
%   convFact       - Spatial conversion factor (cm/px). If not provide,
%                   normalize maze size.
%   saveMat        - default true
%   forceReload    - default false
%   optitrack      - if optitrack instead of basler was used. Default false
%
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   position.x               - x position in cm/ normalize
%   position.y               - y position in cm/ normalize
%   timestamps      - in seconds, if Basler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%       only for OptiTrack
%   position.z 
%   orientation.x
%   orientation.y
%   orientation.z

%   HISTORY:
%     - Manuel Valero 2019
%     - Added OptiTrack support: 5/20, AntonioFR (STILL NEEDS TESTING)


% HISTORY: 
%     - Pablo Abad 2021
%     - Added AnyMaze support

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'convFact',[],@isnumeric); % 0.1149
addParameter(p,'roiTracking',[],@ismatrix);
addParameter(p,'roiLED',[],@ismatrix);
addParameter(p,'roisPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'optitrack',false,@islogical)
addParameter(p,'anymaze',true,@islogical)

parse(p,varargin{:});
basepath = p.Results.basepath;
convFact = p.Results.convFact;
roiTracking = p.Results.roiTracking;
roiLED = p.Results.roiLED;
roisPath = p.Results.roisPath;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
optitrack = p.Results.optitrack;
anymaze = p.Results.anymaze;

%% In case tracking already exists 
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) || forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

%% Basler tracking 
if ~(optitrack) && ~(anymaze)
    cd(basepath); cd ..; upBasepath = pwd; cd(basepath);
    if isempty(roisPath)
        if exist([basepath filesep 'roiTracking.mat'],'file') || ...
            exist([basepath filesep 'roiLED.mat'],'file')
                roisPath = basepath;
                try load([roisPath filesep 'roiLED.mat'],'roiLED'); end
                load([roisPath filesep 'roiTracking.mat'],'roiTracking');
        elseif exist([upBasepath filesep 'roiTracking.mat'],'file') || ...
            exist([upBasepath filesep 'roiLED.mat'],'file')
                roisPath = upBasepath;
                try load([roisPath filesep 'roiLED.mat'],'roiLED'); end
                load([roisPath filesep 'roiTracking.mat'],'roiTracking');
        end   
    end

    %% Find subfolder recordings
    cd(basepath);
    [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    %C = strsplit(sessionInfo.session.name,'_');
    %sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
    if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        count = 1;
        for ii = 1:size(MergePoints.foldernames,2)
            %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
             if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*Basler*avi']))   
                cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
                fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});
                 tempTracking{count}= LED2Tracking([],'convFact',convFact,'roiTracking',...
                     roiTracking,'roiLED',roiLED,'forceReload',forceReload); % computing trajectory
                 trackFolder(count) = ii; 
                count = count + 1;
            end
        end
        cd(basepath);
    else
        error('missing MergePoints, quiting...');
    end
    
    %% Concatenate and sync timestamps
    ts = []; subSessions = []; maskSessions = [];
    if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        for ii = 1:length(trackFolder)
            if strcmpi(MergePoints.foldernames{trackFolder(ii)},tempTracking{ii}.folder)
                sumTs = tempTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
                subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
                maskSessions = [maskSessions; ones(size(sumTs))*ii];
                ts = [ts; sumTs];
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

    % Concatenating tracking fields...
    x = []; y = []; folder = []; samplingRate = []; description = [];
    for ii = 1:size(tempTracking,2) 
        x = [x; tempTracking{ii}.position.x]; 
        y = [y; tempTracking{ii}.position.y]; 
        folder{ii} = tempTracking{ii}.folder; 
        samplingRate = [samplingRate; tempTracking{ii}.samplingRate];  
        description{ii} = tempTracking{ii}.description;
    end

    tracking.position.x = x;
    tracking.position.y = y;
    tracking.folders = folder;
    tracking.samplingRate = samplingRate;
    tracking.timestamps = ts;
    tracking.events.subSessions =  subSessions;
    tracking.events.subSessionsMask = maskSessions;
end

%% OptiTrack 
if optitrack
    warning('this option has not yet been properly tested');
   % it requieres that csv files have been generated with OptiTrack software and 
   % saved inside sub-session folder 
   
    %% Get tracking in subfolder recordings    
    cd(basepath);
    [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    %C = strsplit(sessionInfo.session.name,'_');
    %sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
    if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        count = 1;
        for ii = 1:size(MergePoints.foldernames,2)
            %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
             if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*Basler*avi']))   
                cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
                fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});
                tempTracking{count} = bz_processConvertOptitrack2Behav(basename,'syncSampFq',sessionInfo.samplingRate);
                trackFolder(count) = ii; 
                count = count + 1;
            end
        end
        cd(basepath);
    else
        error('missing MergePoints, quiting...');
    end    
    
        
    %% Concatenate and sync timestamps
    ts = []; subSessions = []; maskSessions = [];
    if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        for ii = 1:length(trackFolder)
            if strcmpi(MergePoints.foldernames{trackFolder(ii)},tempTracking{ii}.folder)
                sumTs = tempTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
                subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
                maskSessions = [maskSessions; ones(size(sumTs))*ii];
                ts = [ts; sumTs];
            else
                error('Folders names do not match!!');
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

    %% Concatenating tracking fields...
    x = []; y = []; folder = []; samplingRate = []; description = [];
    for ii = 1:size(tempTracking,2) 
        x = [x; tempTracking{ii}.position.x]; 
        y = [y; tempTracking{ii}.position.y]; 
        z = [z; tempTracking{ii}.position.z]; 
        rx = [x; tempTracking{ii}.orientation.rx]; 
        ry = [y; tempTracking{ii}.orientation.ry]; 
        rz = [z; tempTracking{ii}.orientation.rz]; 
        folder{ii} = tempTracking{ii}.folder; 
        samplingRate = [samplingRate; tempTracking{ii}.samplingRate];  
        description{ii} = tempTracking{ii}.description;
    end

    tracking.position.x = x;
    tracking.position.y = y;
    tracking.position.z = z;
    tracking.orientation.rx = rx;
    tracking.orientation.y = ry;
    tracking.orientation.z = rz;    
    tracking.folders = folder;
    tracking.samplingRate = samplingRate;
    tracking.timestamps = ts;
    tracking.events.subSessions =  subSessions;
    tracking.events.subSessionsMask = maskSessions;

end

if anymaze
    %% Get tracking in subfolder recordings
    
    %% Find subfolder recordings
    cd(basepath);
    [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    
     if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        count = 1;
        for ii = 1:size(MergePoints.foldernames,2)
            %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
             if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*.csv*']))   
                cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
                fprintf('Computing tracking in %s folder \n',MergePoints.foldernames{ii});                 
                 tempTracking{count} = anyMazeTracking([],[]);
                 trackFolder(count) = ii;
                count = count + 1;
            end
        end
        cd(basepath);
    else
        error('missing MergePoints, quiting...');
     end
    
     
     %% Concatenate and sync timestamps
    ts = []; subSessions = []; maskSessions = [];
    if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        for ii = 1:length(trackFolder)
            if strcmpi(MergePoints.foldernames{trackFolder(ii)},tempTracking{ii}.folder)
                sumTs = tempTracking{ii}.timestamps + MergePoints.timestamps(trackFolder(ii),1);
                subSessions = [subSessions; MergePoints.timestamps(trackFolder(ii),1:2)];
                maskSessions = [maskSessions; ones(size(sumTs))*ii];
                ts = [ts; sumTs];
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

    % Concatenating tracking fields...
    x = []; y = []; folder = []; speed = []; ace = []; mspeed = []; mace = []; vx = []; vy = []; ax = []; ay = [];
    samplingRate = []; description = []; apparatus = [];
    
    
    x_head = []; y_head = []; speed_head = []; ace_head = []; mspeed_head = []; mace_head = []; vx_head = []; vy_head = []; ax_head = []; ay_head = [];
    x_tail = []; y_tail = []; speed_tail = []; ace_tail = []; mspeed_tail = []; mace_tail = []; vx_tail = []; vy_tail = []; ax_tail = []; ay_tail = [];
    zone = [];
    pixelsmetre = [];
    distance = [];
    
    for ii = 1:size(tempTracking,2) 
%         x{ii} = [x; tempTracking{ii}.position.x]; % Hecho por Pablo Abad
%         y{ii} = [y; tempTracking{ii}.position.y];
        x = [x; tempTracking{ii}.position.x]; 
        y = [y; tempTracking{ii}.position.y];
        vx = [vx; tempTracking{ii}.position.vx];
        vy = [vy; tempTracking{ii}.position.vy];
        ax = [ax; tempTracking{ii}.position.ax];
        ay = [ay; tempTracking{ii}.position.ay];
        ace = [ace; tempTracking{ii}.position.ace];
        speed = [speed; tempTracking{ii}.position.speed];
        mace = [mace; tempTracking{ii}.position.mace];
        mspeed = [mspeed; tempTracking{ii}.position.mspeed];
        
        if isfield(tempTracking{ii},'headposition')
            x_head = [x_head; tempTracking{ii}.headposition.x];
            y_head = [y_head; tempTracking{ii}.headposition.y];
            vx_head = [vx_head; tempTracking{ii}.headposition.vx];
            vy_head = [vy_head; tempTracking{ii}.headposition.vy];
            ax_head = [ax_head; tempTracking{ii}.headposition.ax];
            ay_head = [ay_head; tempTracking{ii}.headposition.ay];
            ace_head = [ace_head; tempTracking{ii}.headposition.ace];
            speed_head = [speed_head; tempTracking{ii}.headposition.speed];
            mace_head = [mace_head; tempTracking{ii}.headposition.mace];
            mspeed_head = [mspeed_head; tempTracking{ii}.headposition.mspeed];
        
        end
        if isfield(tempTracking{ii},'tailposition')
            x_tail = [x_tail; tempTracking{ii}.tailposition.x];
            y_tail = [y_tail; tempTracking{ii}.tailposition.y];
            vx_tail = [vx_tail; tempTracking{ii}.tailposition.vx];
            vy_tail = [vy_tail; tempTracking{ii}.tailposition.vy];
            ax_tail = [ax_tail; tempTracking{ii}.tailposition.ax];
            ay_tail = [ay_tail; tempTracking{ii}.tailposition.ay];
            ace_tail = [ace_tail; tempTracking{ii}.tailposition.ace];
            speed_tail = [speed_tail; tempTracking{ii}.tailposition.speed];
            mace_tail = [mace_tail; tempTracking{ii}.tailposition.mace];
            mspeed_tail = [mspeed_tail; tempTracking{ii}.tailposition.mspeed];
        end
        
        if isfield(tempTracking{ii},'zone')
            zone{ii} = tempTracking{ii}.zone;
        else
            zone{ii} = [];
        end
        
        folder{ii} = tempTracking{ii}.folder; 
        samplingRate = [samplingRate; tempTracking{ii}.samplingRate];  
        description{ii} = tempTracking{ii}.description;
        apparatus{ii} = tempTracking{ii}.apparatus;
        pixelsmetre{ii} = tempTracking{ii}.pixelsmetre;        
        distance = [distance; tempTracking{ii}.distance];

    end

    tracking.position.x = x;
    tracking.position.y = y;
    tracking.position.vx = vx;
    tracking.position.vy = vy;
    tracking.position.ax = ax;
    tracking.position.ay = ay;
    tracking.position.ace = ace;
    tracking.position.speed = speed;
    tracking.position.mspeed = mspeed;
    tracking.position.mace = mace;
    if exist('x_head','var') && ~isempty(x_head)
        tracking.headposition.x = x_head;
        tracking.headposition.y = y_head;
        tracking.headposition.vx = vx_head;
        tracking.headposition.vy = vy_head;
        tracking.headposition.ax = ax_head;
        tracking.headposition.ay = ay_head;
        tracking.headposition.ace = ace_head;
        tracking.headposition.speed = speed_head;
        tracking.headposition.mspeed = mspeed_head;
        tracking.headposition.mace = mace_head;
    end
    if exist('x_tail','var') && ~isempty(x_tail)
        tracking.tailposition.x = x_tail;
        tracking.tailposition.y = y_tail;
        tracking.tailposition.vx = vx_tail;
        tracking.tailposition.vy = vy_tail;
        tracking.tailposition.ax = ax_tail;
        tracking.tailposition.ay = ay_tail;
        tracking.tailposition.ace = ace_tail;
        tracking.tailposition.speed = speed_tail;
        tracking.tailposition.mspeed = mspeed_tail;
        tracking.tailposition.mace = mace_tail;
    
    end
    if exist('zone','var') && ~isempty(zone)
        tracking.zone = zone;
    end
    
    tracking.folders = folder;
    tracking.samplingRate = samplingRate;
    tracking.timestamps = ts;
    tracking.events.subSessions =  subSessions;
    tracking.events.subSessionsMask = maskSessions;
    tracking.apparatus = apparatus;
    tracking.pixelsmetre = pixelsmetre;
    tracking.distance = distance;
    
    
end



%% save tracking 
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
if saveMat
    save([basepath filesep sessionInfo.FileName '.Tracking.Behavior.mat'],'tracking');
end

end

