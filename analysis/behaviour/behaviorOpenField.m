function [behavior] = behaviorOpenField(varargin)
% Creates behavior as linearizeLinearMaze to run the following scripts, no
% linearization
%
% USAGE
%
%   [behavior] = linearizeLinearMaze(varargin)
%
% INPUTS
% (OPTIONAL)
% basePath            -(default: pwd) basePath for the recording file, in
%                        buzcode format. 
% tracking            - Tracking structure, with a timestamps field and a position field that
%                        contains x (1xC) and y (1xC) subfields. By default, runs LED2Tracking 
%                        to get it.
% digitalIn           - DigitalIn structure with T maze convention:
%                                 1. Basler,            2. maze LEd, 
%                                 3. Left arm,          4.Righ arm
% editLOI             - Edit loaded Line of interest (LOI). 
% saveMat             - Default true
% forceReload         - Default false
% verbose             - Default true
% 
% OUTPUT
%                     - Behavior structure with the following fields updated:
% 
% behavior.timestamps                Total behavioral timestamps
% behavior.position.lin              Linearized position in cm
% behavior.position.x                X coordinates of tracking, in cm/norm
% behavior.position.y                Y coordinates, in cm/norm 
% behavior.masks.arm                 Code for map maze arms (ej, 0 is left, 1 is arm)
% behavior.maps                      Cell array as [time position], one cell/map
% behavior.description               
% behavior.events
% behavior.trials.startPoint         Trial epochs, defined as epochs
%                                       between one side to the other side
% behavior.trials.endDelay           Trial epochs, defnied as delays door openings.
% behavior.trials.arm                (1x#trials). Trial's arm (ej 0 left, 1 right)
% 
%   Manu Valero 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with inputs
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'digitalIn',[],@isstruct);
addParameter(p,'editLOI',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'verbose',false,@islogical);

parse(p,varargin{:});
tracking = p.Results.tracking;
basepath = p.Results.basepath;
digitalIn = p.Results.digitalIn;
editLOI = p.Results.editLOI;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
verbose = p.Results.verbose;

if ~isempty(dir('*OpenField.Behavior.mat')) && ~forceReload 
    disp('OpenField already computed! Loading file.');
    file = dir('*OpenField.Behavior.mat');
    load(file.name);
    return
end

cd(basepath);
if isempty(tracking)
%     tracking = LED2Tracking;
    tracking = anyMazeTracking([],[]);
    % Insertar mi función
end

if isempty(digitalIn)
    digitalIn = bz_getDigitalIn;
%     digitalIn = findThresholdLinearArm;
end

if isempty(tracking) || isempty(digitalIn)
    warning('Missing components. No behaviour performed?');
    return
end

% get components
average_frame  = tracking.avFrame.r;        % get average frames
xMaze = tracking.avFrame.xSize;
yMaze = tracking.avFrame.ySize;
x = tracking.position.x;
y = tracking.position.y;
t = tracking.timestamps;

if isfield(tracking,'headposition')
    x_head = tracking.headposition.x;
    y_head = tracking.headposition.y;
end

if isfield(tracking,'tailposition')
    x_tail = tracking.tailposition.x;
    y_tail = tracking.tailposition.y;
end

if isfield(tracking,'zone')
    zone = tracking.zone;
end

% generate events
maps = [];
armList = 1;
for ii = 1:length(armList)
%     maps{ii}(:,1) = tracking.timestamps(arm==armList(ii));
%     maps{ii}(:,2) = linCont(arm==armList(ii));
    maps{ii}(:,1) = tracking.timestamps;
    maps{ii}(:,2) = tracking.position.x;
    maps{ii}(:,3) = tracking.position.y;
end

trialMask = 2*ones(size(tracking.timestamps));
direction = 2*ones(size(tracking.timestamps));


behavior.timestamps = tracking.timestamps;

behavior.position.lin = [];
behavior.position.x = tracking.position.x;
behavior.position.y = tracking.position.y;

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

    
behavior.masks.arm = [];
behavior.masks.direction = direction;
behavior.masks.trials = trialMask;
behavior.masks.trialsDirection = [];

behavior.maps = maps;

behavior.description = 'Open Field';

behavior.events.startPoint = [];
behavior.events.rReward = [];
behavior.events.lReward = [];
behavior.events.startDelay = [];
behavior.events.endDelay = [];
behavior.events.intersection = [];

behavior.trials.startPoint = [];
behavior.trials.endDelay = [];
behavior.trials.visitedArm = [];
behavior.trials.choice = [];
behavior.trials.expectedArm = [];

if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.OpenField.Behavior.mat'], 'behavior');
end

end

