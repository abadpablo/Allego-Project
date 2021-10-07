function [transitions] = analysisTransitions_ObjectRecognition(transition, varargin)
%
% Analysis of Performance in Object Recognition
% USAGE
%   [transitions] = analysisTransitions_ObjectRecognition(transition,varargin)
 
% INPUTS
% transition: variable with ROI and time for each timestamp
% tracking: tracking structure obtained with getSessionTracking()
%
%
% OUTPUT
% transitions with the following fields:
% performance(%)
% transitions(nonzero)
% frames of the transitions(nonzero)
% right transitions
% frames of right transitions
% wrong transitions
% frames of wrong transitions
% entrances in stem
% entrances in left
% entrances in right
% entrances in center
% transitions out (with 0's to analyse centerMaze)
% transitions frames out (with 0's to analyse centerMaze)
% frames stem
% frames left
% frames right
% Ballistic 
% Doubted
% Ballistic (%)
% Doubted (%)
% time in Center
% time in Stem
% time in Left
% time in Right
% 
%   HISTORY:
%       - Pablo Abad 2021

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'savemat',true,@islogical);
addParameter(p,'forceReload',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
tracking = p.Results.tracking;
savemat = p.Results.savemat;
forceReload = p.Results.forceReload;

if isempty(tracking)
    tracking = getSessionTracking();
end

fr = tracking.samplingRate;
filename = strsplit(basepath,filesep);
filename = filename{end};

%% In case performance already exists
if ~isempty(dir([basepath filesep filename, '.*Performance.Behavior.mat'])) || forceReload
    disp('Performance already detected ! Loading file.');
    file = dir([basepath filesep filename, '*Performance.Behavior.mat' ]);
    load(file.name);
    return
end


%% Post-processing analysis of the YMaze

frames_inObject1 = find(transition(1,:) == 1);
frames_inObject2 = find(transition(1,:) == 2);
frames_inOut = find(transition(1,:) == 3);

%% Convert outputROI (num) to transitions (string)

for i = 1:length(transition)-1
    if i== 1
        transitions_(1,1) = num2str(transition(1,1));
%         transitions(2,1) = num2str(i);
%         transitions(3,:) = num2str(tstep_position(str2num(transitions(2,1))));
        
        transitions_frame(1) = find(transition == transition(1,1),1);
        transitions_timestamp(1) = tracking.timestamps(transitions_frame(1));
    end
    
    if ~isequal(transition(1,i+1),transition(1,i))
        transitions_(1,i+1) = num2str(transition(1,i+1));        
        transitions_frame(i+1) = i+1;
        transitions_timestamp(i+1) = tracking.timestamps(transitions_frame(i+1));
    end
end

transitions_(1,length(transition)) = num2str(transition(1,end));
transitions_frame(1,length(transition)) = length(transition);
transitions_timestamp(1,length(transition)) = tracking.timestamps(length(transition));

outROI = transition(1,:);
%%
%Clean of variable transitions to oberve all the transitions
% 1 -> inObject1
% 2 -> inObject2
% 3 -> inOut
tran = transitions_;
tran(tran == 0) = [];
tran_frame = transitions_frame;
tran_frame(tran_frame == 0) = [];
tran_timestamp = transitions_timestamp;
tran_timestamp(tran_timestamp == 0) = [];
tran_timestamp = [tracking.timestamps(1) tran_timestamp];

%% Compute final output analysis
timeObject1 = length(frames_inObject1) / fr;
timeObject2 = length(frames_inObject2) / fr;
timeOut = length(frames_inOut) / fr;

scoreObject1 = timeObject1 / (timeObject1+timeObject2);
scoreObject2 = timeObject2 / (timeObject1+timeObject2);


entrances_Object1 = length(find(tran == '1'));
entrances_Object2 = length(find(tran == '2'));


for i=1:length(tran)
    tran_(i) = str2double(tran(i));
end
tr = [tran_; tran_frame; tran_timestamp];

% Let's take times when animal enters eht ezone ob both objects to be able
% to perform psth


ts_Object1 = tr(3,find(tr(1,:) == 1));
ts_Object2 = tr(3,find(tr(1,:) == 2));

%%
transitions = [];

transitions.entrances.Object1 = entrances_Object1;
transitions.entrances.Object2 = entrances_Object2;

transitions.score.Object1 = scoreObject1;
transitions.score.Object2 = scoreObject2;

transitions.transition = tr;

transitions.time.Object1 = timeObject1;
transitions.time.Object2 = timeObject2;
transitions.time.out = timeOut;

transitions.ts.Object1 = ts_Object1;
transitions.ts.Object2 = ts_Object2;


% if savemat
%     
%     folder = pwd;
%     folder = strsplit(folder,filesep);
%     folder = folder{end};
%     
%     save([folder,'.Performance.Behavior.mat'],'transitions')
%     
% end
end

