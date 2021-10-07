function [transitions] = analysisTransitions(transition, varargin)
%
% Goes one-by-one step defining triades and
% correct vs uncorrect trials for YMaze paradigm
% USAGE
%   [transitions] = analysisTransitions(transition,varargin)
 
% INPUTS
% transition: variable with ROI and time for each timestamp
% tracking: tracking structure obtained with getSessionTracking()
%
%
%
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

frames_stemArm = find(transition(1,:) == 1);
frames_leftArm = find(transition(1,:) == 3);
frames_rightArm = find(transition(1,:) == 2);
frames_centerMaze = find(transition(1,:) == 4);

frames_noROI = find(transition(1,:) == 5);

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
% 1 -> Stem Arm
% 2 -> Right Arm
% 3 -> Left Arm
% 4 -> Center Maze
% 5 -> Zone of transition between ROIs
tran = transitions_;
tran(tran == 0) = [];

% Now we change 5 for 0 (noZone)
tran(tran == num2str(5)) = num2str(0);
%% Let's compute those ballistic epochs and doubted epochs when the animal goes out of some ROI but then goes to the same ROI

[ballistic, ballistic_perc, doubted, doubted_perc] = bz_computeBallisticVsDoubted(tran);

%%
% Now lets clean those transitions where the mouse gout ut and got in on
% the same ROI (cleaning of the doubted trials and we only take into
% account when the animal finally takes the decission
cleanTransitions = {'1010','2020','3030','4040'};
for i = 1:length(tran)-3
    cj{i} = contains(cleanTransitions,tran(i:i+4-1));
end

for i=1:length(cj)
    
    if (find(cj{i} == 1))
        a(i) = i;
    end

end

if exist('a')
    a(a == 0) = [];
    b = a+1;
    a_ = [a,b];
    a_ = sort(a_);
    transitions_out = tran;
    transitions_out(a_) = [];
transitions_fr = transitions_frame;
transitions_fr(transitions_fr == 0) = [];
epochs_eliminated = transitions_fr(a_);
epochs_doubted = epochs_eliminated;
transitions_fr_out = transitions_fr;
for i=1:length(epochs_eliminated)
    a(i) = find(transitions_fr_out == epochs_eliminated(i));
end
transitions_fr_out(a) = [];

else
    transitions_fr = transitions_frame;
    transitions_fr(transitions_fr == 0) = [];
    transitions_out = tran;
    transitions_fr_out = transitions_fr;
end
% we need to keep track of those epochs eliminated





% End of variable transition erros

% Let's eliminate the value '0' from the transitions variable
transitions_out_nonzero = transitions_out;
transitions_out_nonzero(transitions_out_nonzero == '0') = [];

z = find(transitions_out == '0')
transitions_fr_out_nonzero = transitions_fr_out;
transitions_fr_out_nonzero(z) = [];

%% Compute final output analysis
rightTransitions = {'24143';'24341';'34142';'34241';'14243';'14342'};
index = length(rightTransitions{1});
score = 0;
i = 1;
if transitions_out_nonzero(1) == '4'
    transitions_out_nonzero(1) = [];
    transitions_fr_out_nonzero(1) = [];
end
%%
while i < length(transitions_out_nonzero)-index+1
    
    index_epoch = i:i+index-1
    transitions_epoch{i} = transitions_out_nonzero(i:i+index-1)
    jt = contains(rightTransitions,transitions_out_nonzero(i:i+index-1)) 
    a = find(jt == 1)
    % Added this feature to know which has been the transition
    if ~isempty(a)
        which_transitions{i} = rightTransitions(jt)
    end
%Right transitions
        if ~isempty(a)
%             right_transitions_fr{i} = transitions_fr_out(i:i+index-1)
            % Here we pick until the animal goes to the next ROI, because
            % otherwise we only keep until when the animal first enters the
            % ROI
            right_transitions_fr{i} = transitions_fr_out_nonzero(i:i+index)
            transitions_out_nonzero(i:i+index-1)
            score = score + 1
            output_score(i) = score
%             right_transitions_epoch{i} = transitions_out(i:i+index-1)
            right_transitions_epoch{i} = transitions_out_nonzero(i:i+index-1)
            i = i + 2 % we pick next value (omitting the '0' value)
        else %if its a wrong transition, then we pick the next value,
            transitions_out_nonzero(i:i+index-1)
            wrong_transitions_fr{i} = transitions_fr_out_nonzero(i:i+index)
            wrong_transitions_epoch{i} = transitions_out_nonzero(i:i+index-1)

              i = i+2

        end
    
end

if exist('output_score')
%%
output_score = max(output_score);
else
    output_score = 0;
end

%Let's compute the entrances in each ROI

entrances = transitions_out_nonzero;    
stemROI_output.entrances = length(find(entrances == '1')); % Number of entrances in StemArm
leftROI_output.entrances = length(find(entrances == '3')); % Number of entrances in LeftArm
rightROI_output.entrances = length(find(entrances == '2')); % Number of entrances in RightArm
centerROI_output.entrances = length(find(entrances == '4')); % Number of entrances in CenterMaze

figure,
subplot(2,2,1)
plot_numberEntrances(stemROI_output.entrances,leftROI_output.entrances,rightROI_output.entrances, centerROI_output.entrances);
subplot(2,2,2)
% timeInROIs = plot_timeSpentROIs(frames_stemArm, frames_leftArm, frames_rightArm, frames_centerMaze);
timeInROIs = bz_timeSpentROIs(frames_stemArm,frames_leftArm,frames_rightArm,frames_centerMaze);
% Heat map of the position of the animal during the trial
position = [tracking.position.x, tracking.position.y];
subplot(2,2,3)
hist3(position,'Nbins', [20 20],'CdataMode', 'auto')
colorbar
view(0,-90), xlabel('X position'), ylabel('Y position'), title('Position heat map')

saveas(gcf,'Number Entrances.png')
%% Let's save all the variables we are interested in in a struct

% Performance in spontaneus alternation
% SA(%) = (spontaneous alternations / total number of arm entries - 2)*100;
% We need to take into account that the last transition is repeated

lastEntry = transitions_out_nonzero(end);
lastEntry = str2num(lastEntry);

switch (lastEntry)
    case 1 % stemArm
        stemROI_output.entrances = stemROI_output.entrances - 1;
    case 3 % leftArm
        leftROI_output.entrances = leftROI_output.entrances -1;
    case 2 % right Arm
        rightROI_output.entrances = rightROI_output.entrances -1;
    case 4 % centrMaze
        centerROI_output.entrances = centerROI_output.entrances -1;
end

armEntries = leftROI_output.entrances + rightROI_output.entrances + stemROI_output.entrances;

performance = (output_score/(armEntries-2))* 100;

if ~exist('wrong_transitions_epoch','var')
    wrong_transitions_epoch = NaN;
    wrong_transitions_fr = NaN;
end

if ~exist('right_transitions_epoch','var')
    right_transitions_epoch = NaN;
    right_transitions_fr = NaN;
end

transitions = [];

transitions.performance = performance;
transitions.transitions_out_nonzero = transitions_out_nonzero;
transitions.transitions_fr_out_nonzero = transitions_fr_out_nonzero;
transitions.transitions_times_out_nonzero = tracking.timestamps(transitions_fr_out_nonzero)';
transitions.right_transitions_epoch = right_transitions_epoch;
transitions.right_transitions_fr = right_transitions_fr;
transitions.wrong_transitions_epoch = wrong_transitions_epoch;
transitions.wrong_transitions_fr = wrong_transitions_fr;
transitions.entrances.stem = stemROI_output.entrances;
transitions.entrances.left = leftROI_output.entrances;
transitions.entrances.right = rightROI_output.entrances;
transitions.entrances.center = centerROI_output.entrances;
transitions.transitions_out = transitions_out;
transitions.transitions_fr_out = transitions_fr_out;
transitions.transitions_times_out = tracking.timestamps(transitions_fr_out)';
transitions.frames_stem = frames_stemArm;
transitions.frames_left = frames_leftArm;
transitions.frames_right = frames_rightArm;
transitions.times_stem = tracking.timestamps(frames_stemArm);
transitions.times_left = tracking.timestamps(frames_leftArm);
transitions.times_right = tracking.timestamps(frames_rightArm);
transitions.ballistic = ballistic;
transitions.doubted = doubted;
transitions.ballistic_perc = ballistic_perc;
transitions.doubted_perc = doubted_perc;
transitions.times.stem = timeInROIs{1};
transitions.times.left = timeInROIs{2};
transitions.times.right = timeInROIs{3};
transitions.times.center = timeInROIs{4};

% Need to change right transitions and wrong transitions and add the times
count = 1;
for i=1:length(right_transitions_epoch)
    if ~isempty(right_transitions_epoch{i})
        right_epoch{count} = right_transitions_epoch{i};
        right_fr{count} = right_transitions_fr{i};
        right_times{count} = tracking.timestamps(right_transitions_fr{i});
        count = count+1;
    end
end

    
count = 1;
for i=1:length(wrong_transitions_epoch)
    if ~isempty(wrong_transitions_epoch{i})
        wrong_epoch{count} = wrong_transitions_epoch{i};
        wrong_fr{count} = wrong_transitions_fr{i};
        wrong_times{count} = tracking.timestamps(wrong_transitions_fr{i});
        count = count+1;
    end
end

transitions.right_epoch = right_epoch;
transitions.right_fr = right_fr;
transitions.right_times = right_times;
transitions.wrong_epoch = wrong_epoch;
transitions.wrong_fr = wrong_fr;
transitions.wrong_times = wrong_times;
transitions.rightTriades = length(right_epoch);
transitions.wrongTriades = length(wrong_epoch);

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