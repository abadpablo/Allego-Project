function [performance] = analysisTMaze_txt(txt, varargin)
%
% Analysis of the txt file created by Arduino
% USAGE
%   [performance] = analysisTMaze_txt(txt,varargin)
 
% INPUTS
% txt: txt file obtained with the importdata function. Fields:
%       1. TrialSCNumber (Trial Sample-Choice, same number for the pair)
%       2. TrialNumber (Trial Sample or Choice, different for sample and
%       choice)
%       3. forcedArm (1- right / 2 - left). Arm forced during the sample
%       4. chosenArm (1 - right / 2 - left). Arm chosen during the choice
%       5. TrialCorrect (1 - Correct / 2 - Incorrect). Only in choice
%       6. performance (Total Score of the experiment). Only last value
%       important
%       7. time trial (time between start and end of the trial)
%       8. time sample (time between forcedArm and start). Only sample
%       9. time choice (time between chosenArm and start). Only choice
%       10. time decision (time between chosenArm and forcedArm). This time
%       includes the 5 seconds animal is waiting between sample and choice
%
% OUTPUT
% performance with the following fields:
% performance(%)

% 
%   HISTORY:
%       - Pablo Abad 2021

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'savemat',true,@islogical);
addParameter(p,'forceReload',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
savemat = p.Results.savemat;
forceReload = p.Results.forceReload;

if isempty(txt)    
    if ~isempty(dir([basepath filesep '*TMaze.txt']))
        disp('txt file detected! Loading file.')
        file = dir([basepath filesep '*TMaze.txt']);
        txt = importdata(file.name);
    end
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

%% 
trialSCNumber_col = 1;
trialNumber_col = 2;
forcedArm_col = 3;
chosenArm_col = 4;
trialCorrect_col = 5;
performance_col = 6;
timeTrial_col = 7;
timeSample_col = 8;
timeChoice_col = 9;
timeDecision_col = 10;

%% Post-processing analysis of the YMaze

% Let's eliminate the first row (trial number 0)
txt(1,:) = [];
score = txt(size(txt,1),performance_col);
% Odd trials are the sample and even trials are the choice

%% RIGHT TRIALS VS WRONG TRIALS
forcedArm = [];
chosenArm = [];
trialCorrect = [];
timeTrial_sample = [];
timeTrial_choice = [];
timeForced = [];
timeChoice = [];
timeDecision = [];
trialCorrect_index = [];

timeTrial = [];

for i=1:size(txt,1)
    if mod(i,2) ~= 0 % odd trial / sample
        forcedArm = [forcedArm; txt(i,forcedArm_col)];
        timeTrial_sample = [timeTrial_sample; txt(i,timeTrial_col)];
        timeForced = [timeForced; txt(i,timeSample_col)];
    end
    if mod(i,2) == 0 % even trial / choice
        chosenArm = [chosenArm; txt(i,chosenArm_col)];
        trialCorrect = [trialCorrect; txt(i,trialCorrect_col)];
        trialCorrect_index = [trialCorrect_index; i];
        timeTrial_choice = [timeTrial_choice; txt(i,timeTrial_col)];
        timeChoice = [timeChoice; txt(i,timeChoice_col)];
        timeDecision = [timeDecision; txt(i,timeDecision_col)];
    end
end

timeTrial = [timeTrial_sample'; timeTrial_choice'];
    
%%    
right_index = find(trialCorrect == 1);
wrong_index = find(trialCorrect == 0);

%% RIGHT TRIALS
right.sample.forcedArm = forcedArm(right_index);
right.sample.timeTrial = timeTrial_sample(right_index);
right.sample.timeForced = timeForced(right_index);
right.choice.chosenArm = chosenArm(right_index);
right.choice.trialCorrect = trialCorrect(right_index);
right.choice.timeTrial = timeTrial_choice(right_index);
right.choice.timeChoice = timeChoice(right_index);
right.choice.timeDecision = timeDecision(right_index);
right.trials = length(right_index);

%% WRONG TRIALS

wrong.sample.forcedArm = forcedArm(wrong_index);
wrong.sample.timeTrial = timeTrial_sample(wrong_index);
wrong.sample.timeForced = timeForced(wrong_index);
wrong.choice.chosenArm = chosenArm(wrong_index);
wrong.choice.trialCorrect = trialCorrect(wrong_index);
wrong.choice.timeTrial = timeTrial_choice(wrong_index);
wrong.choice.timeChoice = timeChoice(wrong_index);
wrong.choice.timeDecision = timeDecision(wrong_index);
wrong.trials = length(wrong_index);

performance.right = right;
performance.wrong = wrong;
performance.score = score;


end