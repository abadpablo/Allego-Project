function [sampleVSchoice] = bz_computeSampleChoice(varargin)
%
% Computes sample vs choice and assigns every triade (sample and choice) to
% one of the following groups:
%   group1: sample right - choice right
%   group2: sample wrong - choice right
%   group3: sample right - choice wrong
%   group4: sample wrong - choice wrong
%
%   USAGE
%       [sampleVSchoice] = bz_computeSampleChoice(tran,tran_fr)
%
%
%   INPUTS
%       performance: struct containing all info about the performance
%
%   OUTPUT
%       sampleVSchoice with fields:
%
%
%
%
%
%
%
%
%   HISTORY:
%       - Pablo Abad 2021
%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'performance',[],@isstruct);
addParameter(p,'showfig',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
performance = p.Results.performance;
showfig = p.Results.showfig;

filename = strsplit(basepath,filesep);
filename = filename{end};

if isempty(performance)
    disp('Loading performance from recording folder')
    try
        file = dir([basepath filesep filename, '*Performance.Behavior.mat']);
        load(file.name)
    catch
        error('Error: There is no Performance file and not entered to sampleVschoice function')
    end
end

right_transitions = {'14243', '14342', '24341', '24143', '34241', '34142'}
index = 5;
counter1 = 0;
counter2 = 0;
counter3 = 0;
counter4 = 0;
previous_choice = 'right';
tran = performance.transitions_out_nonzero;
tran_fr = performance.transitions_fr_out_nonzero;
tran_times = performance.transitions_times_out_nonzero;

for i=1:2:length(tran)-5
    
    if i==1
        
        epoch = tran(1:5)
        epoch_fr = tran_fr(1:5)
        epoch_times = tran_times(1:5)
    else
        epoch = tran(i:i+index-1)
        epoch_fr = tran_fr(i:i+index-1)
        epoch_times = tran_times(i:i+index-1)
    end
    
    jt = contains(epoch,right_transitions);
    a = find(jt == 1);
    
    
    if ~isempty(a)
        

        if previous_choice == 'right'
            
            counter1 = counter1 + 1;
            
            sRight_cRight.sample{counter1}(1,:) = [str2double(epoch(1)), str2double(epoch(2)), str2double(epoch(3))];
            sRight_cRight.sample{counter1}(2,:) = [epoch_fr(1), epoch_fr(2), epoch_fr(3)];
            sRight_cRight.sample{counter1}(3,:) = [epoch_times(1), epoch_times(2), epoch_times(3)];
            sRight_cRight.choice{counter1}(1,:) = [str2double(epoch(3)), str2double(epoch(4)), str2double(epoch(5))];
            sRight_cRight.choice{counter1}(2,:) = [epoch_fr(3), epoch_fr(4), epoch_fr(5)];
            sRight_cRight.choice{counter1}(3,:) = [epoch_times(3), epoch_times(4), epoch_times(5)];

        elseif previous_choice == 'wrong'
            
            counter2 = counter2 + 1;
            
            sWrong_cRight.sample{counter2}(1,:) = [str2double(epoch(1)), str2double(epoch(2)), str2double(epoch(3))];
            sWrong_cRight.sample{counter2}(2,:) = [epoch_fr(1), epoch_fr(2), epoch_fr(3)];
            sWrong_cRight.sample{counter2}(3,:) = [epoch_times(1), epoch_times(2), epoch_times(3)];
            sWrong_cRight.choice{counter2}(1,:) = [str2double(epoch(3)), str2double(epoch(4)), str2double(epoch(5))];
            sWrong_cRight.choice{counter2}(2,:) = [epoch_fr(3), epoch_fr(4), epoch_fr(5)];
            sWrong_cRight.choice{counter2}(3,:) = [epoch_times(3), epoch_times(4), epoch_times(5)]; 
        end
        
            previous_choice = 'right';

    else
        
        
        if previous_choice == 'right'
            
            counter3 = counter3 + 1;
            
            sRight_cWrong.sample{counter3}(1,:) = [str2double(epoch(1)), str2double(epoch(2)), str2double(epoch(3))];
            sRight_cWrong.sample{counter3}(2,:) = [epoch_fr(1), epoch_fr(2), epoch_fr(3)];
            sRight_cWrong.sample{counter3}(3,:) = [epoch_times(1), epoch_times(2), epoch_times(3)];
            sRight_cWrong.choice{counter3}(1,:) = [str2double(epoch(3)), str2double(epoch(4)), str2double(epoch(5))];
            sRight_cWrong.choice{counter3}(2,:) = [epoch_fr(3), epoch_fr(4), epoch_fr(5)];
            sRight_cWrong.choice{counter3}(3,:) = [epoch_times(3), epoch_times(4), epoch_times(5)];
        elseif previous_choice == 'wrong'
            
            counter4 = counter4 + 1;
            
            sWrong_cWrong.sample{counter4}(1,:) = [str2double(epoch(1)), str2double(epoch(2)), str2double(epoch(3))];
            sWrong_cWrong.sample{counter4}(2,:) = [epoch_fr(1), epoch_fr(2), epoch_fr(3)];
            sWrong_cWrong.sample{counter4}(3,:) = [epoch_times(1), epoch_times(2), epoch_times(3)];
            sWrong_cWrong.choice{counter4}(1,:) = [str2double(epoch(3)), str2double(epoch(4)), str2double(epoch(5))];
            sWrong_cWrong.choice{counter4}(2,:) = [epoch_fr(3), epoch_fr(4), epoch_fr(5)];
            sWrong_cWrong.choice{counter4}(3,:) = [epoch_times(3), epoch_times(4), epoch_times(5)];
        end

            previous_choice = 'wrong';
        
    end
    
end

if counter1 == 0
    sRight_cRight.sample = [];
    sRight_cRight.choice = [];
end

if counter2 == 0
    sWrong_cRight.sample = [];
    sWrong_cRight.choice = [];
end

if counter3 == 0
    sRight_cWrong = [];
    sRight_cWrong = [];
end

if counter4 == 0
    sWrong_cWrong.sample = [];
    sWrong_cWrong.choice = [];
end


sampleVSchoice = [];
sampleVSchoice.group1 = sRight_cRight;
sampleVSchoice.group2 = sWrong_cRight;
sampleVSchoice.group3 = sRight_cWrong;
sampleVSchoice.group4 = sWrong_cWrong;
sampleVSchoice.group1.trials = counter1;
sampleVSchoice.group2.trials = counter2;
sampleVSchoice.group3.trials = counter3;
sampleVSchoice.group4.trials = counter4;

if showfig
    counter = [counter1 counter2 counter3 counter4];

    c = categorical({'Right-Right', 'Wrong-Right', 'Right-Wrong' , 'Wrong-Wrong'});
    figure,
    bar(c,counter)
    title('Number of trials for each group (Sample-Choice)')
    saveas(gcf,'SampleVSchoice.png')
end
end
