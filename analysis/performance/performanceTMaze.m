function performance = performanceTMaze(varargin)

% Need to include the zones for the tracking

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'fs',30,@isnumeric);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'digitalIn',[],@isstruct);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'useDigitalIn',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
fs = p.Results.fs;
forceReload = p.Results.forceReload;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
digitalIn = p.Results.digitalIn;
tracking = p.Results.tracking;
useDigitalIn = p.Results.useDigitalIn;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*Performance.Behavior.mat'])) || forceReload
    disp('Performance already detected! Loading file.');
    file = dir([basepath filesep '*Performance.Behavior.mat']);
    load(file.name);
    return
end

if isempty(digitalIn)
    digitalIn = bz_getDigitalIn('all');
end

if isempty(tracking)
    tracking = getSessionTracking();
end


% Let's try to take profit of the txt file generated in Arduino

if ~isempty(dir([basepath filesep '*TMaze.txt']))
    disp('txt file detected! Loading file.')
    file = dir([basepath filesep '*TMaze.txt']);
    txt = importdata(file.name);
end


% Let's analyze the txt file
tmaze = analysisTMaze_txt(txt);


%% Let's try with the digital Inputs

if length(digitalIn.timestampsOn) >= 6 && useDigitalIn
    disp('Performance can be computed with digital Inputs')
    
    
else
    disp('Peformance can NOT be computed with digital Inputs')
    
    inLeftReward = 1*tracking.zone.inLeftReward;
    inRightReward = 2*tracking.zone.inRightReward;
    inDecision = 3*tracking.zone.inDecision;
    inStarting = 4*tracking.zone.inStarting;
    transition = zeros(2,length(inLeftArm));
    
    for i=1:length(inLeftArm)
        if inLeftArm(i) > 0
            transition(1,i) = inLeftArm(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inRightArm(i) > 0
            transition(1,i) = inRightArm(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inDecision(i) > 0
            transition(1,i) = inDecision(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inStarting(i) > 0
            transition(1,i) = inStarting(i);
            transition(2,i) = tracking.timestamps(i);
        end  
    end
    
    transitions = analysisTransitions(transition,'tracking',tracking);   
end
    
    
performance = transitions;
sampleVSchoice = bz_computeSampleChoice('performance',performance);
performance.sampleVSchoice = sampleVSchoice;

filename = strsplit(basepath,filesep);
filename = filename{end};

performance.folder = filename;
performance.paradigm = 'TMaze';
performance.timestamps = tracking.timestamps;

if saveMat
    
    save([basepath filesep filename '.Performance.Behavior.mat'],'performance');
end

end
