function performance = performanceYMaze(varargin)

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

%% Let's try with the digital Inputs

if length(digitalIn.timestampsOn) >= 6 && useDigitalIn
    disp('Performance can be computed with digital Inputs')
    
    
else
    disp('Peformance can NOT be computed with digital Inputs')
    
    inArm1 = 1*tracking.zone.inArm1;
    inArm2 = 2*tracking.zone.inArm2;
    inArm3 = 3*tracking.zone.inArm3;
    inCenter = 4*tracking.zone.inCenter;
    inNoZone = 5*tracking.zone.inNoZone;
    transition = zeros(2,length(inArm1));
    
    for i=1:length(inArm1)
        if inArm1(i) > 0
            transition(1,i) = inArm1(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inArm2(i) > 0
            transition(1,i) = inArm2(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inArm3(i) > 0
            transition(1,i) = inArm3(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inCenter(i) > 0
            transition(1,i) = inCenter(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inNoZone(i) > 0
            transition(1,i) = inNoZone(i);
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
performance.paradigm = 'YMaze';
performance.timestamps = tracking.timestamps;

if saveMat
    
    save([basepath filesep filename '.Performance.Behavior.mat'],'performance');
end

end
