function performance = performanceObjectRecognition(varargin)

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

if isfield(digitalIn,'timestampsOn') && length(digitalIn.timestampsOn) >= 6 && useDigitalIn
    disp('Performance can be computed with digital Inputs')
    
    
else
    disp('Peformance can NOT be computed with digital Inputs')
    
    inObject1 = 1*tracking.zone.inObject1;
    inObject2 = 2*tracking.zone.inObject2;
    inOut = 3*tracking.zone.inOut;
    
    for i=1:length(inObject1)
        if inObject1(i) > 0
            transition(1,i) = inObject1(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inObject2(i) > 0
            transition(1,i) = inObject2(i);
            transition(2,i) = tracking.timestamps(i);
        elseif inOut(i) > 0
            transition(1,i) = inOut(i);
            transition(2,i) = tracking.timestamps(i);
        end
    end
    
%     transitions = analysisTransitions(transition,'tracking',tracking);   
    transitions = analysisTransitions_ObjectRecognition(transition,'tracking',tracking);
    
end
    
    
performance = transitions;

filename = strsplit(basepath,filesep);
filename = filename{end};

performance.folder = filename;
performance.paradigm = 'Object Recognition';
performance.timestamps = tracking.timestamps;

if saveMat
    
    save([basepath filesep filename '.Performance.Behavior.mat'],'performance');
end

end

