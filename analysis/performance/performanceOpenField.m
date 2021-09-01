function performance = performanceOpenField(varargin)

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
    
performance.performance = [];
performance.transitions_out_nonzero = [];
performance.transitions_fr_out_nonzero = [];
performance.transitions_times_out_nonero = [];
performance.right_transitions_epoch = [];
performance.right_transitions_fr = [];
performance.wrong_transitions_epoch = [];
performance.wrong_transitions_fr = [];
performance.entrances = [];
performance.transitions_out = [];
performance.transitions_fr_out = [];
performance.transitions_times_out = [];
performance.frames_stem = [];
performance.frames_left = [];
performance.frames_right = [];
performance.times_stem = [];
performance.times_left = [];
performance.times_right = [];
performance.ballistic = [];
performance.doubted = [];
performance.ballistic_perc = [];
performance.doubted_perc = [];
performance.times = [];
performance.right_epoch = [];
performance.right_fr = [];
performance.righ_times = [];
performance.wrong_epoch = [];
performance.wrong_fr = [];
performance.wrong_times = [];
performance.rightTriades = [];
performance.wrongTriades = [];
performance.sampleVSchoice = [];

filename = strsplit(basepath,filesep);
filename = filename{end};

performance.folder = filename;
performance.paradigm = 'Open Field';
performance.timestamps = tracking.timestamps;

if saveMat
    
    save([basepath filesep filename '.Performance.Behavior.mat'],'performance');
end

end

end
