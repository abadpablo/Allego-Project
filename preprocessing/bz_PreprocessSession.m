function  bz_PreprocessSession(varargin)

%         bz_PreprocessSession(varargin)

%   Master function to run the basic pre-processing pipeline for an
%   individual sessions. Is based on sessionsPipeline.m but in this case
%   works on an individual session basis no in a folfer with multiple ones.
% 

% INPUTS
%   <options>       optional list of property-value pairs (see table below)
%   basepath        - Basepath for experiment. It contains all session
%                       folders. If not provided takes pwd.
%   analogCh       - List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
%   forceSum       - Force make folder summary (overwrite, if necessary). Default false.
%   cleanArtifacts - Remove artifacts from dat file (false by default). If
%                   true, remove artifacts from all Analog events. It also
%                   accepts a two rows cell with the analog channel
%                   (cleanArtifacts{1}) and the digital channels
%                   (cleanArtifacts{2}) to be used.
%   stateScore     - Run automatic brain state detection with SleepScoreMaster. Default true.
%   spikeSort      - Run automatic spike sorting using Kilosort. Default true.
%   getPos         - get tracking positions. Default true. 
%   pullData       - Path for raw data. Look for not analized session to copy to the main folder basepath. To do...
%   medianSubstr   - Perform median substraction in dat file before
%                       kilosort. Careful!! it would compromises dat file!
%                       (default false). If scalar, perform median
%                       substraction in those channels.
%
%  HISTORY: 
%     - Created based on sessionsPipeline: AntonioFR, 5/20

%  TO DO:
%   - Verify that data format and alysis output are compatible with CellExplorer
%   - Include Kilosort2 support
%   - Improve auto-clustering routine 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder); % by default, current folder
addParameter(p,'analogCh',[],@isnumeric);
addParameter(p,'stateScore',true,@islogical);
addParameter(p,'spikeSort',true,@islogical);
addParameter(p,'getPos',true,@islogical);
addParameter(p,'cleanArtifacts',false);
addParameter(p,'medianSubstr',false);
addParameter(p,'runSummary',false);
addParameter(p,'forceSum',false,@islogical);
% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

basepath = p.Results.basepath;
analogCh = p.Results.analogCh;
stateScore = p.Results.stateScore;
spikeSort = p.Results.spikeSort;
getPos = p.Results.getPos;
cleanArtifacts = p.Results.cleanArtifacts;
medianSubstr = p.Results.medianSubstr;
runSummary = p.Results.runSummary;
forceSum = p.Results.forceSum;

if ~exist('basepath') || isempty(basepath)
    basepath = uigetdir; % select folder
end
cd(basepath);

%% deals with xml.
if strcmp(basepath(end),filesep)
    basepath = basepath(1:end-1);
end
[~,basename] = fileparts(basepath);

disp('Check xml...');

actualPath = pwd;
prevPath = strsplit(actualPath,filesep);
length_prevPath = length(prevPath);

for i=1:length_prevPath-1
    prevPath_aux(i) = strcat(prevPath(i));
end

if length(prevPath_aux) == 1
    prevPath = strcat(prevPath(1));
elseif length(prevPath_aux) == 2
    prevPath = strcat(prevPath(1),'\',prevPath(2));
elseif length(prevPath_aux) == 3
    prevPath = strcat(prevPath(1),'\',prevPath(2),'\',prevPath(3));
elseif length(prevPath_aux) == 4
    prevPath = strcat(prevPath(1),'\',prevPath(2),'\',prevPath(3),'\',prevPath(4));
elseif length(prevPath_aux) == 5
    prevPath = strcat(prevPath(1),'\',prevPath(2),'\',prevPath(3),'\',prevPath(4),'\',prevPath(5));
end

prevPath = prevPath{1};
cd(prevPath);

if ~isempty(dir('global.xml'))
    xmlFile = dir('global.xml');
    cd(basepath)
    allpath = strsplit(genpath(basepath),';'); % all folders
    copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,basename,'.xml'));


elseif isempty(dir([basename '.xml'])) && isempty(dir('global.xml'))
    disp('No xml global file! Looking for it...');
    allpath = strsplit(genpath(basepath),';'); % all folders
    xmlFile = []; ii = 2;
    while isempty(xmlFile) && ii < size(allpath,2)
        disp(ii);
        cd(allpath{ii});
        xmlFile = dir('*global*.xml');
        ii = ii + 1;
    end
    if isempty(xmlFile)    
        [file, path] = uigetfile('*.xml','Select global xml file');
        copyfile(strcat(path,file),[basename '.xml']);
    else
        copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,basename,'.xml'));
    end
    cd(basepath);
end

%% Make SessionInfo
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  
save(strcat(basename,'.sessionInfo.mat'),'sessionInfo');
try
    session = sessionTemplate(pwd,'showGUI',false); % 
    save([basename '.session.mat'],'session');
catch
    warning('it seems that CellExplorer is not on your path');
end

%% Concatenate sessions
cd(basepath);
disp('Concatenate session folders...');
bz_ConcatenateDats(pwd,0,1);



%% Get analog and digital pulses
% if  ~isempty(analogCh)
%     [pulses] = bz_getAnalogPulses('analogCh',analogCh,'manualThr');
% end
% if ~isempty(dir('*digitalIn.mat'))
%     digitalIn = bz_getDigitalIn('all','fs',session.extracellular.sr); 
% end
%% Remove stimulation artifacts
if iscell(cleanArtifacts) || cleanArtifacts
    if iscell(cleanArtifacts)
        pulArtifacts_analog = bz_getAnalogPulses('analogCh',cleanArtifacts{1});
        pulArtifacts_dig = [];
        for ii = cleanArtifacts{2}
            disp(ii);
            pulArtifacts_dig = [pulArtifacts_dig; digitalIn.ints{ii}(:)];
        end
        pulArtifacts = sort([pulArtifacts_analog.timestamps(:); pulArtifacts_dig]);
    else
        pulArtifacts = pulses.timestamps(:);
    end
    cleanPulses(pulArtifacts);
end

%% Eliminate artifacts from linearTrack or TMaze

% sessionRemoveArtifacts()

%% Make LFP
if isempty(dir('*.lfp'))
    try 
        disp('Making LFP file ...');
        bz_LFPfromDat(pwd,'outFs',1250,'useGPU',false); % generating lfp, NOTE: THIS FUNCTION WILL GENERATE A SESSIONINFO FILE!! WE NEED TO FIX THIS
        
    catch
        disp('Problems with bz_LFPfromDat, resampling...');
%         sessionFile = dir('*session.mat*');
%         if ~isempty(sessionFile)
%             load(sessionFile.name);
%         end
        ResampleBinary(strcat(basename,'.dat'),strcat(basename,'.lfp'),...
            session.extracellular.nChannels,1, session.extracellular.sr/session.extracellular.srLfp);
    end
end

%% MEDIAN SUBS
if islogical(medianSubstr) && medianSubstr
%     removeNoiseFromDat(pwd);
    removeNoiseFromLFP(pwd,'keepDat',true);
elseif medianSubstr
    removeNoiseFromDat(pwd,'ch',medianSubstr);
end



%% Get brain states
if stateScore 
    try 
        if exist('pulses','var')
            SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
        else
            if ~exist([basename,'.SleepScoreLFP.LFP.mat'])
                SleepScoreMaster(pwd,'noPrompts',true); % try to sleep score
            end
        end
    catch
        disp('Problem with SleepScore skyping...');
    end
end

%% Kilosort concatenated sessions
if spikeSort
    disp('Spike sort concatenated sessions...');
    if  isempty(dir('*Kilosort*')) % if not kilosorted yet
    fprintf(' ** Kilosorting session...');
        if islogical(medianSubstr) && medianSubstr
            KiloSortWrapper_medianSubstraction;
            kilosortFolder = dir('*Kilosort*');
            try PhyAutoClustering(strcat(kilosortFolder.folder,'\',kilosortFolder.name)); % autoclustering
            catch
                warning('PhyAutoClustering not possible!!');
            end
            if exist('phyLink') && ~isempty(phyLink) % move phy link to
                kilosort_path = dir('*Kilosort*');
                try copyfile(phyLink, strcat(kilosort_path.name,filesep,'LaunchPhy')); % copy pulTime to kilosort folder
                end
            end
        else
            KiloSortWrapper;
            kilosortFolder = dir('*Kilosort*');
            try PhyAutoClustering(strcat(kilosortFolder.folder,'\',kilosortFolder.name)); % autoclustering
            catch
                warning('PhyAutoClustering not possible!!');
            end
            if exist('phyLink') && ~isempty(phyLink) % move phy link to
                kilosort_path = dir('*Kilosort*');
                try copyfile(phyLink, strcat(kilosort_path.name,filesep,'LaunchPhy')); % copy pulTime to kilosort folder
                end
            end
            
        end
    end
end
cd(basepath);
%% Get tracking positions 
if getPos
    getSessionTracking;
end

%% Summary
if runSummary
    if forceSum || (~isempty(dir('*Kilosort*')) && (isempty(dir('summ*')))) % is kilosorted but no summ
        disp('running summary analysis...');
%         sessionSummary;
        sessionSummary_abad;
    end
end


end

