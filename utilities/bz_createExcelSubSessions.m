function [] = bz_createExcelSubSessions(listOfAnalysis,varargin)
%
% Creates excel file with the properties indicated in inputs (analyseSubSessions, spikes, ripples, theta-gamma, spikeTrain) 
%
% USAGE
%
%   [excel] = bz_createExcel(varargin)
%
% INPUTS
%   basePath       -(default: pwd) basePath for the recording file, in buzcode format:
%   analyzeSubSessions     
%   
%
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   position.x               - x position in cm/ normalize
%   position.y               - y position in cm/ normalize
%   timestamps      - in seconds, if Basler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%       only for OptiTrack
%   position.z 
%   orientation.x
%   orientation.y
%   orientation.z

%   HISTORY:
%     - Pablo Abad 2021

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'analyzeSubSessions',true, @islogical);
addParameter(p,'diffLFPs',true, @islogical); % To compute phase-locking with the channel where the amplitude of the unit is maximum and not with a reference channel for all units
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'getWaveformsFromDat',true,@islogical);
addParameter(p,'showWaveforms',true,@islogical);
addParameter(p,'firingMaps',[],@isstruct);
addParameter(p,'ripples',[],@isstruct);
addParameter(p,'rippleChannels',[],@isstruct);
addParameter(p,'spikeTrain',[],@isstruct);
addParameter(p,'behaviour',[],@isstruct);
addParameter(p,'performance',[], @isstruct);
addParameter(p,'thetaGamma',[], @isstruct);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'pathExcel',[],@isstr);
addParameter(p,'nameExcel',[],@isstr);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'PhaseLockingData',[],@isstruct);
addParameter(p,'PhaseLockingData_sg',[],@issstruct);
addParameter(p,'powerProfile_theta',[],@isstruct);
addParameter(p,'powerProfile_sg',[],@isstruct);
addParameter(p,'powerProfile_hg',[],@isstruct);
addParameter(p,'thetaFreq',[],@isnumeric);
addParameter(p,'sgFreq',[],@isnumeric);
addParameter(p,'hgFreq',[],@isnumeric);
addParameter(p,'hfoFreq',[],@isnumeric);

% LFP
addParameter(p,'phaseFreq',[4 12], @isnumeric);
% addParameter(p,'ampFreq',[30 80; 80 150; 150 200; 1 200], @isnumeric);
addParameter(p,'CFCPhaseAmp',[],@isstruct);
addParameter(p,'coherence_Shanks',[],@isstruct);
addParameter(p,'coherogram',[],@isstruct);
addParameter(p,'GMI',[],@isstruct);
addParameter(p,'PhaseAmpCouplingByAmp',[],@isstruct);
addParameter(p,'MI',[],@isstruct);


parse(p,varargin{:});
basepath = p.Results.basepath;
analyzeSubSessions = p.Results.analyzeSubSessions;
diffLFPs = p.Results.diffLFPs;
spikes = p.Results.spikes;
getWaveformsFromDat = p.Results.getWaveformsFromDat;
showWaveforms = p.Results.showWaveforms;
firingMaps = p.Results.firingMaps;
ripples = p.Results.ripples;
rippleChannels = p.Results.rippleChannels;
spikeTrain = p.Results.spikeTrain;
behaviour = p.Results.behaviour;
performance = p.Results.performance;
thetaGamma = p.Results.thetaGamma;
forceReload = p.Results.forceReload;
saveMat = p.Results.saveMat;
pathExcel = p.Results.pathExcel;
nameExcel = p.Results.nameExcel;
tracking = p.Results.tracking;
PhaseLockingData = p.Results.PhaseLockingData;
PhaseLockingData_sg = p.Results.PhaseLockingData_sg;
powerProfile_theta = p.Results.powerProfile_theta;
powerProfile_sg = p.Results.powerProfile_sg;
powerProfile_hg = p.Results.powerProfile_hg;
thetaFreq = p.Results.thetaFreq;
sgFreq = p.Results.sgFreq;
hgFreq = p.Results.hgFreq;
hfoFreq = p.Results.hfoFreq;

% LFP
phaseFreq = p.Results.phaseFreq;
% ampFreq = p.Results.ampFreq;
CFCPhaseAmp = p.Results.CFCPhaseAmp;
coherence_Shanks = p.Results.coherence_Shanks;
coherogram = p.Results.coherogram;
GMI = p.Results.GMI;
PhaseAmpCouplingByAmp = p.Results.PhaseAmpCouplingByAmp;
MI = p.Results.MI;


% ampFreq = [thetaFreq; sgFreq; hgFreq; hfoFreq];
ampFreq = [sgFreq; hgFreq; hfoFreq];

%% In case Excel already exists 
if ~isempty(dir([basepath filesep '*Excel.mat'])) || forceReload
    disp('Excel already detected! Loading file.');
    file = dir([basepath filesep '*Excel.Behavior.mat']);
    load(file.name);
    return
end

%% LOAD SESSIONINFO
sessionInfo = bz_getSessionInfo();
session = loadSession();

%% LOAD MERGEPOINTS
if ~isempty(dir([basepath filesep '*MergePoints.events.mat']))
    disp('Loading MergePoints...')
    file = dir([basepath filesep '*MergePoints.events.mat']);
    load(file.name)
end

%% 1 - LOAD SPIKES
if any(ismember(listOfAnalysis,'spikes'))
    if isempty(spikes)
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'showWaveforms',showWaveforms);
    end
end
%% 2 - RIPPLES
if any(ismember(listOfAnalysis,'ripples'))
    if isempty(ripples)
        if ~isempty(dir([basepath filesep '*ripples.SubSession.events.mat']))
            disp('Loading Ripples for SubSessions ...')
            file = dir([basepath filesep '*ripples.SubSession.events.mat']);
            load(file.name)
        else
            warning('It is not possible to load Ripples SubSessions ...')
        end
    end 

%     if isempty(rippleChannels)
%         if ~isempty(dir([basepath filesep '*.channelInfo.SubSession.ripples.mat']))
%             disp('Loading rippleChannels for SubSessions ...')
%             file = dir([basepath filesep '*channelInfo.SubSession.ripples.mat']);
%             load(file.name)
%         else
%             warning('It is not possible to load rippleChannel SubSessions...');
%         end
%     end
    
    if isempty(rippleChannels)
        if ~isempty(dir([basepath filesep '*.channelInfo.ripples.mat']))
            disp('Loading rippleChannels  ...')
            file = dir([basepath filesep '*channelInfo.ripples.mat']);
            load(file.name)
        else
            warning('It is not possible to load rippleChannel SubSessions...');
        end
    end  
end
%% 3 - POWER PROFILE
if any(ismember(listOfAnalysis,'powerSpectrumProfile'))
    if isempty(powerProfile_theta)
        if ~isempty(dir([basepath filesep '*PowerSpectrumProfile_',num2str(thetaFreq(1)),'_',num2str(thetaFreq(2)),'.SubSession.channelinfo.mat']))
            disp('Loading Power Profile Theta for SubSessions...')
            file = dir([basepath filesep '*PowerSpectrumProfile_',num2str(thetaFreq(1)),'_',num2str(thetaFreq(2)),'.SubSession.channelInfo.mat'])
            load(file.name);
        end
    end
    
    if isempty(powerProfile_sg)
        if ~isempty(dir([basepath filesep '*PowerSpectrumProfile_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.channelinfo.mat']))
            disp('Loading Power Profile Slow Gamma for SubSessions ...')
            file = dir([basepath filesep '*PowerSpectrumProfile_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.channelinfo.mat'])
            load(file.name);
        end
    end
    
    if isempty(powerProfile_hg)
        if ~isempty(dir([basepath filesep '*PowerSpectrumProfile_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.channelinfo.mat']))
            disp('Loading Power Profile High Gamma for SubSessions ...')
            file = dir([basepath filesep '*PowerSpectrumProfile_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.channelinfo.mat'])
            load(file.name);
        end
    end
end

%% 3 - PHASELOCKING THETA-GAMMA
if any(ismember(listOfAnalysis,'thetaModulation'))
    if ~isempty(dir([basepath filesep '*PhaseLockingData.SubSession.cellinfo.mat']))
        disp('Loading SubSession Phase Locking Theta...')
        file = dir([basepath filesep '*PhaseLockingData.SubSession.cellinfo.mat'])
        load(file.name)
    end

    if ~isempty(dir([basepath filesep '*PhaseLockingData_sg.SubSession.cellinfo.mat']))
        disp('Loading SubSession Phase Locking Slow Gamma...')
        file = dir([basepath filesep '*PhaseLockingData_sg.SubSession.cellinfo.mat'])
        load(file.name) 
    end
end

%% 4 - LOAD TRACKING 
if any(ismember(listOfAnalysis,'behaviour'))
    if isempty(tracking)
        tracking = getSessionTracking();
    end
    tracking_sr = tracking.samplingRate(1);

    %  LOAD BEHAVIOUR
    if isempty(behaviour)
        behaviour = getSessionBehaviour_v2();
    end
end

%% 4 - LOAD FIRINGMAPS
if any(ismember(listOfAnalysis,'placeCells'))
    if isempty(firingMaps)
        if ~isempty(dir([basepath filesep '*firingMapsAvg.cellinfo.mat']))
            disp('Firing Maps already detected! Loading file.');
            file = dir([basepath filesep '*firingMapsAvg.cellinfo.mat']);
            load(file.name);
        end
    end
end

%% 5 - LOAD PERFORMANCE
if any(ismember(listOfAnalysis,'performance'))
    if isempty(performance)
        try
            performance = getSessionPerformance('tracking',tracking); 
        catch     
            warning('It is not possible to load Performance...');
            performance = [];
        end
    end
end

%% 6 - SPIKETRAIN
if any(ismember(listOfAnalysis,'spikeTrain'))
    if isempty(spikeTrain)
        spikeTrain = bz_SpikeTrain(spikes,'analyzeSubSessions',true);
    end 
end

%% LOAD ALL LFP VARIABLES
if any(ismember(listOfAnalysis,'lfp_analysis'))

    % Coherence_Shanks
    if isempty(coherence_Shanks)
       if ~isempty(dir([basepath filesep sessionInfo.FileName '.Coherence_Shanks.SubSession.lfp.mat']))
           disp('Coherence_Shanks.SubSession.lfp.mat already detected. Loading file !');
           file = dir([basepath filesep sessionInfo.FileName '.Coherence_Shanks.SubSession.lfp.mat']);
           load(file.name);
       end   
    end

    % Coherogram
    if isempty(coherogram)
        if ~isempty(dir([basepath filesep sessionInfo.FileName '.Coherogram.SubSession.lfp.mat']))
            disp('Coherogram.SubSession.lfp.mat already detected. Loading file !');
            file = dir([basepath filesep sessionInfo.FileName '.Coherogram.SubSession.lfp.mat']);
            load(file.name);
        end
    end

    % CFCPhaseAmp
    if isempty(CFCPhaseAmp)
        for i=1:size(ampFreq,1)
            if ~isempty(dir([basepath filesep sessionInfo.FileName '.CFCPhaseAmp_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat']))
                disp(['CFCPhaseAmp_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat already detected !. Loading file...' ]);
                file = dir([basepath filesep sessionInfo.FileName '.CFCPhaseAmp_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat']);
%                 CFCPhaseAmp{i} = load(file.name);
%                 CFCPhaseAmp{i}.ampFreq = ampFreq(i,:);
                load(file.name);
            end       
        end 
    end
    
    count = 0;
    for i=1:size(ampFreq,1)
        count = count+1;
        switch count
            case 1
                if exist('CFCPhaseAmp_lg','var')
                    CFCPhaseAmp_lg.ampFreq = ampFreq(i,:);
                end
            case 2
                if exist('CFCPhaseAmp_hg','var')
                    CFCPhaseAmp_hg.ampFreq = ampFreq(i,:);
                end
            case 3
                if exist('CFCPhaseAmp_hfo','var')
                    CFCPhaseAmp_hfo.ampFreq = ampFreq(i,:);
                end
        end
    end
    
    
        

    % PhaseAmpCouplingByAmp
    if isempty(PhaseAmpCouplingByAmp)
        for i=1:size(ampFreq,1)
            if ~isempty(dir([basepath filesep sessionInfo.FileName '.PhaseAmpCouplingbyAmp_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat']))
                disp(['PhaseAmpCouplingByAmp',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat already detected! Loading file...']);
                file = dir([basepath filesep sessionInfo.FileName '.PhaseAmpCouplingByAmp_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat']);
%                 PhaseAmpCouplingByAmp{i} = load(file.name);  
                load(file.name)
            end           
        end
    end
    
    count = 0;
    for i=1:size(ampFreq,1)
        count = count+1;
        switch count
            case 1
                if exist('PhaseAmpCoupling_lg','var')
                    PhaseAmpCoupling_lg.ampFreq = ampFreq(i,:);
                end
            case 2
                if exist('PhaseAmpCoupling_hg','var')
                    PhaseAmpCoupling_hg.ampFreq = ampFreq(i,:);
                end
            case 3
                if exist('PhaseAmpCoupling_hfo','var')
                    PhaseAmpCoupling_hfo.ampFreq = ampFreq(i,:);
                end
        end
    end

    % GMI
    if isempty(GMI)
        for i=1:size(ampFreq,1)
            if ~isempty(dir([basepath filesep sessionInfo.FileName '.GMI_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat']))
                disp(['GMI_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat already detected !. Loading file...' ]);
                file = dir([basepath filesep sessionInfo.FileName '.GMI_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat']);
                load(file.name);
            end            
        end  
    end
    
    count = 0;
    for i=1:size(ampFreq,1)
        count = count+1;
        switch count
            case 1
                if exist('GMI_lg','var')
                    GMI_lg.ampFreq = ampFreq(i,:);
                end
            case 2
                if exist('GMI_hg','var')
                    GMI_hg.ampFreq = ampFreq(i,:);
                end
            case 3
                if exist('GMI_hfo','var')
                    GMI_hfo.ampFreq = ampFreq(i,:);
                end
        end
    end
    
    
    % MI
    if isempty(MI)
        for i=1:size(ampFreq,1)
            if ~isempty(dir([basepath filesep sessionInfo.FileName '.GMI_Tort_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat']))
                disp(['MI Tort_', num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat already detected! Loading file...']);
                file = dir([basepath filesep sessionInfo.FileName '.GMI_Tort_',num2str(ampFreq(i,1)),'_',num2str(ampFreq(i,end)),'.SubSession.lfp.mat']);
                load(file.name);
            end       
        end    
    end
    
    count = 0;
    for i=1:size(ampFreq,1)
        count = count+1;
        switch count
            case 1
                if exist('MI_lg','var')
                    MI_lg.ampFreq = ampFreq(i,:);
                end
            case 2
                if exist('MI_hg','var')
                    MI_hg.ampFreq = ampFreq(i,:);
                end
            case 3
                if exist('MI_hfo','var')
                    MI_hfo.ampFreq = ampFreq(i,:);
                end
        end
    end
    
end

%% EXTRACT THE VARIABLES RELATED TO POWER SPECTRUM AND COHERENCE
if ~isempty(coherence_Shanks)
    coherence_Shanks_mean = bz_extractFrequencyBands(coherence_Shanks,'variable','coherence_Shanks','analyzeSubSessions',analyzeSubSessions);
end
if ~isempty(coherogram)
    coherogram_mean = bz_extractFrequencyBands(coherogram,'variable','coherogram','analyzeSubSessions',analyzeSubSessions);
end

%% ===================================================================================
% --------------- CREATING EXCEL VARIABLES  ----------------------------------
% ====================================================================================

%% SPIKES
if any(ismember(listOfAnalysis,'spikes'))
    sheet_spikes = 'spikes';
    Header_spikes = {'foldername' 'shankID' 'ch' 'totalSpikes' 'amplitude' }
    for j=1:spikes.numcells
        Output_spikes{j} = {spikes.basename spikes.shankID(j) spikes.maxWaveformCh(j) spikes.total(j) mean(spikes.amplitudes{j})};
    end
end
%% RIPPLES
if any(ismember(listOfAnalysis,'ripples'))
    sheet_ripples = 'ripples';
    Header_ripples = {'foldername' 'rippleChannel' 'SharpWaveChannel' 'NoiseChannel' 'numRipples' 'peakFrequency' 'peakAmplitude' 'duration' };

%     for i=1:length(ripples)
%         Output_ripples{i} = {ripples{i}.foldername rippleChannels{i}.Ripple_Channel rippleChannels{i}.Sharpwave_Channel rippleChannels{i}.Noise_Channel size(ripples{i}.timestamps,1) mean(ripples{i}.data.peakFrequency) mean(ripples{i}.data.peakAmplitude) mean(ripples{i}.data.duration)};
%     end
    for i=1:length(ripples)
        if isfield(ripples{i},'timestamps')
            if ~isempty(ripples{i}.timestamps) && ~isempty(ripples{i}.data)
                Output_ripples{i} = {ripples{i}.foldername rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel size(ripples{i}.timestamps,1) mean(ripples{i}.data.peakFrequency) mean(ripples{i}.data.peakAmplitude) mean(ripples{i}.data.duration)};
            else
                Output_ripples{i} = {ripples{i}.foldername rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel [] [] [] []};
            end
                
        else
            Output_ripples{i} = {ripples{i}.foldername rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel rippleChannels.Noise_Channel [] [] [] []};
        end
    end
end

%% LFP
if any(ismember(listOfAnalysis,'lfp_analysis'))
    % CFCPhaseAmp
    sheet_CFCPhaseAmp = 'CFCPhaseAmp';
    Header_CFCPhaseAmp = {'foldername','Channel1','Channel2','freqAmp(1)','freqAmp(2)','maxValue','phase_maxValue','amp_maxValue'};

    for i=1:size(ampFreq,1)
        for j = 1:length(MergePoints.foldernames)
            if i == 1
                Output_CFC{i}{j} = {MergePoints.foldernames{j} CFCPhaseAmp_lg.(MergePoints.foldernames{j}).phaseCh CFCPhaseAmp_lg.(MergePoints.foldernames{j}).ampChans CFCPhaseAmp_lg.ampFreq(1) CFCPhaseAmp_lg.ampFreq(2) CFCPhaseAmp_lg.(MergePoints.foldernames{j}).maxValue CFCPhaseAmp_lg.(MergePoints.foldernames{j}).phase_maxValue CFCPhaseAmp_lg.(MergePoints.foldernames{j}).amp_maxValue};                
            elseif i == 2
                Output_CFC{i}{j} = {MergePoints.foldernames{j} CFCPhaseAmp_hg.(MergePoints.foldernames{j}).phaseCh CFCPhaseAmp_hg.(MergePoints.foldernames{j}).ampChans CFCPhaseAmp_hg.ampFreq(1) CFCPhaseAmp_hg.ampFreq(2) CFCPhaseAmp_hg.(MergePoints.foldernames{j}).maxValue CFCPhaseAmp_hg.(MergePoints.foldernames{j}).phase_maxValue CFCPhaseAmp_hg.(MergePoints.foldernames{j}).amp_maxValue};                
            elseif i == 3
                Output_CFC{i}{j} = {MergePoints.foldernames{j} CFCPhaseAmp_hfo.(MergePoints.foldernames{j}).phaseCh CFCPhaseAmp_hfo.(MergePoints.foldernames{j}).ampChans CFCPhaseAmp_hfo.ampFreq(1) CFCPhaseAmp_hfo.ampFreq(2) CFCPhaseAmp_hfo.(MergePoints.foldernames{j}).maxValue CFCPhaseAmp_hfo.(MergePoints.foldernames{j}).phase_maxValue CFCPhaseAmp_hfo.(MergePoints.foldernames{j}).amp_maxValue};                
            end
        end
    end

    if ~isempty(coherence_Shanks)
        % Coherence_Shanks
        sheet_coherenceShanks = 'coherence_Shanks';
        Header_coherenceShanks = {'foldername', 'shankID 1', 'Channel 1', 'shankID 2', 'Channel 2', ...
            'Coherence Delta (1-3Hz)', 'Coherence Theta (4-12Hz)', 'Coherence Alpha (13-16Hz)', 'Coherence Beta (17-29Hz)', 'Coherence LG (30-65Hz)', 'Coherence HG (66-130)', 'Coherence HFO (150 185)',...
            'Phase Delta (1-3Hz)', 'Phase Theta (4-12Hz)', 'Phase Alpha (13-16Hz)', 'Phase Beta (17-29Hz)', 'Phase LG (30-65Hz)', 'Phase HG (66-130)', 'Phase HFO (150 185)',...
            'Power Delta (1-3Hz) Ch1 ', 'Power Theta  (4-12Hz) Ch1 ', 'Power Alpha (13-16Hz) Ch1 ', 'Power Beta (17-29Hz) Ch1 ', 'Power LG (30-65Hz) Ch1 ', 'Power HG (66-130Hz) Ch1 ', 'Power HFO (150-185Hz) Ch1 ',...
            'Power Delta (1-3Hz) Ch2 ', 'Power Theta  (4-12Hz) Ch2 ', 'Power Alpha (13-16Hz) Ch2 ', 'Power Beta (17-29Hz) Ch2 ', 'Power LG (30-65Hz) Ch2 ', 'Power HG (66-130Hz) Ch2 ', 'Power HFO (150-185Hz) Ch2 '};

        for i=1:length(coherence_Shanks_mean.coherence_mean)
            for j=1:length(coherence_Shanks_mean.coherence_mean{i})
                Output_coherenceShanks{i}{j} = {sessionInfo.FileName i coherence_Shanks.ch(i) j coherence_Shanks.ch(j) ...
                    coherence_Shanks_mean.coherence_mean{i}{j}(1) coherence_Shanks_mean.coherence_mean{i}{j}(2) coherence_Shanks_mean.coherence_mean{i}{j}(3) coherence_Shanks_mean.coherence_mean{i}{j}(4) coherence_Shanks_mean.coherence_mean{i}{j}(5) coherence_Shanks_mean.coherence_mean{i}{j}(6) coherence_Shanks_mean.coherence_mean{i}{j}(7)...
                    coherence_Shanks_mean.phase_mean{i}{j}(1) coherence_Shanks_mean.phase_mean{i}{j}(2) coherence_Shanks_mean.phase_mean{i}{j}(3) coherence_Shanks_mean.phase_mean{i}{j}(4) coherence_Shanks_mean.phase_mean{i}{j}(5) coherence_Shanks_mean.phase_mean{i}{j}(6) coherence_Shanks_mean.phase_mean{i}{j}(7)...
                    coherence_Shanks_mean.S1_mean{i}{j}(1) coherence_Shanks_mean.S1_mean{i}{j}(2) coherence_Shanks_mean.S1_mean{i}{j}(3) coherence_Shanks_mean.S1_mean{i}{j}(4) coherence_Shanks_mean.S1_mean{i}{j}(5) coherence_Shanks_mean.S1_mean{i}{j}(6) coherence_Shanks_mean.S1_mean{i}{j}(7)...
                    coherence_Shanks_mean.S2_mean{i}{j}(1) coherence_Shanks_mean.S2_mean{i}{j}(2) coherence_Shanks_mean.S2_mean{i}{j}(3) coherence_Shanks_mean.S2_mean{i}{j}(4) coherence_Shanks_mean.S2_mean{i}{j}(5) coherence_Shanks_mean.S2_mean{i}{j}(6) coherence_Shanks_mean.S2_mean{i}{j}(7)};
            end
        end
    end

    % Coherogram
    sheet_coherogram = 'coherogram';
    Header_coherogram = {'foldername' , 'Channel 1', 'Channel 2', ...
        'Coherence Delta (1-3Hz)', 'Coherence Theta (4-12Hz)', 'Coherence Alpha (13-16Hz)', 'Coherence Beta (17-29Hz)', 'Coherence LG (30-65Hz)', 'Coherence HG (66-130)', 'Coherence HFO (150 185)',...
        'Phase Delta (1-3Hz)', 'Phase Theta (4-12Hz)', 'Phase Alpha (13-16Hz)', 'Phase Beta (17-29Hz)', 'Phase LG (30-65Hz)', 'Phase HG (66-130)', 'Phase HFO (150 185)',...
        'Power Delta (1-3Hz) Ch1 ', 'Power Theta  (4-12Hz) Ch1 ', 'Power Alpha (13-16Hz) Ch1 ', 'Power Beta (17-29Hz) Ch1 ', 'Power LG (30-65Hz) Ch1 ', 'Power HG (66-130Hz) Ch1 ', 'Power HFO (150-185Hz) Ch1 ',...
        'Power Delta (1-3Hz) Ch2 ', 'Power Theta  (4-12Hz) Ch2 ', 'Power Alpha (13-16Hz) Ch2 ', 'Power Beta (17-29Hz) Ch2 ', 'Power LG (30-65Hz) Ch2 ', 'Power HG (66-130Hz) Ch2 ', 'Power HFO (150-185Hz) Ch2 '};
    for i = 1:length(fields(coherogram_mean))
        Output_coherogram{i} = {MergePoints.foldernames{i} coherogram_mean.(MergePoints.foldernames{i}).ch1 coherogram_mean.(MergePoints.foldernames{i}).ch2...
            coherogram_mean.(MergePoints.foldernames{i}).coherence_mean(1) coherogram_mean.(MergePoints.foldernames{i}).coherence_mean(2) coherogram_mean.(MergePoints.foldernames{i}).coherence_mean(3) coherogram_mean.(MergePoints.foldernames{i}).coherence_mean(4) coherogram_mean.(MergePoints.foldernames{i}).coherence_mean(5) coherogram_mean.(MergePoints.foldernames{i}).coherence_mean(6) coherogram_mean.(MergePoints.foldernames{i}).coherence_mean(7)...
            coherogram_mean.(MergePoints.foldernames{i}).phase_mean(1) coherogram_mean.(MergePoints.foldernames{i}).phase_mean(2) coherogram_mean.(MergePoints.foldernames{i}).phase_mean(3) coherogram_mean.(MergePoints.foldernames{i}).phase_mean(4) coherogram_mean.(MergePoints.foldernames{i}).phase_mean(5) coherogram_mean.(MergePoints.foldernames{i}).phase_mean(6) coherogram_mean.(MergePoints.foldernames{i}).phase_mean(7)...
            coherogram_mean.(MergePoints.foldernames{i}).S1_mean(1) coherogram_mean.(MergePoints.foldernames{i}).S1_mean(2) coherogram_mean.(MergePoints.foldernames{i}).S1_mean(3) coherogram_mean.(MergePoints.foldernames{i}).S1_mean(4) coherogram_mean.(MergePoints.foldernames{i}).S1_mean(5) coherogram_mean.(MergePoints.foldernames{i}).S1_mean(6) coherogram_mean.(MergePoints.foldernames{i}).S1_mean(7)...
            coherogram_mean.(MergePoints.foldernames{i}).S2_mean(1) coherogram_mean.(MergePoints.foldernames{i}).S2_mean(2) coherogram_mean.(MergePoints.foldernames{i}).S2_mean(3) coherogram_mean.(MergePoints.foldernames{i}).S2_mean(4) coherogram_mean.(MergePoints.foldernames{i}).S2_mean(5) coherogram_mean.(MergePoints.foldernames{i}).S2_mean(6) coherogram_mean.(MergePoints.foldernames{i}).S2_mean(7)};
    end

    % GMI
    shanks = sessionInfo.AnatGrps;
    sheet_GMI = 'GMI';
    Header_GMI = {'foldername' 'shankID' 'Channel' 'ampfreq(1)' 'ampfreq(2)' 'GMI' 'GFMI' };
    numChan = length(GMI_lg.(MergePoints.foldernames{1}).Channels.ampCh);
    if numChan == 1
        for i=1:length(shanks)
            a = any(ismember(shanks(i).Channels,GMI_lg.(MergePoints.foldernames{1}).Channels.ampCh));
            if a == 1
                numShank = i;
            end
        end
    end
    
    for i=1:size(ampFreq,1)
        for j = 1:length(MergePoints.foldernames)
            if numChan == 1
                if i == 1
                    Output_GMI{i}{j} = {GMI_lg.(MergePoints.foldernames{j}).folder numShank GMI_lg.(MergePoints.foldernames{j}).Channels.ampCh GMI_lg.ampFreq(1) GMI_lg.ampFreq(2) GMI_lg.(MergePoints.foldernames{j}).GMI{1} GMI_lg.(MergePoints.foldernames{j}).GFMI{1}};
                elseif i == 2
                    Output_GMI{i}{j} = {GMI_hg.(MergePoints.foldernames{j}).folder numShank GMI_hg.(MergePoints.foldernames{j}).Channels.ampCh GMI_hg.ampFreq(1) GMI_hg.ampFreq(2) GMI_hg.(MergePoints.foldernames{j}).GMI{1} GMI_hg.(MergePoints.foldernames{j}).GFMI{1}};
                elseif i == 3
                    Output_GMI{i}{j} = {GMI_hfo.(MergePoints.foldernames{j}).folder numShank GMI_hfo.(MergePoints.foldernames{j}).Channels.ampCh GMI_hfo.ampFreq(1) GMI_hfo.ampFreq(2) GMI_hfo.(MergePoints.foldernames{j}).GMI{1} GMI_hfo.(MergePoints.foldernames{j}).GFMI{1}};
                end
            else
                for j = 1:size(shanks,2)
                    for k = 1:length(shanks(j).Channels)
                        Output_GMI{i}{j}{k} = {sessionInfo.FileName j shanks(j).Channels(k) GMI{i}.ampFreq(1) GMI{i}.ampFreq(end) GMI{i}.GMI.GMI{shanks(j).Channels(k)+1} GMI{i}.GMI.GFMI{shanks(j).Channels(k)+1}};
                    end
                end
            end
        end
    end

    % PhaseAmpCouplingByAmp
%     sheet_phaseAmpCouplingByAmp = 'PhaseAmpCouplingByAmp';
%     Header_phaseAmpCouplingByAmp = {'foldername' 'phaseCh' 'ampCh'};

    % MI
    shanks = sessionInfo.AnatGrps;
    sheet_MI = 'MI';
    Header_MI = {'foldername' 'shankID' 'Channel' 'ampfreq(1)' 'ampfreq(2)' 'GMI'};
    numChan = length(MI_lg.(MergePoints.foldernames{1}).Channels.ampCh);
    if numChan == 1
        for i=1:length(shanks)
            a = any(ismember(shanks(i).Channels,MI_lg.(MergePoints.foldernames{1}).Channels.ampCh));
            if a == 1
                numShank = i;
            end
        end
    end
    
    for i=1:size(ampFreq,1)
        for j = 1:length(MergePoints.foldernames)
            if numChan == 1
                if i == 1
                    Output_MI{i}{j} = {MI_lg.(MergePoints.foldernames{j}).foldername numShank MI_lg.(MergePoints.foldernames{j}).Channels.ampCh MI_lg.ampFreq(1) MI_lg.ampFreq(2) MI_lg.(MergePoints.foldernames{j}).MI{1} };
                elseif i == 2
                    Output_MI{i}{j} = {MI_hg.(MergePoints.foldernames{j}).foldername numShank MI_hg.(MergePoints.foldernames{j}).Channels.ampCh MI_hg.ampFreq(1) MI_hg.ampFreq(2) MI_hg.(MergePoints.foldernames{j}).MI{1} };
                elseif i == 3
                    Output_MI{i}{j} = {MI_hfo.(MergePoints.foldernames{j}).foldername numShank MI_hfo.(MergePoints.foldernames{j}).Channels.ampCh MI_hfo.ampFreq(1) MI_hfo.ampFreq(2) MI_hfo.(MergePoints.foldernames{j}).MI{1} };
                end
            else
                for j = 1:size(shanks,2)
                    for k = 1:length(shanks(j).Channels)
                        Output_MI{i}{j}{k} = {sessionInfo.FileName j shanks(j).Channels(k) GMI{i}.ampFreq(1) GMI{i}.ampFreq(end) GMI{i}.GMI.GMI{shanks(j).Channels(k)+1} GMI{i}.GMI.GFMI{shanks(j).Channels(k)+1}};
                    end
                end
            end
        end
    end
    
    
end

%% FIRING MAPS

if any(ismember(listOfAnalysis,'placeCells'))
    
    sheet_firingMaps = 'firingMaps';
    Header_firingMaps = {'foldername' 'paradigm' 'shankID' 'ch' 'peak' 'mean' 'specificity' 'm' 'r' 'mode' 'k' 'spatialCorrelation_r' 'spatialCorrelation_p' 'spatialCorrelation_sc'  ...
                            'skaggs_bitsPerSec' 'skaggs_bitsPerSpike' 'skaggs_bitsPerSec_Uns' 'skaggs_bitsPerSpike_Uns' ...
                            'firingField_numFF' 'firingField_areaFF' 'firingField_areaTotFF' 'firingField_patchs' 'firingField_patchsArea' 'firingField_patchsAreaTot' 'firingField_FFArevspatchAr' 'firingField_TotSizeRat' 'firingField_MaxFfre'...
                            'borderIndex_west' 'borderIndex_east' 'borderIndex_north' 'borderIndex_south' 'borderIndex_maxBorder'};

    for i=1:length(firingMaps)
        for j=1:spikes.numcells
            Output_firingMaps{i}{j} = {tracking.folders{i} behaviour.description{i} spikes.shankID(j) spikes.maxWaveformCh(j) ... %ch
                                        firingMaps.stats{j}{1}.peak firingMaps.stats{j}{1}.mean firingMaps.stats{j}{1}.specificity firingMaps.stats{j}{1}.m firingMaps.stats{j}{1}.r firingMaps.stats{j}{1}.mode firingMaps.stats{j}{1}.k ... %k
                                        firingMaps.stats{j}{1}.spatialCorr.valor_r firingMaps.stats{j}{1}.spatialCorr.valor_p firingMaps.stats{j}{1}.spatialCorr.sc...
                                        firingMaps.stats{j}{1}.skaggs.bitsPerSec firingMaps.stats{j}{1}.skaggs.bitsPerSpike firingMaps.stats{j}{1}.skaggs.bitsPerSec_Uns firingMaps.stats{j}{1}.skaggs.bitsPerSpike_Uns ...
                                        firingMaps.stats{j}{1}.firingField.numFF firingMaps.stats{j}{1}.firingField.areaFF firingMaps.stats{j}{1}.firingField.areaTotFF firingMaps.stats{j}{1}.firingField.patchs firingMaps.stats{j}{1}.firingField.patchsArea ...
                                        firingMaps.stats{j}{1}.firingField.patchsAreaTot firingMaps.stats{j}{1}.firingField.FFArevspatchAr firingMaps.stats{j}{1}.firingField.TotSizeRat firingMaps.stats{j}{1}.firingField.MaxFfre...
                                        firingMaps.stats{j}{1}.borderIndex.west firingMaps.stats{j}{1}.borderIndex.east firingMaps.stats{j}{1}.borderIndex.north firingMaps.stats{j}{1}.borderIndex.south firingMaps.stats{j}{1}.borderIndex.maxBorder};                                
        end

    end
end

%% SPIKE TRAIN
if any(ismember(listOfAnalysis,'spikeTrain'))
    
    sheet_spikeTrain = 'spikeTrain';
    Header_spikeTrain = {'foldername' 'shankID' 'ch' 'expected' 'lb' 'ub' 'RPV' 'meanF' 'MeanAutoc' 'ProbBurst' 'BurstIndex' 'BurstIndexA' 'num_burst' 'n_spike' 'prom_isi' 'prom_duration' 'numBurstNorm' 'TMI' ...
                            'TMI2' 'TMI3' 'S_theta' 'Th_Peak' 'S_Gamma' 'S_Gamma2' 'G_Peak' 'thet_Mod' 'BackG' 'S_thetaFFT' 'Th_PeakFF' 'S_GammaFFT' 'G_PeakFF' 'BackGFF' };
    field = fields(spikeTrain);
    for i=1:length(fields(spikeTrain))
        for j=1:spikes.numcells
            Output_spikeTrain{i}{j} = {field{i} spikes.shankID(j) spikes.maxWaveformCh(j) spikeTrain.(field{i}){j}.expected spikeTrain.(field{i}){j}.lb spikeTrain.(field{i}){j}.ub spikeTrain.(field{i}){j}.RPV spikeTrain.(field{i}){j}.meanF spikeTrain.(field{i}){j}.MeanAutoc ... %MeanAutoc
                                        spikeTrain.(field{i}){j}.ProbBurst spikeTrain.(field{i}){j}.BurstIndex spikeTrain.(field{i}){j}.BurstIndexA spikeTrain.(field{i}){j}.num_burst spikeTrain.(field{i}){j}.n_spike spikeTrain.(field{i}){j}.prom_isi... %prom_isi
                                        spikeTrain.(field{i}){j}.prom_duration spikeTrain.(field{i}){j}.numBurstNorm spikeTrain.(field{i}){j}.TMI spikeTrain.(field{i}){j}.TMI2 spikeTrain.(field{i}){j}.TMI3 spikeTrain.(field{i}){j}.S_theta... %S_theta
                                        spikeTrain.(field{i}){j}.Th_Peak spikeTrain.(field{i}){j}.S_Gamma spikeTrain.(field{i}){j}.S_Gamma2 spikeTrain.(field{i}){j}.G_Peak spikeTrain.(field{i}){j}.thet_Mod spikeTrain.(field{i}){j}.BackG... %BackG
                                        spikeTrain.(field{i}){j}.S_thetaFFT spikeTrain.(field{i}){j}.Th_PeakFF spikeTrain.(field{i}){j}.S_GammaFFT spikeTrain.(field{i}){j}.G_PeakFF spikeTrain.(field{i}){j}.BackGFF};
        end
    end
end


%% PERFORMANCE
if any(ismember(listOfAnalysis,'performance'))
    sheet_performance = 'performance';
    Header_performance = {'foldername' 'paradigm' 'score'  'entrances_center' 'entrances_left' 'entrances_right' 'entrances_stem' ...
                            'times_center' 'times_left' 'times_right' 'times_stem'...
                            'group1 trials' 'group2 trials' 'group3 trials' 'group4 trials'  'ballistic' 'doubted' 'ballistic_perc' 'doubted_perc' };
                        
    for i=1:length(performance.score)
        Output_performance{i} = {performance.folder{i} performance.paradigm{i} performance.score(i) performance.entrances{i}.center performance.entrances{i}.left performance.entrances{i}.right performance.entrances{i}.stem...
                                    performance.times{i}.center performance.times{i}.left performance.times{i}.right performance.times{i}.stem ...
                                    performance.sampleVSchoice{i}.group1.trials performance.sampleVSchoice{i}.group2.trials performance.sampleVSchoice{i}.group3.trials performance.sampleVSchoice{i}.group4.trials ...
                                    performance.ballistic{i} performance.doubted{i} performance.ballistic_perc{i} performance.doubted_perc{i}};
    end

end
%% PHASELOCKINGDATA (THETA)
if any(ismember(listOfAnalysis,'thetaModulation'))
    sheet_phaseLockingTheta = 'phaseLockingTheta';
    Header_phaseLockingTheta = {'folder' 'shankID' 'ch' 'm' 'r' 'k' 'p' 'mode' };
    
    for j=1:spikes.numcells
        Output_phaseLockingTheta{j} = {PhaseLockingData{j}.sessionName spikes.shankID(j) spikes.maxWaveformCh(j) PhaseLockingData{j}.phasestats.m PhaseLockingData{j}.phasestats.r PhaseLockingData{j}.phasestats.k ...
                                        PhaseLockingData{j}.phasestats.p PhaseLockingData{j}.phasestats.mode};

    end


    %% PHASELOCKINGDATA (SLOW GAMMA)
    sheet_phaseLockingSG = 'phaseLockingSG';
    Header_phaseLockingSG = {'folder' 'shankID' 'ch' 'm' 'r' 'k' 'p' 'mode'};

    for j=1:spikes.numcells
        Output_phaseLockingSG{j} = {PhaseLockingData_sg{j}.sessionName spikes.shankID(j) spikes.maxWaveformCh(j) PhaseLockingData_sg{j}.phasestats.m PhaseLockingData_sg{j}.phasestats.r PhaseLockingData_sg{j}.phasestats.k ...
                                        PhaseLockingData_sg{j}.phasestats.p PhaseLockingData_sg{j}.phasestats.mode};
    end
end

%% POWER SPECTRUM PROFILE THETA
if any(ismember(listOfAnalysis,'powerSpectrumProfile'))
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
    for i=1:length(sessionInfo.channels)
        for j=1:length(sessionInfo.AnatGrps)
            if ismember(sessionInfo.channels(i),sessionInfo.AnatGrps(j).Channels)
                shank_ch(i) = j;
            end
        end
    end

    sheet_powerProfile_theta = 'powerProfile_theta';
    Header_powerProfile_theta = {'folder' 'ch' 'shankID' 'mean' 'std' 'ic95' 'median'};

    for i=1:length(powerProfile_theta)
        for j=1:sessionInfo.nChannels
            Output_powerProfile_theta{i}{j} = {MergePoints.foldernames{i} powerProfile_theta{i}.channels(j)+1 shank_ch(j) powerProfile_theta{i}.mean(j) powerProfile_theta{i}.std(j)...
                                                powerProfile_theta{i}.ic95(j) powerProfile_theta{i}.median(j)};
        end 
    end

    %% POWER SPECTRUM PROFILE SLOW GAMMA

    sheet_powerProfile_sg = 'powerProfile_sg';
    Header_powerProfile_sg = {'folder' 'ch' 'shankID' 'mean' 'std' 'ic95' 'median'};

    for i=1:length(powerProfile_sg)
        for j=1:sessionInfo.nChannels
            Output_powerProfile_sg{i}{j} = {MergePoints.foldernames{i} powerProfile_sg{i}.channels(j)+1 shank_ch(j) powerProfile_sg{i}.mean(j) powerProfile_sg{i}.std(j)...
                                                powerProfile_sg{i}.ic95(j) powerProfile_sg{i}.median(j)};
        end 
    end

    %% POWER SPECTRUM PROFILE HIGH GAMMA
    sheet_powerProfile_hg = 'powerProfile_hg';
    Header_powerProfile_hg = {'folder' 'ch' 'shankID' 'mean' 'std' 'ic95' 'median'};

    for i=1:length(powerProfile_hg)
        for j=1:sessionInfo.nChannels
            Output_powerProfile_hg{i}{j} = {MergePoints.foldernames{i} powerProfile_hg{i}.channels(j)+1 shank_ch(j) powerProfile_hg{i}.mean(j) powerProfile_hg{i}.std(j)...
                                                powerProfile_hg{i}.ic95(j) powerProfile_hg{i}.median(j)};
        end 
    end

end
%% ====================================================================================
%% -------- WRITING EXCEL FILE ----------------------------------------
% =====================================================================================
%% WRITING EXCEL FILE
try

    NoFile = isfile([pathExcel,'\',nameExcel]);

    cd(pathExcel)
    if ~NoFile
        % Sheet Spikes
        if any(ismember(listOfAnalysis,'spikes'))
            xlswrite([nameExcel], Header_spikes, [sheet_spikes], ['A',num2str(1),':','E',num2str(1)])
        end
        % Sheet Ripples
        if any(ismember(listOfAnalysis,'ripples'))
            xlswrite([nameExcel], Header_ripples, [sheet_ripples],['A',num2str(1),':','H',num2str(1)])
        end
        % Sheet FiringMaps
        if any(ismember(listOfAnalysis,'placeCells'))
            xlswrite([nameExcel], Header_firingMaps, [sheet_firingMaps], ['A',num2str(1),':','AF',num2str(1)])
        end
        % Sheet spikeTrain
        if any(ismember(listOfAnalysis,'spikeTrain'))
            xlswrite([nameExcel], Header_spikeTrain,[sheet_spikeTrain],['A',num2str(1),':','AF',num2str(1)])
        end
        % Sheet Performance
        if any(ismember(listOfAnalysis,'performance'))
            xlswrite([nameExcel],Header_performance,[sheet_performance],['A',num2str(1),':','S',num2str(1)])
        end
        % Sheet PhaseLockingTheta
        if any(ismember(listOfAnalysis,'thetaModulation'))
            xlswrite([nameExcel], Header_phaseLockingTheta, [sheet_phaseLockingTheta], ['A',num2str(1),':','H',num2str(1)])
            %Sheet PhaseLockingSG
            xlswrite([nameExcel], Header_phaseLockingSG, [sheet_phaseLockingSG], ['A',num2str(1),':','H',num2str(1)])
        end
        % Sheet powerProfile_theta
        if any(ismember(listOfAnalysis,'powerSpectrumProfile'))
            xlswrite([nameExcel], Header_powerProfile_theta, [sheet_powerProfile_theta], ['A',num2str(1),':','G',num2str(1)])
            % Sheet powerProfile_sg
            xlswrite([nameExcel], Header_powerProfile_sg, [sheet_powerProfile_sg], ['A',num2str(1),':','G',num2str(1)])
            % Sheet powerProfile_hfo
            xlswrite([nameExcel], Header_powerProfile_hg, [sheet_powerProfile_hg], ['A',num2str(1),':','G',num2str(1)]) 
        end
        % Sheet lfp_analysis
        if any(ismember(listOfAnalysis,'lfp_analysis'))
            if ~isempty(Output_CFC)
                xlswrite([nameExcel],Header_CFCPhaseAmp,[sheet_CFCPhaseAmp],['A',num2str(1),':','H',num2str(1)])
            end
            if ~isempty(coherence_Shanks)
                xlswrite([nameExcel],Header_coherenceShanks,[sheet_coherenceShanks],['A',num2str(1),':','AG',num2str(1)])
            end
            if ~isempty(Output_coherogram)
                xlswrite([nameExcel],Header_coherogram,[sheet_coherogram],['A',num2str(1),':','AE',num2str(1)])
            end
    %         if ~isempty(Output_phaseAmpCoupling)
    %             xlswrite([nameExcel],Header_phaseAmpCouplingByAmp,[sheet_phaseAmpCouplingByAmp],['A',num2str(1),':','C',num2str(1)])
    %         end
            if ~isempty(Output_GMI)
                xlswrite([nameExcel],Header_GMI,[sheet_GMI],['A',num2str(1),':','G',num2str(1)])
            end
            if ~isempty(Output_MI)
                xlswrite([nameExcel],Header_MI,[sheet_MI],['A',num2str(1),':','F',num2str(1)])
            end
        end
    end

    % SPIKES
    if any(ismember(listOfAnalysis,'spikes'))
        for i=1:length(Output_spikes)
            file_excel = openExcel(nameExcel, sheet_spikes);
            length_fileExcel = numel(fieldnames((file_excel)));
            if length_fileExcel ==1
                lineToWrite = 2;
            else
                line = size(file_excel.data(:,:),1);
                lineToWrite = line+2;
            end
            xlswrite([nameExcel], Output_spikes{i}, [sheet_spikes], ['A',num2str(lineToWrite),':','E',num2str(lineToWrite)]);
        end
    end

    % RIPPLE
    if any(ismember(listOfAnalysis,'ripples'))
        for i=1:length(Output_ripples)
            file_excel = openExcel(nameExcel,sheet_ripples);
            length_fileExcel = numel(fieldnames((file_excel))); % length == 1 means that there is only Header in xls file

            if length_fileExcel == 1
                lineToWrite = 2;
            else
                line = size(file_excel.data(:,:),1);
                lineToWrite = line+2;
            end

            xlswrite([nameExcel], Output_ripples{i},[sheet_ripples], ['A',num2str(lineToWrite),':','H',num2str(lineToWrite)]);
        end
    end

    % FIRINGMAPS
    if any(ismember(listOfAnalysis,'placeCells'))
        for i=1:length(Output_firingMaps)
            for j=1:spikes.numcells
                file_excel = openExcel(nameExcel, sheet_firingMaps);
                length_fileExcel = numel(fieldnames((file_excel)));
                if length_fileExcel == 1
                    lineToWrite = 2;
                else
                    line = size(file_excel.data(:,:),1);
                    lineToWrite = line+2;
                end
                xlswrite([nameExcel], Output_firingMaps{i}{j}, [sheet_firingMaps], ['A', num2str(lineToWrite), ':', 'AF', num2str(lineToWrite)])
            end

        end
    end

    % SPIKETRAIN
    if any(ismember(listOfAnalysis,'spikeTrain'))
        for i=1:length(Output_spikeTrain)
            for j=1:spikes.numcells

                file_excel = openExcel(nameExcel,sheet_spikeTrain);
                length_fileExcel = numel(fieldnames((file_excel))); % length == 1 means that there is only Header in excel file

                if length_fileExcel ==1
                    lineToWrite = 2;
                else
                    line = size(file_excel.data(:,:),1);
                    lineToWrite = line+2;
                end

                xlswrite([nameExcel], Output_spikeTrain{i}{j}, [sheet_spikeTrain], ['A',num2str(lineToWrite), ':', 'AF', num2str(lineToWrite)]);
            end
        end
    end

    % PERFORMANCE
    if any(ismember(listOfAnalysis,'performance'))
        for i=1:length(Output_performance)
            file_excel = openExcel(nameExcel,sheet_performance);
            length_fileExcel = numel(fieldnames((file_excel)));

            if length_fileExcel == 1
                lineToWrite = 2;
            else
                line = size(file_excel.data(:,:),1);
                lineToWrite = line+2;
            end

            xlswrite([nameExcel],Output_performance{i},[sheet_performance],['A',num2str(lineToWrite),':','S',num2str(lineToWrite)]);
        end
    end

    % PHASELOCKING THETA
    if any(ismember(listOfAnalysis,'thetaModulation'))
        for i=1:length(Output_phaseLockingTheta)    
            file_excel = openExcel(nameExcel,sheet_phaseLockingTheta);
            length_fileExcel = numel(fieldnames((file_excel)));
            if length_fileExcel == 1
                lineToWrite = 2;
            else
                line = size(file_excel.data(:,:),1);
                lineToWrite = line+2;
            end

            xlswrite([nameExcel], Output_phaseLockingTheta{i},[sheet_phaseLockingTheta], ['A',num2str(lineToWrite),':','H', num2str(lineToWrite)])   
        end

        % PHASELOCKING SG
        for i=1:length(Output_phaseLockingSG)    
            file_excel = openExcel(nameExcel,sheet_phaseLockingSG);
            length_fileExcel = numel(fieldnames((file_excel)));
            if length_fileExcel == 1
                lineToWrite = 2;
            else
                line = size(file_excel.data(:,:),1);
                lineToWrite = line+2;
            end

            xlswrite([nameExcel], Output_phaseLockingSG{i},[sheet_phaseLockingSG], ['A',num2str(lineToWrite),':','H', num2str(lineToWrite)])   
        end
    end

    % POWER PROFILE THETA
    if any(ismember(listOfAnalysis,'powerSpectrumProfile'))
        for i=1:length(Output_powerProfile_theta)
            for j=1:sessionInfo.nChannels
                file_excel = openExcel(nameExcel,sheet_powerProfile_theta);
                length_fileExcel = numel(fieldnames((file_excel)));
                if length_fileExcel == 1
                    lineToWrite = 2;
                else
                    line = size(file_excel.data(:,:),1);
                    lineToWrite = line+2;
                end    
                xlswrite([nameExcel], Output_powerProfile_theta{i}{j}, [sheet_powerProfile_theta], ['A',num2str(lineToWrite),':','G',num2str(lineToWrite)])
            end
        end

        % POWER PROFILE SLOW GAMMA
        for i=1:length(Output_powerProfile_sg)

            for j=1:sessionInfo.nChannels
                file_excel = openExcel(nameExcel,sheet_powerProfile_sg);
                length_fileExcel = numel(fieldnames((file_excel)));
                if length_fileExcel == 1
                    lineToWrite = 2;
                else
                    line = size(file_excel.data(:,:),1);
                    lineToWrite = line+2;
                end    
                xlswrite([nameExcel], Output_powerProfile_sg{i}{j}, [sheet_powerProfile_sg], ['A',num2str(lineToWrite),':','G',num2str(lineToWrite)])
            end
        end

        % POWER PROFILE HIGH FREQUENCY OSCILLATIONS
        for i=1:length(Output_powerProfile_hg)
            for j=1:sessionInfo.nChannels
                file_excel = openExcel(nameExcel,sheet_powerProfile_hg);
                length_fileExcel = numel(fieldnames((file_excel)));
                if length_fileExcel == 1
                    lineToWrite = 2;
                else
                    line = size(file_excel.data(:,:),1);
                    lineToWrite = line+2;
                end    
                xlswrite([nameExcel], Output_powerProfile_hg{i}{j}, [sheet_powerProfile_hg], ['A',num2str(lineToWrite),':','G',num2str(lineToWrite)])
            end
        end
    end

    % LFP ANALYSIS
    if any(ismember(listOfAnalysis,'lfp_analysis'))
        % CFC PHASE AMP
        for j = 1:length(MergePoints.foldernames)
            for i=1:length(Output_CFC)
                file_excel = openExcel(nameExcel,sheet_CFCPhaseAmp);
                length_fileExcel = numel(fieldnames((file_excel)));
                if length_fileExcel == 1
                    lineToWrite = 2;
                else
                    line = size(file_excel.data(:,:),1);
                    lineToWrite = line+2;
                end

                xlswrite([nameExcel],Output_CFC{i}{j},[sheet_CFCPhaseAmp],['A',num2str(lineToWrite),':','H',num2str(lineToWrite)])
            end
        end

        % COHERENCE SHANKS
        if ~isempty(coherence_Shanks)
            for i = 1:length(Output_coherenceShanks)
                for j=1:length(Output_coherenceShanks{i})

                    file_excel = openExcel(nameExcel,sheet_coherenceShanks);
                    length_fileExcel = numel(fieldnames((file_excel)));
                    if length_fileExcel == 1
                        lineToWrite = 2;
                    else
                        line = size(file_excel.data(:,:),1);
                        lineToWrite = line+2;
                    end

                    xlswrite([nameExcel],Output_coherenceShanks{j}{j},[sheet_coherenceShanks],['A',num2str(lineToWrite),':','AG',num2str(lineToWrite)])
                end
            end
        end

        % COHEROGRAM
        for i=1:length(Output_coherogram)
            file_excel = openExcel(nameExcel,sheet_coherogram);
            length_fileExcel = numel(fieldnames((file_excel)));
            if length_fileExcel == 1
                lineToWrite = 2;
            else
                line = size(file_excel.data(:,:),1);
                lineToWrite = line+2;
            end

            xlswrite([nameExcel],Output_coherogram{i},[sheet_coherogram],['A',num2str(lineToWrite),':','AE',num2str(lineToWrite)])
        end

        % GMI   
        for j=1:length(MergePoints.foldernames)
            for i=1:length(Output_GMI)
                file_excel = openExcel(nameExcel,sheet_GMI);
                length_fileExcel = numel(fieldnames((file_excel)));
                if length_fileExcel == 1
                    lineToWrite = 2;
                else
                    line = size(file_excel.data(:,:),1);
                    lineToWrite = line+2;
                end

                xlswrite([nameExcel],Output_GMI{i}{j},[sheet_GMI],['A',num2str(lineToWrite),':','G',num2str(lineToWrite)])
            end
        end

        % MI
        for j = 1:length(MergePoints.foldernames)
            for i=1:length(Output_MI)
                file_excel = openExcel(nameExcel,sheet_MI);
                length_fileExcel = numel(fieldnames((file_excel)));
                if length_fileExcel == 1
                    lineToWrite = 2;
                else
                    line = size(file_excel.data(:,:),1);
                    lineToWrite = line+2;
                end

                xlswrite([nameExcel],Output_MI{i}{j},[sheet_MI],['A',num2str(lineToWrite),':','F',num2str(lineToWrite)])
            end
        end
    end

    cd(basepath)
    disp('Excel has been written successfully...')

catch
    warning('Problems creating Excel...')
end

%% save Excel
% [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
% if saveMat
%     save([basepath filesep sessionInfo.FileName '.Excel.mat'],'excel');
% end

end


