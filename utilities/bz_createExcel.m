function [excel] = bz_createExcel(varargin)
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
addParameter(p,'analyzeSubSessions',false, @islogical);
addParameter(p,'diffLFPs',true, @islogical); % To compute phase-locking with the channel where the amplitude of the unit is maximum and not with a reference channel for all units
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'getWaveformsFromDat',true,@islogical);
addParameter(p,'showWaveforms',true,@islogical);
addParameter(p,'firingMaps',[],@isstruct);
addParameter(p,'ripples',[],@isstruct);
addParameter(p,'spikeTrain',[],@isstruct);
addParameter(p,'behaviour',[],@isstruct);
addParameter(p,'performance',[], @isstruct);
addParameter(p,'thetaGamma',[], @isstruct);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'pathExcel','F:\data',@isstr);
addParameter(p,'nameExcel','Excel Prueba',@isstr);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'PhaseLockingData',[],@isstruct);
addParameter(p,'PhaseLockingData_sg',[],@issstruct);
addParameter(p,'powerProfile_theta',[],@isstruct);
addParameter(p,'powerProfile_sg',[],@isstruct);
addParameter(p,'powerProfile_hfo',[],@isstruct);

parse(p,varargin{:});
basepath = p.Results.basepath;
analyzeSubSessions = p.Results.analyzeSubSessions;
diffLFPs = p.Results.diffLFPs;
spikes = p.Results.spikes;
getWaveformsFromDat = p.Results.getWaveformsFromDat;
showWaveforms = p.Results.showWaveforms;
firingMaps = p.Results.firingMaps;
ripples = p.Results.ripples;
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
powerProfile_hfo = p.Results.powerProfile_hfo;

%% In case Excel already exists 
if ~isempty(dir([basepath filesep '*Excel.mat'])) || forceReload
    disp('Excel already detected! Loading file.');
    file = dir([basepath filesep '*Excel.Behavior.mat']);
    load(file.name);
    return
end

%% LOAD TRACKING

if isempty(tracking)
    tracking = getSessionTracking();
end

%% LOAD BEHAVIOUR
if isempty(behaviour)
    behaviour = getSessionBehaviour_v2();
end
%% LOAD MERGEPOINTS
if ~isempty(dir([basepath filesep '*MergePoints.events.mat']))
    disp('Loading MergePoints...')
    file = dir([basepath filesep '*MergePoints.events.mat']);
    load(file.name)
end
    
%% LOAD SPIKES
if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'showWaveforms',showWaveforms);
end

%% LOAD FIRINGMAPS
if isempty(firingMaps)
    if ~isempty(dir([basepath filesep '*firingMapsAvg.cellinfo.mat']))
        disp('Firing Maps already detected! Loading file.');
        file = dir([basepath filesep '*firingMapsAvg.cellinfo.mat']);
        load(file.name);
    end
end

%% LOAD PERFORMANCE

if isempty(performance)
    try
        performance = getSessionPerformance('tracking',tracking);
        
    catch
        
        warning('It is not possible to load Performance...');
        performance = [];
        
    end
end

%% PHASELOCKING THETA-GAMMA
if diffLFPs
    if ~isempty(dir([basepath filesep '*PhaseLockingData.diffLFPs.cellinfo.mat']))
        disp('Loading SubSession and different LFPs Phase Locking Theta...')
        file = dir([basepath filesep '*PhaseLockingData.diffLFPs.cellinfo.mat'])
        load(file.name)
    end

    if ~isempty(dir([basepath filesep '*PhaseLockingData_sg.diffLFPs.cellinfo.mat']))
        disp('Loading SubSession and different LFPs Phase Locking Slow Gamma...')
        file = dir([basepath filesep '*PhaseLockingData_sg.diffLFPs.cellinfo.mat'])
        load(file.name) 
    end
else
            
end
        
%% LOAD RIPPLES AND SPIKETRAIN AND POWERPROFILE
if analyzeSubSessions
    % RIPPLES
    if isempty(ripples)
        if ~isempty(dir([basepath filesep '*ripples.SubSession.events.mat']))
            disp('Loading SubSession Ripples...')
            file = dir([basepath filesep '*ripples.SubSession.events.mat']);
            load(file.name)
        else
            warning('It is not possible to load ripple SubSessions ...')
        end
    end 
    
    % SPIKETRAIN
    if isempty(spikeTrain)
        spikeTrain = bz_SpikeTrain(spikes,'analyzeSubSessions',analyzeSubSessions);
    end 
    
    % POWER PROFILE
    if isempty(powerProfile_theta)
        
        if ~isempty(dir([basepath filesep '*PowerSpectrumProfile_6_12.SubSession.channelinfo.mat']))
            disp('Loading SubSession Power Profile Theta...')
            file = dir([basepath filesep '*PowerSpectrumProfile_6_12.SubSession.channelInfo.mat'])
            for i=1:length(file)
                powerProfile_theta{i} = load(file(i).name);
            end
            
        end
        
    end
    
    
    if isempty(powerProfile_sg)
        if ~isempty(dir([basepath filesep '*PowerSpectrumProfile_30_60.SubSession.channelinfo.mat']))
            disp('Loading SubSession Power Profile Slow Gamma...')
            file = dir([basepath filesep '*PowerSpectrumProfile_30_60.SubSession.channelinfo.mat'])
            for i=1:length(file)
                powerProfile_sg{i} = load(file(i).name);
            end
            
        end
    end
    
    if isempty(powerProfile_hfo)
        if ~isempty(dir([basepath filesep '*PowerSpectrumProfile_120_250.SubSession.channelinfo.mat']))
            disp('Loading SubSession Power Profile High Frequency Oscillations ...')
            file = dir([basepath filesep '*PowerSpectrumProfile_120_250.SubSession.channelinfo.mat'])
            for i=1:length(file)
                powerProfile_hfo{i} = load(file(i).name);
            end
            
        end
    end
    
    
    
end
%% RIPPLES

sheet_ripples = 'ripples';
Header_ripples = {'foldername' 'numRipples' 'peakFrequency' 'peakAmplitude' 'duration' };

for i=1:length(ripples)
    Output_ripples{i} = {ripples{i}.foldername size(ripples{i}.timestamps,1) mean(ripples{i}.data.peakFrequency) mean(ripples{i}.data.peakAmplitude) mean(ripples{i}.data.duration)};
end



%% SPIKES
sheet_spikes = 'spikes';
Header_spikes = {'foldername' 'shankID' 'ch' 'totalSpikes' 'amplitude' }


for j=1:spikes.numcells
    Output_spikes{j} = {spikes.basename spikes.shankID(j) spikes.maxWaveformCh(j) spikes.total(j) mean(spikes.amplitudes{j})};
end
%% FIRING MAPS

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

%% SPIKE TRAIN
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


%% PERFORMANCE


%% PHASELOCKINGDATA (THETA)
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

%% POWER SPECTRUM PROFILE THETA
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
        Output_powerProfile_theta{i}{j} = {MergePoints.foldernames{i} powerProfile_theta{i}.powerProfile.channels(j)+1 shank_ch(j) powerProfile_theta{i}.powerProfile.mean(j) powerProfile_theta{i}.powerProfile.std(j)...
                                            powerProfile_theta{i}.powerProfile.ic95(j) powerProfile_theta{i}.powerProfile.median(j)};
    end 
end

%% POWER SPECTRUM PROFILE SLOW GAMMA

sheet_powerProfile_sg = 'powerProfile_sg';
Header_powerProfile_sg = {'folder' 'ch' 'shankID' 'mean' 'std' 'ic95' 'median'};

for i=1:length(powerProfile_sg)
    for j=1:sessionInfo.nChannels
        Output_powerProfile_sg{i}{j} = {MergePoints.foldernames{i} powerProfile_sg{i}.powerProfile.channels(j)+1 shank_ch(j) powerProfile_sg{i}.powerProfile.mean(j) powerProfile_sg{i}.powerProfile.std(j)...
                                            powerProfile_sg{i}.powerProfile.ic95(j) powerProfile_sg{i}.powerProfile.median(j)};
    end 
end

%% POWER SPECTRUM PROFILE HIGH FREQUENCY OSCILLATIONS
sheet_powerProfile_hfo = 'powerProfile_hfo';
Header_powerProfile_hfo = {'folder' 'ch' 'shankID' 'mean' 'std' 'ic95' 'median'};

for i=1:length(powerProfile_hfo)
    for j=1:sessionInfo.nChannels
        Output_powerProfile_hfo{i}{j} = {MergePoints.foldernames{i} powerProfile_hfo{i}.powerProfile.channels(j)+1 shank_ch(j) powerProfile_hfo{i}.powerProfile.mean(j) powerProfile_hfo{i}.powerProfile.std(j)...
                                            powerProfile_hfo{i}.powerProfile.ic95(j) powerProfile_hfo{i}.powerProfile.median(j)};
    end 
end



%% WRITING EXCEL FILE

NoFile = isfile([pathExcel,'\',nameExcel]);

cd(pathExcel)
if ~NoFile
    % Sheet Spikes
    xlswrite([nameExcel], Header_spikes, [sheet_spikes], ['A',num2str(1),':','E',num2str(1)])
    % Sheet FiringMaps
    xlswrite([nameExcel], Header_firingMaps, [sheet_firingMaps], ['A',num2str(1),':','AF',num2str(1)])
    % Sheet Ripples
    xlswrite([nameExcel], Header_ripples, [sheet_ripples],['A',num2str(1),':','E',num2str(1)])
    % Sheet spikeTrain
    xlswrite([nameExcel], Header_spikeTrain,[sheet_spikeTrain],['A',num2str(1),':','AF',num2str(1)])
    % Sheet Performance
    
    % Sheet PhaseLockingTheta
    xlswrite([nameExcel], Header_phaseLockingTheta, [sheet_phaseLockingTheta], ['A',num2str(1),':','H',num2str(1)])
    %Sheet PhaseLockingSG
    xlswrite([nameExcel], Header_phaseLockingSG, [sheet_phaseLockingSG], ['A',num2str(1),':','H',num2str(1)])
    % Sheet powerProfile_theta
    xlswrite([nameExcel], Header_powerProfile_theta, [sheet_powerProfile_theta], ['A',num2str(1),':','G',num2str(1)])
    % Sheet powerProfile_sg
    xlswrite([nameExcel], Header_powerProfile_sg, [sheet_powerProfile_sg], ['A',num2str(1),':','G',num2str(1)])
    % Sheet powerProfile_hfo
    xlswrite([nameExcel], Header_powerProfile_hfo, [sheet_powerProfile_hfo], ['A',num2str(1),':','G',num2str(1)])
    
    
end

% SPIKES
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

% FIRINGMAPS
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

% RIPPLE
for i=1:length(Output_ripples)
    file_excel = openExcel(nameExcel,sheet_ripples);
    length_fileExcel = numel(fieldnames((file_excel))); % length == 1 means that there is only Header in xls file

    if length_fileExcel == 1
        lineToWrite = 2;
    else
        line = size(file_excel.data(:,:),1);
        lineToWrite = line+2;
    end

    xlswrite([nameExcel], Output_ripples{i},[sheet_ripples], ['A',num2str(lineToWrite),':','E',num2str(lineToWrite)]);
end

% SPIKETRAIN
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

% PHASELOCKING THETA

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

% POWER PROFILE THETA

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

for i=1:length(Output_powerProfile_hfo)
    
    for j=1:sessionInfo.nChannels
        file_excel = openExcel(nameExcel,sheet_powerProfile_hfo);
        length_fileExcel = numel(fieldnames((file_excel)));
        if length_fileExcel == 1
            lineToWrite = 2;
        else
            line = size(file_excel.data(:,:),1);
            lineToWrite = line+2;
        end    
        xlswrite([nameExcel], Output_powerProfile_hfo{i}{j}, [sheet_powerProfile_hfo], ['A',num2str(lineToWrite),':','G',num2str(lineToWrite)])
    end
end


cd(basepath)

%% save Excel
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
if saveMat
    save([basepath filesep sessionInfo.FileName '.Excel.mat'],'excel');
end

end
