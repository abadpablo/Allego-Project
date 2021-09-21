function [outputArg1,outputArg2] = bz_analyzeSession(varargin)
% Runs a list of descriptive analysis
%
% USAGE:
%   bz_analyzeSession(varargin)
%
% INPUT:
%   basepath              By default pwd
%   listOfAnalysis        Summary analysis that will be computed(as a cell
%                           with strings). Default 'all'. Possible values:
%                           'spikes','analogPulses','digitalPulses','ripples',
%                           'tMazeBehaviour','linearMazeBehaviour','OpenFieldBehaviour',
%                           'YMazeBehaviour','thetaModulation'
% exclude                 Cell with strings with the list of analysis to
%                         exluce
%
% (specific analysis options)
%  exlcludeShanks          Default []
%  getWaveformsFromDat     From 'spikes' summary, default false
%  analogChannelsList      Array of Channel to perform 'analogPulses' psth
%                           and csd. Default 'all'
%  digitalChannelList      Array of channel to perform 'digitalPulses' psth
%                           and csd. Default 'all'.
% analyzeSubSessions       Default false. Indicate if the analysis will be
%                           performed in individual subsessions instead of the whole session
%
% Pablo Abad 2021. Adapted from Manu Valero (computeSessionSummary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'listOfAnalysis','all',@iscellstr);
addParameter(p,'exclude',[],@iscellstr);
addParameter(p,'excludeShanks',[],@isnmueric);
addParameter(p,'getWaveformsFromDat',false,@islogical);
addParameter(p,'analogChannelsList','all',@isnumeric);
addParameter(p,'digitalChannelsList','all',@isnumeric);
addParameter(p,'analyzeSubSessions',false,@islogical);
addParameter(p,'showWaveforms',true,@islogical);
addParameter(p,'forceReloadRipples',false,@islogical);
addParameter(p,'diffLFPs',false,@islogical);
addParameter(p,'cellClassification',false,@islogical);
addParameter(p,'thetaFreq',[4 12],@isnumeric);
addParameter(p,'sgFreq',[30 65], @isnumeric);
addParameter(p,'hgFreq',[66 130],@isnumeric);
addParameter(p,'hfoFreq',[150 185], @isnumeric);
addParameter(p,'showFig',false,@isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
listOfAnalysis = p.Results.listOfAnalysis;
exclude = p.Results.exclude;
excludeShanks = p.Results.excludeShanks;
getWaveformsFromDat = p.Results.getWaveformsFromDat;
analogChannelsList = p.Results.analogChannelsList;
digitalChannelsList = p.Results.digitalChannelsList;
analyzeSubSessions = p.Results.analyzeSubSessions;
showWaveforms = p.Results.showWaveforms;
diffLFPs = p.Results.diffLFPs;
thetaFreq = p.Results.thetaFreq;
sgFreq = p.Results.sgFreq;
hgFreq = p.Results.hgFreq;
hfoFreq = p.Results.hfoFreq;
showFig = p.Results.showFig;

prevPath = pwd;
cd(basepath);

if ischar(listOfAnalysis) && strcmpi(listOfAnalysis,'all')
    listOfAnalysis = {'spikes','digitalPulses','ripples','tMazeBehaviour','linearMazeBehaviour','YMazeBehaviour','OpenFieldBehaviour','thetaModulation','behaviour','spikeTrain','performance','excel','CellExplorer','lfp_analysis'};
end

if ~isempty(exclude)
    listOfAnalysis(ismember(listOfAnalysis,exclude)) = [];
end

session = loadSession;
excludeChannels = [];

for ii=1:length(excludeShanks)
    excludeChannels = session.extracellular.electrodeGroups.channels(excludeShanks(ii));
end

mkdir('SummaryFigures'); % create folder
mkdir('PhaseModulationFig')
close all

%% 1 - Spikes Summary

if any(ismember(listOfAnalysis,'spikes'))
    try
          
        disp('Spike-waveform, ACG and cluster location...');
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);

        % plot spikes summary
        disp('Plotting spikes summary...');
        if exist('chanMap.mat','file')
            load('chanMap.mat','xcoords','ycoords','kcoords');
            xcoords = xcoords - min(xcoords); xcoords = xcoords/max(xcoords);
            ycoords = ycoords - min(ycoords); ycoords = ycoords/max(ycoords); 
            xcoords(excludeChannels) = [];
            ycoords(excludeChannels) = [];
            kcoords(excludeChannels) = [];
            
        else
            xcoords = NaN;
            ycoords = NaN;
            kcoords = NaN;
        end
        ccg=CCG(spikes.times,[],'binSize',0.001,'duration',0.08);
        dur = 0.08;
        xt = linspace(-dur/2*1000,dur/2*1000,size(ccg,1));
        figure
        set(gcf,'Position',[100 -100 2500 1200])
        for jj = 1:size(spikes.UID,2)
            fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
            area(xt,ccg(:,jj,jj),'LineStyle','none');
            try 
                xt2 = linspace(dur/3*1000+5,dur/2*1000+5+80,length(spikes.filtWaveform{jj})); % waveform
                spk = spikes.filtWaveform{jj} - min(spikes.filtWaveform{jj}); spk = spk/max(spk) * max(ccg(:,jj,jj));
                hold on
                    plot(xt2,spk,'color','.8 .2 .2');

                plot((xcoords*30) + dur/2*1000+5+60, ycoords*max(ccg(:,jj,jj)),'.','color',[.8 .8 .8],'MarkerSize',5); % plotting xml
                plot((xcoords(spikes.maxWaveformCh1(jj))*30) + dur/2*1000+5+60,...
                    ycoords(spikes.maxWaveformCh1(jj))*max(ccg(:,jj,jj)),'.','color',[.1 .1 .1],'MarkerSize',10); % plotting xml
            end
            title(num2str(jj),'FontWeight','normal','FontSize',10);
            if jj == 1
                ylabel('Counts/ norm');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (100 ms /1.5ms)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,'SummaryFigures\spikes AutoCorr.png');

        win = [-0.3 0.3];
        disp('Plotting CCG...');
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        [allCcg, t] = CCG(spikes.times,[],'binSize',0.005,'duration',0.6);
        indCell = [1:size(allCcg,2)];
        for jj = 1:size(spikes.UID,2)
            fprintf(' **CCG from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
            subplot(7,ceil(size(spikes.UID,2)/7),jj);
            cc = zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2); % get crosscorr
            imagesc(t,1:max(indCell)-1,cc)
            set(gca,'YDir','normal'); colormap jet; caxis([-abs(max(cc(:))) abs(max(cc(:)))])
            hold on
            zmean = mean(zscore(squeeze(allCcg(:,jj,indCell(indCell~=jj)))',[],2));
            zmean = zmean - min(zmean); zmean = zmean/max(zmean) * (max(indCell)-1) * std(zmean);
            plot(t, zmean,'k','LineWidth',2);
            xlim([win(1) win(2)]); ylim([0 max(indCell)-1]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);

            if jj == 1
                ylabel('Cell');
            elseif jj == size(spikes.UID,2)
                xlabel('Time (s)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,'SummaryFigures\spikes CrossCorr.png'); 
    catch
        warning('Error on Spike-waveform, autocorrelogram and cluster location! ');
    end
end


%% 2 - Ripples CSD and PSTH

if any(ismember(listOfAnalysis,'ripples'))
    disp('Ripples CSD and PSTH...');
    [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
%     lfp = bz_GetLFP(sessionInfo.channels,'noPrompts',true);
    spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);
    forceReloadRipples = true;
    if isempty(dir('*ripples.events.mat')) || forceReloadRipples
        if analyzeSubSessions  
            try
                disp('Ripples CSD and PSTH for SubSession...')
            if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
                load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
                count = 1;
                for ii=1:size(MergePoints.foldernames,2)
                    timestamps = MergePoints.timestamps(ii,:);
                    foldername = MergePoints.foldernames{ii};
                    lfp = bz_GetLFP(sessionInfo.channels,'restrict',timestamps);
                    
                    rippleChannels{ii} = computeRippleChannel('discardShanks',excludeShanks,'timestamps_subSession',timestamps,'foldername',foldername,'saveMat',false);
                    if ~isempty(rippleChannels{ii}.Ripple_Channel)
                        ripples{ii} = bz_FindRipples(basepath,rippleChannels{ii}.Ripple_Channel,'timestamps_subSession',timestamps,'foldername',foldername);
                    
                        shanks = sessionInfo.AnatGrps;
                        shanks(excludeShanks) = [];
                        passband = [130 200];
                        filtered = bz_Filter(lfp,'channels',rippleChannels{ii}.Ripple_Channel,'filter','butter','passband',passband,'order',3);
                        [maps,data,stats] = bz_RippleStats(filtered.data,lfp.timestamps,ripples{ii});
                        ripples{ii}.maps = maps;
                        ripples{ii}.data = data;
                        ripples{ii}.stats = stats;
                        ripples{ii}.foldername = foldername;
                    
                        % CSD
                        twin = 0.1;
                        evs = ripples{ii}.peaks;
                        figure
                        set(gcf,'Position',[100 100 1400 600])
                        for jj = 1:size(shanks,2)
                            lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts', true);
                            [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                            taxis = linspace(-twin,twin,size(csd.data,1));
                            cmax = max(max(csd.data)); 
                            subplot(1,size(shanks,2),jj);
                            contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                            set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('RIPPLES, Shank #',num2str(jj)),'FontWeight','normal'); 
                            colormap jet; caxis([-cmax cmax]);
                            hold on
                            for kk = 1:size(lfpAvg.data,2)
                                plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                            end
                        end

                        saveas(gcf,['SummaryFigures\',strcat('ripplesCSD.',foldername),'.png']);  
                    
                        % PSTH
                        st = ripples{ii}.peaks;
                        spikeResponse = [];
                        win = [-0.2 0.2];
                        figure
                        set(gcf,'Position',[100 -100 2500 1200])

                        for jj = 1:size(spikes.UID,2)
                            fprintf(' **Ripple from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                            rast_x = []; rast_y = [];
                            for kk = 1:length(st)
                                temp_rast = spikes.times{jj} - st(kk);
                                temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                                rast_x = [rast_x temp_rast'];
                                rast_y = [rast_y kk*ones(size(temp_rast))'];
                            end
                            [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.005,'duration',1);
                            spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,1:end-1)))'];
                            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                            plot(rast_x, rast_y,'.','MarkerSize',1)
                            hold on
                            plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),2,1) * kk/max(stccg(:,2,1))/2,'k','LineWidth',2);
                            xlim([win(1) win(2)]); ylim([0 kk]);
                            title(num2str(jj),'FontWeight','normal','FontSize',10);

                            if jj == 1
                                ylabel('Trial');
                            elseif jj == size(spikes.UID,2)
                                xlabel('Time (s)');
                            else
                                set(gca,'YTick',[],'XTick',[]);
                            end
                        end
                        saveas(gcf,['SummaryFigures\',strcat('ripplesRaster.',foldername),'.png']); 

                        figure
                        imagesc([t(1) t(end)],[1 size(spikeResponse,2)], spikeResponse); caxis([-3 3]); colormap(jet);
                        xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');
                        saveas(gcf,['SummaryFigures\',strcat('ripplesPsth.',foldername),'.png']); title('Ripples');
                    else
                        ripples{ii}.maps = [];
                        ripples{ii}.data = [];
                        ripples{ii}.stats = [];
                        ripples{ii}.foldername = foldername;
                    end
                end
                save(fullfile(basepath, [sessionInfo.session.name, '.ripples.SubSession.events.mat']),'ripples')
                % Maybe saving also the rippleChannels the same way
                save(fullfile(basepath,[sessionInfo.session.name,'.channelInfo.SubSession.ripples.mat']),'rippleChannels')
            end
            catch
                warning('Error on CSD and PSTH from ripples SubSessions !');
            end            
        else  
            try
                disp('Ripples CSD and PSTH for whole Session ...')
                
                rippleChannels = computeRippleChannel('discardShanks',excludeShanks);
                ripples = bz_FindRipples(basepath,rippleChannels.Ripple_Channel,'saveMat',true,'show','');

                shanks = sessionInfo.AnatGrps;
                shanks(excludeShanks) = [];
                passband = [130 200];
                lfp = bz_GetLFP(sessionInfo.channels,'noPrompts',true);
                filtered = bz_Filter(lfp,'channels',rippleChannels.Ripple_Channel,'filter','butter','passband',passband,'order',3);
                [maps,data,stats] = bz_RippleStats(filtered.data,lfp.timestamps,ripples);
                ripples.maps = maps;
                ripples.data = data;
                ripples.stats = stats;
                
                % CSD
                twin = 0.1;
                evs = ripples.peaks;
                figure
                set(gcf,'Position',[100 100 1400 600])
                for jj = 1:size(shanks,2)
                    lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts', true);
                    [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    cmax = max(max(csd.data)); 
                    subplot(1,size(shanks,2),jj);
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('RIPPLES, Shank #',num2str(jj)),'FontWeight','normal'); 
                    colormap jet; caxis([-cmax cmax]);
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                    end
                end
                    
                saveas(gcf,['SummaryFigures\ripplesCSD.png']);  
                    
                    % PSTH
                    st = ripples.peaks;
                    spikeResponse = [];
                    win = [-0.2 0.2];
                    figure
                    set(gcf,'Position',[100 -100 2500 1200])
                    
                    for jj = 1:size(spikes.UID,2)
                        fprintf(' **Ripple from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                        rast_x = []; rast_y = [];
                        for kk = 1:length(st)
                            temp_rast = spikes.times{jj} - st(kk);
                            temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                            rast_x = [rast_x temp_rast'];
                            rast_y = [rast_y kk*ones(size(temp_rast))'];
                        end
                        [stccg, t] = CCG({spikes.times{jj} st},[],'binSize',0.005,'duration',1);
                        spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,1:end-1)))'];
                        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                        plot(rast_x, rast_y,'.','MarkerSize',1)
                        hold on
                        plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),2,1) * kk/max(stccg(:,2,1))/2,'k','LineWidth',2);
                        xlim([win(1) win(2)]); ylim([0 kk]);
                        title(num2str(jj),'FontWeight','normal','FontSize',10);

                        if jj == 1
                            ylabel('Trial');
                        elseif jj == size(spikes.UID,2)
                            xlabel('Time (s)');
                        else
                            set(gca,'YTick',[],'XTick',[]);
                        end
                    end
                    saveas(gcf,['SummaryFigures\ripplesRaster.png']); 
                    
                    figure
                    imagesc([t(1) t(end)],[1 size(spikeResponse,2)], spikeResponse); caxis([-3 3]); colormap(jet);
                    xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');
                    saveas(gcf,['SummaryFigures\ripplesPsth.png']); title('Ripples');
                    save(fullfile(basepath, [sessionInfo.session.name '.ripples.events.mat']),'ripples')
            catch
                warning('Error on CSD and PSTH from ripples Whole Session !');
            end   
        end
    end
end


%% 3 - THETA AND GAMMA PHASE MODULATION BY SPIKES and POWER SPECTRUM PROFILE

if any(ismember(listOfAnalysis,'thetaModulation'))
    
    if analyzeSubSessions
        try
            disp('Theta, gamma and HFO modulation for SubSessions ...')
            [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
            spikes = loadSpikes('getWaveformsFromDat',false);
            if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
                load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
                count = 1;
                for i = 1:size(MergePoints.foldernames,2)
                    timestamps = MergePoints.timestamps(i,:);
                    foldername = MergePoints.foldernames{i};
                    disp(['Theta, gamma and HFO modultaion for folder: ', foldername])
                    % Theta Profile
                    powerProfile_theta{i} = bz_PowerSpectrumProfile([thetaFreq(1) thetaFreq(2)], 'useparfor',false,'timestampsSubSession',timestamps,'foldername',foldername,'channels', sessionInfo.channels, 'showfig',true,'saveMat',false)
                    % Low Gamma Profile                  
                    powerProfile_sg{i} = bz_PowerSpectrumProfile([sgFreq(1) sgFreq(2)], 'useparfor',false,'timestampsSubSession',timestamps,'foldername',foldername,'channels',sessionInfo.channels,'showfig',true,'saveMat',false);
                    % HFO Profile
                    powerProfile_hg{i} = bz_PowerSpectrumProfile([hgFreq(1) hgFreq(2)], 'useparfor',false,'timestampsSubSession',timestamps,'foldername',foldername,'channels',sessionInfo.channels,'showfig',true,'saveMat',false);
                    % get channels of interest % max theta power above pyr layer
                    [~, a] = max(powerProfile_theta{i}.mean);
                    region{i}.CA1sp = powerProfile_theta{i}.channels(a);

                    for ii = 1:size(sessionInfo.AnatGrps,2)
                        if ismember(region{i}.CA1sp, sessionInfo.AnatGrps(ii).Channels)
                            p_th{i} = powerProfile_theta{i}.mean(sessionInfo.AnatGrps(ii).Channels + 1);                
                            p_sg{i} = powerProfile_sg{i}.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                            p_hg{i} = powerProfile_hg{i}.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                            channels = sessionInfo.AnatGrps(ii).Channels;
                            [~,a] = max(p_th{i}(1:find(channels == region{i}.CA1sp)));
                            region{i}.CA1so = channels(a);
                            [~,a] = max(p_th{i}(find(channels == region{i}.CA1sp):end));
                            region{i}.CA1slm = channels(a - 1+ find(channels == region{i}.CA1sp));

                            figure,
                            hold on
                            plot(1:length(channels), zscore(p_th{i}));
                            plot(1:length(channels), zscore(p_hg{i}));
                            plot(1:length(channels), zscore(p_sg{i}));
                            xlim([1 length(channels)]);
                            set(gca,'XTick',1:length(channels),'XTickLabel',channels,'XTickLabelRotation',45);
                            ax = axis;
                            xlabel(strcat('Channels (neuroscope)-Shank',num2str(ii))); ylabel('power (z)');
                            plot(find(channels == region{i}.CA1sp)*ones(2,1),ax([3 4]),'-k');
                            plot(find(channels == region{i}.CA1so)*ones(2,1),ax([3 4]),'--k');
                            plot(find(channels == region{i}.CA1slm)*ones(2,1),ax([3 4]),'-.k');
                            legend('4-12Hz', 'hfo', '30-60','pyr','~or','~slm');

                            saveas(gcf,['SummaryFigures\regionDef_',foldername,'.png']);
                        end
                    end

                    % power profile of pyr channel of all session
                    lfpT = bz_GetLFP(region{i}.CA1sp,'restrict',timestamps,'noPrompts',true);
                    params.Fs = lfpT.samplingRate; params.fpass = [2 200]; params.tapers = [3 5]; params.pad = 1;
                    [S,t,f] = mtspecgramc(single(lfpT.data),[2 1], params);
                    S = log10(S); % in Db
                    S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending        
                    figure;
                    subplot(3,4,[1 2 3 5 6 7])
                    imagesc(t,f,S_det',[-1.5 1.5]);
                    set(gca,'XTick',[]); ylabel('Freqs');
                    subplot(3,4,[4 8]);
                    plot(mean(S,1),f);
                    set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
                    saveas(gcf,['SummaryFigures\spectrogram_',foldername,'.png']);

                    % phase modulation by spikes

                    if diffLFPs
                        filename = dir('*PhaseLockingData.SubSession.diffLFPs.mat*'); 
                        if ~isempty(filename)
                            disp(['Phase Locking Data SubSession Different LFPs already detected ! Loading file: ', filename.name]);
                            load(filename.name)
                            PLD = PhaseLockingData;
                        else
                            for ii = 1:spikes.numcells
                                lfp = bz_GetLFP(spikes.maxWaveformCh1(ii)-1,'noPrompts',true);
                                PLD{i}{ii} = bz_PhaseModulationDifferentLFPs(spikes.times{ii},lfp,[thetaFreq(1) thetaFreq(2)], 'intervals',timestamps,'plotting',false,'method','wavelet');
                            end
                            PhaseLockingData{i} = PLD{i};

                        end
                        if showFig
                            disp('Theta modulation...')
                            set(gcf,'Position',[100 -100 2500 1200]);
                            for jj = 1:size(spikes.UID,2)
                                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                area([PLD{i}{jj}.phasebins; PLD{i}{jj}.phasebins + pi*2],[PLD{i}{jj}.phasedistros; PLD{i}{jj}.phasedistros],'EdgeColor','none');
                                hold on
                                ax = axis;
                                x = 0:.001:4*pi;
                                y = cos(x);
                                y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                                h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                                xlim([0 4*pi]);
                                title(num2str(jj),'FontWeight','normal','FontSize',10);

                                if jj == 1
                                    ylabel('prob');
                                elseif jj == size(spikes.UID,2)
                                    set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                    xlabel('phase (rad)');
                                else
                                    set(gca,'YTick',[],'XTick',[]);
                                end 
                            end
                            saveas(gcf,['PhaseModulationFig\thetaModDiffLFPS_',foldername,'.png']);
                        end

                        % gamma modulation
                        disp('Gamma modulation...');
                        filename = dir('*PhaseLockingData_sg.diffLFPs.SubSession.cellinfo.mat*');
                        if ~isempty(filename)
                            disp(['Phase Locking Data_sg SubSession Different LFPs already detected ! Loading file: ', filename.name]);
                            load(filename.name)
                            PLD_sg = PhaseLockingData_sg;
                        else
                            for ii=1:spikes.numcells
                                lfp = bz_GetLFP(spikes.maxWaveformCh1(ii)-1,'noPrompts',true,'restrict',timestamps);
                                PLD_sg{i}{ii} = bz_PhaseModulationDifferentLFPs(spikes.times{ii},lfp,[sgFreq(1) sgFreq(2)],'intervals',timestamps,'plotting',false,'method','wavelet');                              
                            end
                            PhaseLockingData_sg{i} = PLD_sg{i};
                        end
                        
                        if showFig
                            figure
                            set(gcf,'Position',[100 -100 2500 1200]);
                            for jj = 1:size(spikes.UID,2)
                                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                area([PLD_sg{i}{jj}.phasebins; PLD_sg{i}{jj}.phasebins + pi*2],[PLD_sg{i}{jj}.phasedistros; PLD_sg{i}{jj}.phasedistros],'EdgeColor','none');
                                hold on
                                ax = axis;
                                x = 0:.001:4*pi;
                                y = cos(x);
                                y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                                h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                                xlim([0 4*pi]);
                                title(num2str(jj),'FontWeight','normal','FontSize',10);

                                if jj == 1
                                    ylabel('prob');
                                elseif jj == size(spikes.UID,2)
                                    set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                    xlabel('phase (rad)');
                                else
                                    set(gca,'YTick',[],'XTick',[]);
                                end
                            end
                            saveas(gcf,['PhaseModulationFig\sgModDiffLFPS_',foldername,'.png']);
                        end

                    else
                        PLD{i} = bz_PhaseModulation(spikes,lfpT,[thetaFreq(1) thetaFreq(2)], 'intervals',timestamps,'plotting',true,'method','wavelet','saveMat',true);
                        if showFig
                            disp('Theta modlation...')
                            set(gcf,'Position',[100 -100 2500 1200]);
                            for jj = 1:size(spikes.UID,2)
                                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                area([PLD{i}.phasebins; PLD{i}.phasebins + pi*2],[PLD{i}.phasedistros(:,jj); PLD.phasedistros{i}(:,jj)],'EdgeColor','none');
                                hold on
                                ax = axis;
                                x = 0:.001:4*pi;
                                y = cos(x);
                                y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                                h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                                xlim([0 4*pi]);
                                title(num2str(jj),'FontWeight','normal','FontSize',10);

                                if jj == 1
                                    ylabel('prob');
                                elseif jj == size(spikes.UID,2)
                                    set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                    xlabel('phase (rad)');
                                else
                                    set(gca,'YTick',[],'XTick',[]);
                                end
                            end
                            saveas(gcf,['SummaryFigures\thetaMod_',foldername,'.png']);
                        end

                        % gamma modulation
                        PLD_sg{i} = bz_PhaseModulation(spikes,lfpT,[sgFreq(1) sgFreq(2)],'intervals',timestamps,'plotting',false,'method','wavelet','saveMat',true);     
                        disp('Gamma modulation...');
                        if showFig
                            figure
                            set(gcf,'Position',[100 -100 2500 1200]);
                            for jj = 1:size(spikes.UID,2)
                                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                                area([PLD_sg{i}.phasebins; PLD_sg{i}.phasebins + pi*2],[PLD_sg{i}.phasedistros(:,jj); PLD_sg{i}.phasedistros(:,jj)],'EdgeColor','none');
                                hold on
                                ax = axis;
                                x = 0:.001:4*pi;
                                y = cos(x);
                                y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                                h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                                xlim([0 4*pi]);
                                title(num2str(jj),'FontWeight','normal','FontSize',10);

                                if jj == 1
                                    ylabel('prob');
                                elseif jj == size(spikes.UID,2)
                                    set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                                    xlabel('phase (rad)');
                                else
                                    set(gca,'YTick',[],'XTick',[]);
                                end
                            end
                            saveas([gcf,'SummaryFigures\gammaMod_',foldername,'.png']);
                        end
                    end    
                end
                % save Power Spectrum Profile theta
                save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_theta{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_theta{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_theta');                
                % save Power Spectrum Profile sg
                save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_sg{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_sg{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_sg');
                % save Power Spectrum Profile hg
                save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_hg{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_hg{1}.processinginfo.params.frange(2)),'.SubSession.channelinfo.mat'],'powerProfile_hg');
                save([sessionInfo.FileName,'.region.mat'],'region');
                save([sessionInfo.session.name '.PhaseLockingData.SubSession.cellinfo.mat'],'PhaseLockingData');
                save([sessionInfo.session.name '.PhaseLockingData_sg.SubSession.cellinfo.mat'],'PhaseLockingData_sg');
            end
        catch
            warning('It has not been possible to run theta and gamma mod code for SubSessions...')
        end
    else
        try
            disp('Theta, gamma and HFO modulation...')
            [sessionInfo] = bz_getSessionInfo(pwd,'noPrompts',true);
            % Theta Profile
            powerProfile_theta = bz_PowerSpectrumProfile([thetaFreq(1) thetaFreq(2)], 'channels',sessionInfo.channels,'showfig',true);
            % Low Gamma Profile
            powerProfile_sg = bz_PowerSpectrumProfile([sgFreq(1) sgFreq(2)], 'channels',sessionInfo.channels,'showfig',true);
            % HFO Profile
            powerProfile_hg = bz_PowerSpectrumProfile([hgFreq(1) hgFreq(2)], 'channels', sessionInfo.channels,'showfig',true);
            % get channels of interest % max theta power above pyr layer
            [~, a] = max(powerProfile_theta.mean);
            region.CA1sp = powerProfile_theta.channels(a);

            for ii = 1:size(sessionInfo.AnatGrps,2)
                if ismember(region.CA1sp, sessionInfo.AnatGrps(ii).Channels)
                    p_th = powerProfile_theta.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                    p_hg = powerProfile_hg.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                    p_sg = powerProfile_sg.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                    channels = sessionInfo.AnatGrps(ii).Channels;
                    [~,a] = max(p_th(1:find(channels == region.CA1sp)));
                    region.CA1so = channels(a);
                    [~,a] = max(p_th(find(channels == region.CA1sp):end));
                    region.CA1slm = channels(a - 1+ find(channels == region.CA1sp));

                    figure,
                    hold on
                    plot(1:length(channels), zscore(p_th));
                    plot(1:length(channels), zscore(p_hg));
                    plot(1:length(channels), zscore(p_sg));
                    xlim([1 length(channels)]);
                    set(gca,'XTick',1:length(channels),'XTickLabel',channels,'XTickLabelRotation',45);
                    ax = axis;
                    xlabel(strcat('Channels (neuroscope)-Shank',num2str(ii))); ylabel('power (z)');
                    plot(find(channels == region.CA1sp)*ones(2,1),ax([3 4]),'-k');
                    plot(find(channels == region.CA1so)*ones(2,1),ax([3 4]),'--k');
                    plot(find(channels == region.CA1slm)*ones(2,1),ax([3 4]),'-.k');
                    legend('4-12Hz', 'hfo', '30-60','pyr','~or','~slm');
                    saveas(gcf,'SummaryFigures\regionDef.png');
                end
            end
            save([sessionInfo.FileName,'.region.mat'],'region');

            % power profile of pyr channel of all session
            lfpT = bz_GetLFP(region.CA1sp,'noPrompts',true);
            params.Fs = lfpT.samplingRate; params.fpass = [2 200]; params.tapers = [3 5]; params.pad = 1;
            [S,t,f] = mtspecgramc(single(lfpT.data),[2 1], params);
            S = log10(S); % in Db
            S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending        
            figure;
            subplot(3,4,[1 2 3 5 6 7])
            imagesc(t,f,S_det',[-1.5 1.5]);
            set(gca,'XTick',[]); ylabel('Freqs');
            subplot(3,4,[4 8]);
            plot(mean(S,1),f);
            set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
            saveas(gcf,'SummaryFigures\spectrogramAllSession.png');


        %   bz_plotPower()

            % phase modulation by spikes
            spikes = loadSpikes('getWaveformsFromDat',false);

            if diffLFPs
                filename = dir('*PhaseLockingData.diffLFPs.cellinfo.mat*'); 
                if ~isempty(filename)
                    disp(['Phase Locking Data Different LFPs already detected ! Loading file: ', filename.name]);
                    load(filename.name)
                    PLD = PhaseLockingData;
                else
                    for ii = 1:spikes.numcells
                        lfp = bz_GetLFP(spikes.maxWaveformCh1(ii)-1,'noPrompts',true);
                        PLD{ii} = bz_PhaseModulationDifferentLFPs(spikes.times{ii},lfp,[thetaFreq(1) thetaFreq(2)], 'plotting',false,'method','wavelet');
                    end
                    PhaseLockingData = PLD;
                    save([sessionInfo.session.name '.PhaseLockingData.diffLFPs.cellinfo.mat'],'PhaseLockingData');
                end
                disp('Theta modulation...')
                figure,
                set(gcf,'Position',[100 -100 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    area([PLD{jj}.phasebins; PLD{jj}.phasebins + pi*2],[PLD{jj}.phasedistros; PLD{jj}.phasedistros],'EdgeColor','none');
                    hold on
                    ax = axis;
                    x = 0:.001:4*pi;
                    y = cos(x);
                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                    xlim([0 4*pi]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('prob');
                    elseif jj == size(spikes.UID,2)
                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                        xlabel('phase (rad)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,'PhaseModulationFig\thetaModDiffLFPS.png');

                % gamma modulation
                disp('Gamma modulation...');
                filename = dir('*PhaseLockingData_sg.diffLFPs.cellinfo.mat*');
                if ~isempty(filename)
                    disp(['Phase Locking Data_sg Different LFPs already detected ! Loading file: ', filename.name]);
                    load(filename.name)
                    PLD_sg = PhaseLockingData_sg;
                else
                    for ii=1:spikes.numcells
                        lfp = bz_GetLFP(spikes.maxWaveformCh1(ii)-1,'noPrompts',true);
                        PLD_sg{ii} = bz_PhaseModulationDifferentLFPs(spikes.times{ii},lfp,[sgFreq(1) sgFreq(2)],'plotting',false,'method','wavelet');     
                        PhaseLockingData_sg = PLD_sg;
                        save([sessionInfo.session.name '.PhaseLockingData_sg.diffLFPs.cellinfo.mat'],'PhaseLockingData_sg');
                    end
                end

                figure
                set(gcf,'Position',[100 -100 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    area([PLD_sg{jj}.phasebins; PLD_sg{jj}.phasebins + pi*2],[PLD_sg{jj}.phasedistros; PLD_sg{jj}.phasedistros],'EdgeColor','none');
                    hold on
                    ax = axis;
                    x = 0:.001:4*pi;
                    y = cos(x);
                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                    xlim([0 4*pi]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('prob');
                    elseif jj == size(spikes.UID,2)
                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                        xlabel('phase (rad)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,'PhaseModulationFig\sgModDiffLFPS.png');

            else
                PLD = bz_PhaseModulation(spikes,lfpT,[thetaFreq(1) thetaFreq(2)], 'plotting',true,'method','wavelet','saveMat',true);
                disp('Theta modlation...')
                set(gcf,'Position',[100 -100 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    area([PLD.phasebins; PLD.phasebins + pi*2],[PLD.phasedistros(:,jj); PLD.phasedistros(:,jj)],'EdgeColor','none');
                    hold on
                    ax = axis;
                    x = 0:.001:4*pi;
                    y = cos(x);
                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                    xlim([0 4*pi]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('prob');
                    elseif jj == size(spikes.UID,2)
                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                        xlabel('phase (rad)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,'SummaryFigures\thetaMod.png');

                % gamma modulation
                PLD = bz_PhaseModulation(spikes,lfpT,[sgFreq(1) sgFreq(2)],'plotting',false,'method','wavelet','saveMat',true);     
                disp('Gamma modulation...');
                figure
                set(gcf,'Position',[100 -100 2500 1200]);
                for jj = 1:size(spikes.UID,2)
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    area([PLD.phasebins; PLD.phasebins + pi*2],[PLD.phasedistros(:,jj); PLD.phasedistros(:,jj)],'EdgeColor','none');
                    hold on
                    ax = axis;
                    x = 0:.001:4*pi;
                    y = cos(x);
                    y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
                    h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
                    xlim([0 4*pi]);
                    title(num2str(jj),'FontWeight','normal','FontSize',10);

                    if jj == 1
                        ylabel('prob');
                    elseif jj == size(spikes.UID,2)
                        set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                        xlabel('phase (rad)');
                    else
                        set(gca,'YTick',[],'XTick',[]);
                    end
                end
                saveas(gcf,'SummaryFigures\gammaMod.png');
            end
            % save Power Spectrum Profile theta
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_theta{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_theta{1}.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_theta');                
            % save Power Spectrum Profile sg
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_sg{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_sg{1}.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_sg');
            % save Power Spectrum Profile hg
            save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(powerProfile_hg{1}.processinginfo.params.frange(1)),'_',num2str(powerProfile_hg{1}.processinginfo.params.frange(2)),'.channelinfo.mat'],'powerProfile_hg');
            save([sessionInfo.FileName,'.region.mat'],'region');
            save([sessionInfo.session.name '.PhaseLockingData.cellinfo.mat'],'PhaseLockingData');
            save([sessionInfo.session.name '.PhaseLockingData_sg.cellinfo.mat'],'PhaseLockingData_sg');
            
        catch
            warning('It has not been possible to run theta and gamma mod code...')
        end
    end
end

%% 4- BEHAVIOUR

if any(ismember(listOfAnalysis,'behaviour'))   
    try
        disp('Computing behaviour variables...')
        tracking = getSessionTracking();
        behaviour = getSessionBehaviour_v2();
        
        % PLACE CELLS SUMMARY
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);
        firingMaps = bz_firingMapAvg(behaviour,spikes,'saveMat',true,'speedFilter',true,'periodicAnalysis',false,'spikeShuffling',false);       
    catch
        warning('It has not been possible to run Behaviour code...')
    end
end

%% 5  BEHAVIOUR PERFORMANCE

if any(ismember(listOfAnalysis,'performance'))
    try
        disp('Computing Performance for Experimental Paradigm (YMaze/TMaze) ...');
        tracking = getSessionTracking();
        performance = getSessionPerformance('tracking',tracking);
        
    catch
        warning('It has not been possible to run Performance code...')
    
    end
end
%% 6 - SPIKE TRAIN ANALYSIS

if any(ismember(listOfAnalysis,'spikeTrain'))
    try
        disp('Computing Spike Train Analysis...')
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat,'forceReload',false,'showWaveforms',showWaveforms);
        spikeTrain = bz_SpikeTrain(spikes,'analyzeSubSessions',analyzeSubSessions,'showFigure',true);
    catch
        warning('It has not been possible to run Spike Train Analysis...')
    end
end


%% 7 - CELL EXPLORER

if any(ismember(listOfAnalysis,'CellExplorer'))
    try
        disp('Running CellExplorer code...');
        session = sessionTemplate(basepath,'showGUI',false);        
        validateSessionStruct(session);        
        cell_metrics = ProcessCellMetrics('session', session,'showGUI',true);
        cell_metrics = CellExplorer('metrics',cell_metrics); 
    catch
        disp('It is not possible to run CellExplorer code...')
    end
end

%% 8 - LFP ANALYSIS
if any(ismember(listOfAnalysis,'lfp_analysis'))
    mkdir('lfpAnalysisFigures')
%     analyzeSubSessions = true;
  if analyzeSubSessions   
      try
          disp('Computing lfp_analysis for SubSessions...')
          session = loadSession;
          sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
          shanks = sessionInfo.AnatGrps;
          shanks(excludeShanks) = [];
          
          if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
              load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
              count = 1;
              
              for ii=1:size(MergePoints.foldernames,2)
                  timestamps = MergePoints.timestamps(ii,:);
                  foldername = MergePoints.foldernames{ii};
                  % Pick Channel of Ripples and SW
                  if ~isempty(dir([basepath filesep '*.channelinfo.ripples.mat']))
                     disp('Channelinfo.ripples.mat file detected ! Loading file.')
                     file = dir([basepath filesep '*channelinfo.ripples.mat'])
                     load(file.name);
                  end
                  lfp1 = bz_GetLFP(rippleChannels.Ripple_Channel,'restrict',timestamps);
                  lfp2 = bz_GetLFP(rippleChannels.Sharpwave_Channel,'restrict',timestamps);
                  % Compute Coherence between one channel of each shank
                  coherencePerShank.(foldername) = bz_coherencePerShank('timestampsSubSession',timestamps,'foldername',foldername,'saveMat',false);
                  % Compute Coherence between two channels
                  coherogram.(foldername) = bz_MTCoherogram(lfp1,lfp2,'foldername',foldername,'range',[0 200], 'window',1,'saveMat',false);
                  % Compute Cross-Frequency Coupling
                  lfp = bz_GetLFP('all','restrict',timestamps);
                  CFCPhaseAmp_lg.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[sgFreq(1):1:sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  CFCPhaseAmp_hg.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hgFreq(1):1:hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  CFCPhaseAmp_hfo.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hfoFreq(1):1:hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  CFCPhaseAmp.(foldername) = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[1:1:200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  
                  PhaseAmpCoupling_lg.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  PhaseAmpCoupling_hg.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);
                  PhaseAmpCoupling_hfo.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername,'saveMat',false);                  
                  PhaseAmpCoupling.(foldername) = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [1 200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel','foldername',foldername,'saveMat',false);                  
                  % Power-Power Correlation
                  lfp = bz_GetLFP([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel],'restrict',timestamps);
                  comodulogram.(foldername) = bz_Comodulogram(lfp,[],[],'foldername',foldername,'ch1',rippleChannels.Ripple_Channel, 'ch2',rippleChannels.Sharpwave_Channel,'saveMat',false,'saveMat',false);
                  % Cross-Spectral Matrix ¿How to plot?
                  %bz_MTCrossSpec(lfp);
                  
                  % GMI Index
                  % Maybe we can compute this index over a wide range of
                  % channels (as an input)
                  lfp = bz_GetLFP('all','restrict',timestamps);
                  GMI_lg.(foldername) = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);
%                   bz_GMI(lfp,[4 12], [30 80],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);
                  GMI_hg.(foldername) = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);
%                   bz_GMI(lfp,[4 12], [80 150],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);
                  GMI_hfo.(foldername) = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',fodername,'saveMat',false);

                  GMI.(foldername) = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [1 200],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);
%                   bz_GMI(lfp,[4 12], [1 200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);

                  % GMI Index Tort
                  % Maybe we can compute this index over a wide range of
                  % channels (as an input)
                  lfp = bz_GetLFP('all','restrict',timestamps);
                  MI_lg.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);
                  MI_hg.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)],[hgFreq(1) hgFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'fodlername',foldername,'saveMat',false);
                  MI_hfo.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);
                  MI.(foldername) = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [1 200], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels,'foldername',foldername,'saveMat',false);
              end
              
              % Saving .mat files
              % Coherence per Shank
              save([basepath filesep sessionInfo.FileName,'.Coherence_Shanks.SubSession.lfp.mat'],'coherencePerShank');
              % Coherogram
              save([basepath filesep sessionInfo.FileName,'Coherogram.SubSession.lfp.mat'],'coherogram');
              % CFCPhaseAmp lg
              save([basepath filesep sessionInfo.FileName,'CFCPhaseAmp_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_lg');
              % CFCPhaseAmp hg
              save([basepath filesep sessionInfo.FileName,'CFCPhaseAmp_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_hg');
              % CFCPhaseAmp hfo
              save([basepath filesep sessionInfo.FileName,'CFCPhaseAmp_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'CFCPhaseAmp_hfo');
              % CFCPhaseAmp 
              save([basepath filesep sessionInfo.FileName,'CFCPhaseAmp.SubSession.lfp.mat'],'CFCPhaseAmp');
              % PhaseAmpCouplingByAmp lg
              save([basepath filesep sessionInfo.FileName,'PhaseAmpCouplingByAmp_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_lg');
              % PhaseAmpCouplingByAmp hg
              save([basepath filesep sessionInfo.FileName,'PhaseAmpCouplingByAmp_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_hg');
              % PhaseAmpCouplingByAmp hfo
              save([basepath filesep sessionInfo.FileName,'PhaseAmpCouplingByAmp_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'PhaseAmpCoupling_hfo');
              % PhaseAmpCouplingByAmp
              save([basepath filesep sessionInfo.FileName,'PhaseAmpCouplingByAmp.SubSession.lfp.mat'],'PhaseAmpCoupling');
              % Comodulogram
              %save([basepath filesep sessionInfo.FileName,'Comodulogram'],'comodulogram');
              % GMI lg
              save([basepath filesep sessionInfo.FileName,'GMI_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'GMI_lg');
              % GMI hg
              save([basepath filesep sessionInfo.FileName,'GMI_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'GMI_hg');
              % GMI hfo
              save([basepath filesep sessionInfo.FileName,'GMI_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'GMI_hfo');
              % GMI
              save([basepath filesep sessionInfo.FileName,'GMI.SubSession.lfp.mat'],'GMI');
              % GMI Tort lg
              save([basepath filesep sessionInfo.FileName,'GMI_Tort_',num2str(sgFreq(1)),'_',num2str(sgFreq(2)),'.SubSession.lfp.mat'],'MI_lg');
              % GMI Tort hg
              save([basepath filesep sessionInfo.FileName,'GMI_Tort_',num2str(hgFreq(1)),'_',num2str(hgFreq(2)),'.SubSession.lfp.mat'],'MI_hg');
              % GMI Tort hfo
              save([basepath filesep sessionInfo.FileName,'GMI_Tort_',num2str(hfoFreq(1)),'_',num2str(hfoFreq(2)),'.SubSession.lfp.mat'],'MI_hfo');
              % GMI Tort 
              save([basepath filesep sessionInfo.FileName,'GMI_Tort.SubSession.lfp.mat'],'MI');
          end
          
      catch
          disp('It is not possible to perform lfp_analysis SubSession...')         
      end      
  else
      try
          disp('Computing lfp_analysis for Whole Session...')
          session = loadSession;
          sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
          shanks = sessionInfo.AnatGrps;
          shanks(excludeShanks) = [];
          
         % Pick Channel of Ripples and SW
          if ~isempty(dir([basepath filesep '*.channelinfo.ripples.mat']))
             disp('Channelinfo.ripples.mat file detected ! Loading file.')
             file = dir([basepath filesep '*channelinfo.ripples.mat'])
             load(file.name);
          end
          lfp1 = bz_GetLFP(rippleChannels.Ripple_Channel);
          lfp2 = bz_GetLFP(rippleChannels.Sharpwave_Channel); 
          % Compute Coherence between one channel of each shank
          coherencePerShank = bz_coherencePerShank();
          % Compute Coherence between two channels
          bz_MTCoherogram(lfp1,lfp2,'range',[0 200], 'window',1);
          % Compute Cross-Frequency Coupling
          lfp = bz_GetLFP('all');
          CFCPhaseAmp_lg = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[sgFreq(1):1:sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          CFCPhaseAmp_hg = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hgFreq(1):1:hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          CFCPhaseAmp_hfo = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[hfoFreq(1):1:hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          CFCPhaseAmp = bz_CFCPhaseAmp(lfp,[thetaFreq(1):1:thetaFreq(2)],[1:1:200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);

          PhaseAmpCoupling_lg = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          PhaseAmpCoupling_hg = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);
          PhaseAmpCoupling_hfo = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [hfoFreq(1) hfoFreq(2)],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);                  
          PhaseAmpCoupling = bz_PhaseAmpCouplingByAmp(lfp,[thetaFreq(1) thetaFreq(2)], [1 200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel);                  
          % Power-Power Correlation
          lfp = bz_GetLFP([rippleChannels.Ripple_Channel rippleChannels.Sharpwave_Channel]);
          comodulogram = bz_Comodulogram(lfp,[],[],'ch1',rippleChannels.Ripple_Channel, 'ch2',rippleChannels.Sharpwave_Channel);
                  

          % GMI Index
          % Maybe we can compute this index over a wide range of
          % channels (as an input)
          lfp = bz_GetLFP('all');
          GMI_lg = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
%                   bz_GMI(lfp,[4 12], [30 80],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);
          GMI_hg = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [hgFreq(1) hgFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
%                   bz_GMI(lfp,[4 12], [80 150],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);
          GMI_hfo = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)],[hfoFreq(1) hfoFreq(2)],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
          GMI = bz_GMI(lfp,[thetaFreq(1) thetaFreq(2)], [1 200],'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
%                   bz_GMI(lfp,[4 12], [1 200],'phaseCh',rippleChannels.Ripple_Channel,'ampCh',rippleChannels.Ripple_Channel,'foldername',foldername);

          % GMI Index Tort
          % Maybe we can compute this index over a wide range of
          % channels (as an input)
          lfp = bz_GetLFP('all');
          MI_lg = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [sgFreq(1) sgFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
          MI_hg = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)],[hgFreq(1) hgFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
          MI_hfo = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)],[hfoFreq(1) hfoFreq(2)], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);
          MI = bz_GMITort(lfp,[thetaFreq(1) thetaFreq(2)], [1 200], 'phaseCh',sessionInfo.channels,'ampCh',sessionInfo.channels);    
      catch
          disp('It is not possible to run lfp_analysis for Whole Session..')
      end  
  end    
end

%% 9 - CREATION EXCEL FILE 

if any(ismember(listOfAnalysis,'excel'))
    if analyzeSubSessions
        try
            disp('Creating Excel file for SubSessions ...')
            excel = bz_createExcelUbSessions('basepath',basepath,'diffLFPs',diffLFPs,'thetaFreq',thetaFreq,'sgFreq',sgFreq,'hgFreq',hgFreq);
        catch
            warning('It has not been possible to create Excel file for SubSession...')
        end
    else
        try
            disp('Creating Excel file...')
            excel = bz_createExcel('basepath',basepath,'diffLFPs',diffLFPs,'thetaFreq',thetaFreq,'sgFreq',sgFreq,'hgFreq',hgFreq);
        catch
            warning('It has not been possible to create Excel file...')
        end
    end
end
%% 10 - PLOT PLACE FIELDS
try
    disp('Plotting Place Fields')
%     plot_placeFields(firingMaps,spikes,tracking,cell_metrics);
    plot_placeFields(firingMaps,spikes,tracking);
catch
    disp('It is not possible to run plot Place Fields')
end


