% sessionSummary_abad.m

% This scritp runs preliminary descriptive analysis to get a summary overview
% of the session.

% HISTORY:
% - Based on Manu Valero-BuzsakiLab 2019
% - Some reorganization, pre-processing steps transferred to expPipeline: 5/20, AntonioFR


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir('SummaryFigures'); % create folder
close all

% manually adjust analog pulse detection
if ~isempty(analogCh)
	[pulses] = bz_getAnalogPulses('analogCh',analogCh,'manualThr',false);
else
    pulses.intsPeriods = [];
end

% get SleepStateScoring
targetFile = dir('*SleepState.states.mat');
if isempty(targetFile)
    try 
        SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods); % try to sleep score
    catch
        disp('Problem with SleepScore skyping...');
    end
end

%% 1. Spike-waveform, ACG, spatial position and CCG
if ~exist('SummaryFigures\spikes AutoCorr.png') && ~exist('SummaryFigures\spikes CrossCorr.png')
    try
        disp('Spike-waveform, ACG and cluster location...');
    %     spikes = loadSpikes;
%         spikes = bz_LoadPhy('noPrompts',true,'verbose',true,'nWaveforms',200);
        spikes = loadSpikes('getWaveformsFromDat',true,'forceReload',true);
        % plot spikes summary
        disp('Plotting spikes summary...');
        if exist('chanMap.mat','file')
    %         load('chanMap.mat','xcoords','ycoords');
            load('chanMap.mat');
            xcoords = xcoords - min(xcoords); xcoords = xcoords/max(xcoords);
            ycoords = ycoords - min(ycoords); ycoords = ycoords/max(ycoords); 
        else
            xcoords = NaN;
            ycoords = NaN;
        end
        figure
        [ccg,t] = CCG_CellExplorer(spikes.times,[],'binSize',0.001,'duration',0.08);

        for jj = 1:size(spikes.UID,2)
            fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
            dur = 0.08;
            set(gcf,'Position',[100 -100 2500 1200])
            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
            xt = linspace(-dur/2*1000,dur/2*1000,size(ccg,1));
            area(xt,ccg(:,jj,jj),'LineStyle','none'); % autocorrelogram
            try 
                xt2 = linspace(dur/3*1000+5,dur/2*1000+5+80,size(spikes.filtWaveform{jj},1)); % waveform
                spk = spikes.filtWaveform{jj} - min(spikes.filtWaveform{jj}); spk = spk/max(spk) * max(ccg(:,jj,jj));
                hold on
                plot(xt2,spk)
                hold on

                plot((xcoords*30) + dur/2*1000+5+60, ycoords*max(ccg(:,jj,jj)),'.','color',[.8 .8 .8],'MarkerSize',5); % plotting xml
                plot((xcoords(find(spikes.maxWaveformCh(jj)==sessionInfo.channels))*30) + dur/2*1000+5+60,...
                    ycoords(find(spikes.maxWaveformCh(jj)==sessionInfo.channels))*max(ccg(:,jj,jj)),'.','color',[.1 .1 .1],'MarkerSize',10); % plotting xml
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
        [allCcg, t] = CCG_CellExplorer(spikes.times,[],'binSize',0.005,'duration',0.6);
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
    
    % OverlapTool
%     disp('Computing Overlaping Spikes...')
%     overlap = bz_OverlapSpikes(spikes);
    
end

%% 2. Ripples CSD and PSTH
if ~exist([basename,'.ripples.events.mat'])
    try
        disp('Ripples CSD and PSTH...');
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        lfp = bz_GetLFP(sessionInfo.channels,'noPrompts', true);
        spikes = bz_LoadPhy('noPrompts',true);
        if isempty(dir('*.ripples.events.mat'))
            % ripples = bz_FindRipples(lfp.data,lfp.timestamps,'restrict',pulses.intsPeriods,'saveMat',true);

    %         ripples = bz_FindRipples(lfp.data(:,bz_GetBestRippleChan(lfp)),lfp.timestamps,'restrict',pulses.intsPeriods,...
    %           'noise',lfp.data(:,end),'saveMat',true);

            ripples = bz_FindRipples(lfp.data(:,bz_GetBestRippleChan(lfp)),lfp.timestamps,'show','on','saveMat',true);

            passband = [130 200];
            filtered = bz_Filter(double(lfp.data(:,bz_GetBestRippleChan(lfp))),'filter','butter','passband',passband,'order',3);
            [maps,data,stats] = bz_RippleStats(filtered,lfp.timestamps,ripples);
            ripples.maps = maps;
            ripples.data = data;
            ripples.stats = stats;
            save(fullfile(basepath, [basename '.ripples.events.mat']),'ripples')
    %         bz_SaveRippleEvents
            bz_PlotRippleStats(maps,data,stats);

    %         ripples = restrictEvents(ripples,pulses.intsPeriods);

            try ripples = rmfield(ripples,'detectorinfo'); end % remove field with the complete lfp

            fileRip = dir('*.ripples.events.mat'); save(fileRip.name,'ripples');
        else
            fileRip = dir('*.ripples.events.mat'); 
            load(fileRip.name,'ripples');
            bz_PlotRippleStats(ripples.maps,ripples.data,ripples.stats);

        end

        % CSD
        if exist('ripples') && ~isempty(ripples.timestamps)
            twin = 0.1;
            figure;
            for jj = 1:size(sessionInfo.AnatGrps,2)
                lfp = bz_GetLFP(sessionInfo.AnatGrps(jj).Channels(1:numel(sessionInfo.AnatGrps(jj).Channels)-mod(numel(sessionInfo.AnatGrps(jj).Channels),8)),'noPrompts', true);
                [csd,lfpAvg] = bz_eventCSD(lfp,ripples.peaks,'twin',[twin twin],'plotLFP',true,'plotCSD',true);
                taxis = linspace(-0.2,0.2,size(csd.data,1));
                cmax = max(max(csd.data)); 
                subplot(1,size(sessionInfo.AnatGrps,2),jj);
                contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('RIPPLES, Shank #',num2str(jj)),'FontWeight','normal'); 
                colormap jet; if ~isnan(cmax); caxis([-cmax cmax]); end
                hold on
                for kk = 1:size(lfpAvg.data,2)
                    plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k');
                end
            end
            saveas(gcf,'SummaryFigures\ripples.png');

            % PSTH
            st = ripples.timestamps(:,1);
            spikeResponse = [];
            win = [-0.2 0.2];
            figure
            set(gcf,'Position',[100 -100 2500 1200]);
            for jj = 1:size(spikes.UID,2)
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
            saveas(gcf,'SummaryFigures\psthRip.png');
            figure
            imagesc([t(1) t(end)],[1 size(spikeResponse,2)], spikeResponse); caxis([-3 3]); colormap(jet);
            xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');
            saveas(gcf,['SummaryFigures\ripplesPsth.png']); title('Ripples');
        end     
    catch
        warning('Error on Psth and CSD from analog-in inputs! ');
    end
end
%% 3. THETA AND GAMMA PROFILE/ MODULATION

if ~exist([basename,'.PhaseLockingData.cellinfo.mat'])
    try 
        disp('Theta, gamma and hfo modulation...');
        % Theta profile
        powerProfile_theta = bz_PowerSpectrumProfile([6 12],'channels',sessionInfo.channels,'showfig',true); % [0:63]
        powerProfile_sg = bz_PowerSpectrumProfile([30 60],'channels',sessionInfo.channels,'showfig',true);
        powerProfile_hfo = bz_PowerSpectrumProfile([120 250],'channels',sessionInfo.channels,'showfig',true);
        % get channels of interest % max theta power above pyr layer
        [~, a] = max(powerProfile_hfo.mean);
        region.CA1sp = powerProfile_hfo.channels(a);

        for ii = 1:size(sessionInfo.AnatGrps,2)
            if ismember(region.CA1sp, sessionInfo.AnatGrps(ii).Channels)
                p_th = powerProfile_theta.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                p_hfo = powerProfile_hfo.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                p_gs = powerProfile_sg.mean(sessionInfo.AnatGrps(ii).Channels + 1);
                channels = sessionInfo.AnatGrps(ii).Channels;
                [~,a] = max(p_th(1:find(channels == region.CA1sp)));
                region.CA1so = channels(a);
                [~,a] = max(p_th(find(channels == region.CA1sp):end));
                region.CA1slm = channels(a - 1+ find(channels == region.CA1sp));

                figure,
                hold on
                plot(1:length(channels), zscore(p_th));
                plot(1:length(channels), zscore(p_hfo));
                plot(1:length(channels), zscore(p_gs));
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
        % Spectogram all session
        lfpT = bz_GetLFP(region.CA1so,'noPrompts',true);
        params.Fs = lfpT.samplingRate; params.fpass = [2 120]; params.tapers = [3 5]; params.pad = 1;
        [S,t,f] = mtspecgramc(single(lfpT.data),[2 1],params);
        S = log10(S); % in Db
        S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending

        figure;
        subplot(3,4,[1 2 3 5 6 7])
        imagesc(t,f,S_det',[-1.5 1.5]);
        set(gca,'XTick',[]); ylabel('Freqs');
        subplot(3,4,[4 8]);
        plot(mean(S,1),f);
        set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
        subplot(3,4,[9 10 11]);
        % to do
    %     tracking = getSessionTracking();
    %     plot(tracking.position.x,tracking.position.y);
        xlabel('s'); ylabel('mov (cm/s)');
        saveas(gcf,'SummaryFigures\spectrogramAllSession.png');
        
        % phase modulation by spikes
        spikes = bz_LoadPhy('noPrompts',true);
        for ii = 1:size(spikes.filtWaveform,2)
            try spikes.region{ii} = sessionInfo.region(spikes.maxWaveformCh(ii));
            catch
                spikes.region{ii} = 'NA';
            end
        end
        save([sessionInfo.FileName '.spikes.cellinfo.mat'],'spikes');

        % thetaMod modulation
        PLD = bz_PhaseModulation(spikes,lfpT,[4 12],'plotting',false,'powerThresh',1,'saveMat',true);  
        disp('Theta modulation...');
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
        saveas(gcf,'SummaryFigures\thetaMod.png');

        % gamma modulation
        PLD = bz_PhaseModulation(spikes,lfpT,[30 60],'plotting',false,'powerThresh',2);     
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
    catch
        warning('It has not been possible to run theta and gamma mod code...');
    end
end


%% 6. BEHAVIOUR
try 
    tracking = getSessionTracking(); 
    behavior = getSessionBehaviour('maze','Open Field');
    % PLACE CELLS SUMMARY
    spikes = loadSpikes('getWaveformsFromDat',true,'forceReload',false);
%     spikes = bz_LoadPhy('noPrompts',true,'verbose',true,'nWaveforms',200);
    firingMaps = bz_firingMapAvg(behavior, spikes,'saveMat',true,'speedFilter',true);
    
    spikeTrain = bz_SpikeTrain(spikes);
       
catch
    warning('It has not been possible to run the behaviour code...');
end