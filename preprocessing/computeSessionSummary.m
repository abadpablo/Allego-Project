
function computeSessionSummary(varargin)
% Runs a list of preliminary descriptive analysis 
%
% USAGE 
%   computeSessionSummary(varargin)
%
% INPUT (optional)
%   basepath         By default pwd
%   listOfAnalysis   Summary analysis that will be computed (as a cell with
%                        strings). Default 'all'. Possible values: 'spikes',
%                        'analogPulses', 'digitalPulses',
%                        'downStates','ripples', 'tMazeBehaviour',
%                        'linearMazeBehaviour', 'thetaModulation'
%   exclude         Cell with strings with the list of analysis to exclude
%  
%   (specific analysis optionss)
%   excludeShanks           Default []
%   getWaveformsFromDat     From 'spikes' summary, default false.
%   analogChannelsList      Array of channel to perform 'analogPulses' psth and csd. Default 'all'
%   digitalChannelsList     Array of channel to perform 'digitalPulses' psth and csd. Default 'all'
%
% Manu Valero-BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'listOfAnalysis','all',@iscellstr);
addParameter(p,'exclude',[],@iscellstr);
addParameter(p,'excludeShanks',[],@isnumeric);
addParameter(p,'getWaveformsFromDat',false,@islogical);
addParameter(p,'analogChannelsList','all',@isnumeric);
addParameter(p,'digitalChannelsList','all',@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
listOfAnalysis = p.Results.listOfAnalysis;
exclude = p.Results.exclude;
excludeShanks = p.Results.excludeShanks;
getWaveformsFromDat = p.Results.getWaveformsFromDat;
analogChannelsList = p.Results.analogChannelsList;
digitalChannelsList = p.Results.digitalChannelsList;

prevPath = pwd;
cd(basepath);

if ischar(listOfAnalysis) && strcmpi(listOfAnalysis,'all')
%     listOfAnalysis = {'spikes', 'analogPulses', 'digitalPulses', 'downStates', 'ripples', 'tMazeBehaviour','linearMazeBehaviour','thetaModulation'};
    listOfAnalysis = {'spikes', 'analogPulses', 'digitalPulses', 'downStates', 'ripples', 'tMazeBehaviour','linearMazeBehaviour','YMazeBehaviour','OpenFieldBehaviour','thetaModulation'};
end
if ~isempty(exclude)
    listOfAnalysis(ismember(listOfAnalysis, exclude)) = [];
end
session = loadSession;
excludeChannels = [];

for ii = 1:length(excludeShanks)
    excludeChannels = session.extracellular.electrodeGroups.channels{excludeShanks(ii)};
end

mkdir('SummaryFigures'); % create folder
close all

% SPIKES SUMMARY
if any(ismember(listOfAnalysis,'spikes'))
    try
        disp('Spike-waveform, ACG and cluster location...');
        spikes = loadSpikes('getWaveformsFromDat',getWaveformsFromDat);

        % plot spikes summary
        disp('Plotting spikes summary...');
        if exist('chanMap.mat','file')
            load('chanMap.mat','xcoords','ycoords');
            xcoords = xcoords - min(xcoords); xcoords = xcoords/max(xcoords);
            ycoords = ycoords - min(ycoords); ycoords = ycoords/max(ycoords); 
            xcoords(excludeChannels) = [];
            ycoords(excludeChannels) = [];
        else
            xcoords = NaN;
            ycoords = NaN;
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
        saveas(gcf,'SummaryFigures\spikes.png');

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
        saveas(gcf,'SummaryFigures\CrossCorr.png'); 
    catch
        warning('Error on Spike-waveform, autocorrelogram and cluster location! ');
    end
end

% ANALOG CHANNELS PSTH AND CSD SUMMARY
if any(ismember(listOfAnalysis,'analogPulses'))
    try 
        disp('Psth and CSD from analog-in inputs...');
        % copy NPY file to Kilosort folder, for phy psth option
        pulses = bz_getAnalogPulses('analogCh',analogChannelsList,'manualThr',false);
        writeNPY(pulses.timestamps(:,1),'pulTime.npy'); % if error, transpose!
        kilosort_path = dir('*Kilosort*');
        try copyfile('pulTime.npy', strcat(kilosort_path.name,'\','pulTime.npy')); end% copy pulTime to kilosort folder
        xml = LoadParameters;
        spikes = loadSpikes('getWaveformsFromDat',false);
        
        for mm = 1:length(analogChannelsList)
            fprintf('Stimulus %3.i of %3.i \n',mm, length(analogChannelsList)); %\n
            st = pulses.timestamps(pulses.analogChannel == analogChannelsList(mm),1);
            if isempty(st)
                st = pulses.timestamps(pulses.analogChannel == analogChannelsList(mm) + 1,1);
            end
            % CSD
            shanks = xml.AnatGrps;
            shanks(excludeShanks) = [];
            figure
            set(gcf,'Position',[100 100 1400 600])
            for jj = 1:length(shanks)
                lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts', true);
                twin = 0.02;
                [csd,lfpAvg] = bz_eventCSD(lfp,st,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                taxis = linspace(-twin,twin,size(csd.data,1));
                cmax = max(max(csd.data)); 
                subplot(1,size(xml.AnatGrps,2),jj);
                contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('AnalogCh #',num2str(analogChannelsList(mm))),'FontWeight','normal'); 
                colormap jet; try caxis([-cmax cmax]); end
                hold on
                for kk = 1:size(lfpAvg.data,2)
                    plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                end
            end
            saveas(gcf,['SummaryFigures\analogPulsesCSD_ch',num2str(analogChannelsList(mm)), '.png']);

            % PSTH
            win = [-0.1 0.5];
            if length(st) > 1000
                st = randsample(st, 1000);
                st = sort(st);
            end
            disp('Plotting spikes raster and psth...');
            spikeResponse = [];
            [stccg, t] = CCG([spikes.times st],[],'binSize',0.005,'duration',1);
            figure;
            set(gcf,'Position',[100 -100 2500 1200]);
            for jj = 1:size(spikes.UID,2)
                fprintf(' **Pulses from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                rast_x = []; rast_y = [];
                for kk = 1:length(st)
                    temp_rast = spikes.times{jj} - st(kk);
                    temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                    rast_x = [rast_x temp_rast'];
                    rast_y = [rast_y kk*ones(size(temp_rast))'];
                end
                
                spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,jj)))'];
                subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                plot(rast_x, rast_y,'.','MarkerSize',1)
                hold on
                plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),end,jj) * kk/max(stccg(:,end,jj))/2,'k','LineWidth',2);
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
            saveas(gcf,['SummaryFigures\AnalogPulsesRaster_ch',num2str(analogChannelsList(mm)) ,'ch.png']); 
            
            figure
            imagesc([t(1) t(end)],[1 size(spikeResponse,1)], spikeResponse); caxis([-3 3]); colormap(jet);
            set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([-.1 t(end)]);
            saveas(gcf,['SummaryFigures\AnalogPulsesPsth_',num2str(analogChannelsList(mm)) ,'ch.png']); 
        end  

    catch
        warning('Error on Psth and CSD from analog-in inputs! ');
    end
end

% DIGITAL CHANNELS PSTH AND CSD SUMMARY
if any(ismember(listOfAnalysis,'digitalPulses'))
    try 
        disp('Psth and CSD from digital-in inputs...');
        xml = LoadParameters;
        digitalIn = bz_getDigitalIn('all','fs',xml.rates.wideband);
        spikes = loadSpikes('getWaveformsFromDat',false);
        if ischar(digitalChannelsList) && strcmpi(digitalChannelsList,'all')
            digitalChannelsList = 1:length(digitalIn.timestampsOn);
        end
        
        for mm = 1:length(digitalChannelsList)
            fprintf('Stimulus %3.i of %3.i \n',mm, length(digitalChannelsList)); %\n
            st = digitalIn.timestampsOn{digitalChannelsList(mm)};
            % CSD
            shanks = xml.AnatGrps;
            shanks(excludeShanks) = [];
            figure
            set(gcf,'Position',[100 100 1400 600])
            for jj = 1:length(shanks)
                lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts', true);
                twin = 0.02;
                [csd,lfpAvg] = bz_eventCSD(lfp,st,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                taxis = linspace(-twin,twin,size(csd.data,1));
                cmax = max(max(csd.data)); 
                subplot(1,size(xml.AnatGrps,2),jj);
                contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('DigitalCh #',num2str(analogChannelsList(mm))),'FontWeight','normal'); 
                colormap jet; try caxis([-cmax cmax]); end
                hold on
                for kk = 1:size(lfpAvg.data,2)
                    plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                end
            end
            saveas(gcf,['SummaryFigures\digitalPulsesCSD_ch',num2str(digitalChannelsList(mm)), '.png']);

            % PSTH
            figure;
            set(gcf,'Position',[100 -100 2500 1200])
            
            if ~isempty(st)
                win = [-0.1 0.5];
                if length(st) > 1000
                    st = randsample(st, 1000);
                    st = sort(st);
                end

                disp('Plotting spikes raster and psth...');
                spikeResponse = [];
                [stccg, t] = CCG([spikes.times st],[],'binSize',0.005,'duration',1);

                for jj = 1:size(spikes.UID,2)
                    fprintf(' **Pulses from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
                    rast_x = []; rast_y = [];
                    for kk = 1:length(st)
                        temp_rast = spikes.times{jj} - st(kk);
                        temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
                        rast_x = [rast_x temp_rast'];
                        rast_y = [rast_y kk*ones(size(temp_rast))'];
                    end

                    spikeResponse = [spikeResponse; zscore(squeeze(stccg(:,end,jj)))'];
                    subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
                    plot(rast_x, rast_y,'.','MarkerSize',1)
                    hold on
                    plot(t(t>win(1) & t<win(2)), stccg(t>win(1) & t<win(2),end,jj) * kk/max(stccg(:,end,jj))/2,'k','LineWidth',2);
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
            end
            saveas(gcf,['SummaryFigures\digitalPulsesRaster_ch',num2str(digitalChannelsList(mm)) ,'ch.png']); 
            
            figure
            if ~isempty(st)
                imagesc([t(1) t(end)],[1 size(spikeResponse,1)], spikeResponse); caxis([-3 3]); colormap(jet);
                set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells'); xlim([-.1 t(end)]);
            end
            saveas(gcf,['SummaryFigures\digitalPulsesPsth_',num2str(digitalChannelsList(mm)) ,'ch.png']); 
        end  

    catch
        warning('Error on Psth and CSD from digital-in inputs! ');
    end
end

% DOWN-STATES
if any(ismember(listOfAnalysis,'downStates'))
    try
        disp('Slow-waves CSD and PSTH...');
        UDStates = detectUD;
        xml = LoadParameters;
        shanks = xml.AnatGrps;
        shanks(excludeShanks) = [];
        % CSD
        twin = 0.2;
        evs = UDStates.timestamps.DOWN(:,1);
        figure
        set(gcf,'Position',[100 100 1400 600])
        for jj = 1:size(shanks,2)
            lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts', true);
            [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
            taxis = linspace(-twin,twin,size(csd.data,1));
            cmax = max(max(csd.data)); 
            subplot(1,size(shanks,2),jj);
            contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
            set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('DOWN-UP, Shank #',num2str(jj)),'FontWeight','normal'); 
            colormap jet; caxis([-cmax cmax]);
            hold on
            for kk = 1:size(lfpAvg.data,2)
                plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
            end
        end
        saveas(gcf,'SummaryFigures\downUpStatesCSD.png');

        % PSTH
        st = UDStates.timestamps.DOWN;
        spikeResponse = [];
        win = [-0.2 0.2];
        figure
        set(gcf,'Position',[100 -100 2500 1200])
        for jj = 1:size(spikes.UID,2)
            fprintf(' **UD from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
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
        saveas(gcf,'SummaryFigures\downUpStatesRaster.png'); 
        
        figure
        imagesc([t(1) t(end)],[1 size(spikeResponse,2)], spikeResponse); caxis([-3 3]); colormap(jet);
        set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');
        saveas(gcf,['SummaryFigures\downUpStatesRasterPsth.png']); 
    catch
        warning('Error on Psth and CSD from down-states!');
    end
end

% RIPPLES
if any(ismember(listOfAnalysis,'ripples'))
    try
        disp('Ripples CSD and PSTH...');
        rippleChannels = computeRippleChannel('discardShanks',excludeShanks);
%         ripples = bz_DetectSWR([rippleChannels.Ripple_Channel, rippleChannels.Sharpwave_Channel],'saveMat',true);
        ripples = bz_FindRipples(basepath, rippleChannels.Ripple_Channel);
        
        xml = LoadParameters;
        shanks = xml.AnatGrps;
        shanks(excludeShanks) = [];
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
        saveas(gcf,'SummaryFigures\ripplesCSD.png');

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
        saveas(gcf,'SummaryFigures\ripplesRaster.png'); 
        
        figure
        imagesc([t(1) t(end)],[1 size(spikeResponse,2)], spikeResponse); caxis([-3 3]); colormap(jet);
        xlim([-.2 .2]); set(gca,'TickDir','out'); xlabel('Time'); ylabel('Cells');
        saveas(gcf,['SummaryFigures\ripplesPsth.png']); title('Ripples');
    catch
        warning('Error on Psth and CSD from ripples! ');
    end
end

% TMAZEBEHAVIOUR AND LINEARMAZEBEHAVIOUR
if any(ismember(listOfAnalysis,'tMazeBehaviour')) || any(ismember(listOfAnalysis,'linearMazeBehaviour'))
   try 
%         getSessionTracking('convFact',0.1149,'roiTracking','manual');
        getSessionTracking();
        if any(ismember(listOfAnalysis,'tMazeBehaviour'))
            getSessionArmChoice;
        end
        
%         behaviour = getSessionLinearize; 
        behaviour = getSessionLinearize_v2;        
        % PLACE CELLS SUMMARY
        spikes = loadSpikes('getWaveformsFromDat',false);
        firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
        placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps);
   catch
       warning('It has not been possible to run the behaviour code...');
   end
end



% YMAZE AND OPEN FIELD BEHAVIOUR (NO LINEARIZE AND NO COMPUTING TRIALS)


if any(ismember(listOfAnalysis,'YMazeBehaviour')) || any(ismember(listOfAnalysis,'OpenFieldBehaviour'))
    
    try 
%         getSessionTracking('convFact',0.1149,'roiTracking','manual');
        getSessionTracking();     
        
        behavior = getSessionBehaviour;        
        % PLACE CELLS SUMMARY
        spikes = loadSpikes('getWaveformsFromDat',false);
        firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
%         firingMaps = bz_firingMapAvg_v2(behavior,spikes,'saveMat',true);
        placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps);      
   catch
       warning('It has not been possible to run the behaviour code...');
    end
end


% THETA AND GAMMA PHASE MODULATION
if any(ismember(listOfAnalysis,'thetaModulation'))
    try 
        disp('Theta modulation...');
        % Theta profile
        xml = LoadParameters;
        channels = xml.channels; channels(excludeChannels) = [];
        powerProfile_theta = bz_PowerSpectrumProfile([6 12],'channels',channels,'showfig',true); % [0:63]
        
        % max theta power above pyr layer
        rippleChannels = computeRippleChannel('discardShanks',excludeShanks);
        for ii = 1:length(xml.AnatGrps)
            if any(find(xml.AnatGrps(ii).Channels == rippleChannels.Ripple_Channel))
                rippleShank = ii;
            end
        end
        [~, channels_ripleShank] = intersect(powerProfile_theta.channels, xml.AnatGrps(rippleShank).Channels);
        thetaProfile_rippleShank = powerProfile_theta.mean(channels_ripleShank);
        [~, indx_channel] = max(thetaProfile_rippleShank(1:find(channels_ripleShank == rippleChannels.Ripple_Channel)));
        thetaChannel = channels_ripleShank(indx_channel);

        lfpT = bz_GetLFP(thetaChannel,'noPrompts',true);
        
        % thetaMod modulation
        spikes = loadSpikes('getWaveformsFromDat',false);
        PLD = bz_PhaseModulation(spikes,lfpT,[6 12],'plotting',false,'method','wavelet');  
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
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,'SummaryFigures\thetaPhaseModulation.png');

        % gamma modulation
        PLD = bz_PhaseModulation(spikes,lfpT,[30 60],'plotting',false,'method','wavelet');  
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
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,'SummaryFigures\gamma30_60HzPhaseModulation.png');
        
        % spectrogram
        params.Fs = lfpT.samplingRate; params.fpass = [2 120]; params.tapers = [3 5]; params.pad = 1;
        [S,t,f] = mtspecgramc_fast(single(lfpT.data),[2 1],params);
        S = log10(S); % in Db
        S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending

        figure;
        subplot(1,5,1:4)
        imagesc(t,f,S_det',[-1.5 1.5]);
        set(gca,'XTick',[]); ylabel('Freqs');
        subplot(1,5,5);
        plot(mean(S,1),f);
        set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
        ylim([f(1) f(end)]);
        saveas(gcf,'SummaryFigures\spectrogramAllSession.png');    
        
    catch
        warning('It has not been possible to run theta and gamma mod code...');
    end
end

cd(prevPath);
end