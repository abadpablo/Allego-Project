function plot_placeFields(varargin)

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'firingMaps',[],@isstruct);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'cell_metrics',[],@isstruct);
addParameter(p,'spikeTrain',[],@isstruct);

parse(p,varargin{:})

basepath = p.Results.basepath;
firingMaps = p.Results.firingMaps;
spikes = p.Results.spikes;
tracking = p.Results.tracking;
cell_metrics = p.Results.cell_metrics;
spikeTrain = p.Results.spikeTrain;

% Load sessionInfo and session
sessionInfo = bz_getSessionInfo();
session = loadSession;

mkdir('placeFields')
conditions = size(tracking.events.subSessions,1);

%% Need to include only spikeTrain analysis for folders with tracking

for i=1:length(tracking.folders)
    folder = tracking.folders{i};    
    spkTrain{i} = spikeTrain.(folder);
end

for i=1:length(tracking.folders)
    timestamps{i} = tracking.events.subSessions;
end


if conditions == 1

        
for i=1:length(firingMaps.rateMaps)          
        % Let's choose color for percentiles of spatial coherence and
        % skaggs
        % Spatial Coherence
%         if firingMaps.stats{i}{1}.spatialCoherencev2.valor_r(1,2) >= firingMaps.shuffling{i}{1}.prctile99.sc
%             color_sc2 = [1 0 0]; % red
%         elseif firingMaps.stats{i}{1}.spatialCoherencev2.valor_r(1,2) < firingMaps.shuffling{i}{1}.prctile95.sc
%             color_sc2 = [0 0.5 0.8]; %blue
%         elseif firingMaps.stats{i}{1}.spatialCoherencev2.valor_r(1,2) < firingMaps.shuffling{i}{1}.prctile95.sc
%             color_sc2 = [0.5 0.5 0.5]; % grey
%         end
        
        if firingMaps.stats{i}{1}.spatialCorr.sc >= firingMaps.shuffling{i}{1}.prctile99.sc
            color_sc = [1 0 0]; % red
        elseif firingMaps.stats{i}{1}.spatialCorr.sc < firingMaps.shuffling{i}{1}.prctile95.sc
            color_sc = [0 0.5 0.8]; %blue
        elseif firingMaps.stats{i}{1}.spatialCorr.sc < firingMaps.shuffling{i}{1}.prctile95.sc
            color_sc = [0.5 0.5 0.5]; % grey
        end
        
        % Skaggs
        if firingMaps.stats{i}{1}.skaggs.bitsPerSpike >= firingMaps.shuffling{i}{1}.prctile99.bitsPerSpike
            color_skaggs = [1 0 0]; % red
        elseif firingMaps.stats{i}{1}.skaggs.bitsPerSpike > firingMaps.shuffling{i}{1}.prctile95.bitsPerSpike
            color_skaggs = [0 0.5 0.8]; % blue
        elseif firingMaps.stats{i}{1}.skaggs.bitsPerSpike < firingMaps.shuffling{i}{1}.prctile95.bitsPerSpike
            color_skaggs = [0.5 0.5 0.5]; % grey
        end
        
        
            
        
        figure
        set(gcf,'Position',[100 -100 2500 1200])
        subplot(4,4,1)
        plot(spikes.filtWaveform{i})
        title(['Cell:' , num2str(i), ' Shank:' ,num2str(spikes.shankID(i))])
        xlabel('Filtered Waveform')

        subplot(4,4,2)
        %ACG
        area(cell_metrics.acg.narrow(:,i),'LineStyle','none')
        title(['Cell type: ', num2str(cell_metrics.putativeCellType{i})])
        xlabel('Autocorrelogram')

        subplot(4,4,3)
        % ISI
        area(cell_metrics.isi.log10(:,i),'LineStyle','none')
        plot(cell_metrics.isi.log10(:,i))
        title(['Burst Index: ', num2str(spkTrain{1}{i}.BurstIndex)])
        xlabel('Interspike Interval (msec)')

        subplot(4,4,4)
        % Trilateralization
        plot(cell_metrics.general.chanCoords.x,cell_metrics.general.chanCoords.y,'.k'), hold on
        plot(cell_metrics.trilat_x(i),cell_metrics.trilat_y(i),'ob'), xlabel('x position (µm)'), ylabel('y position (µm)')
        hold off
        title([ ' Shank:' ,num2str(spikes.shankID(i)),])
        
        subplot(4,4,5)
        bz_plotStability(spikes,i,'timestamps',timestamps{1})
        title('Stability')
        
        hax =subplot(4,4,6)
%         bz_plotISI(spikes,i,'timestamps',timestamps{1},'isi_maxlag',0.05)
        % Theta modulation
        if ~isempty(dir([basepath filesep sessionInfo.FileName,'*.channelinfo.ripples.mat']))
            file = dir([basepath filesep sessionInfo.FileName,'*.channelinfo.ripples.mat']);
            load(file.name);
        end
%         lfp = bz_GetLFP(rippleChannels.Ripple_Channel);
%         pld = bz_PhaseModulation(spikes,lfp,[4 12],'intervals',timestamps{1},'unit',i);
%         pld = bz_PhaseModulation_abad(spikes,lfp,[4 12], 'intervals',timestamps{1},'unit',i);
%         pld_gam = bz_PhaseModulation_abad(spikes,lfp,[30 65], 'intervals',timestamps{1},'unit',i);
        lfp = bz_GetLFP(rippleChannels.Ripple_Channel,'restrict',timestamps{1});
        [pld,pld_gam] = bz_SpikeLFP(spikes,lfp,[4 12], 'intervals',timestamps{1},'unit',i,'nBins',90);
        x = [pld.phasebins pld.phasebins+pld.phasebins(end)];
        y = [pld.hist pld.hist];
        bar(x,y)
        hold on;
        x_theta = [0:720];
        y_theta = sin(pi/180*[0:720])*0.05*max(pld.hist)+0.85*max(pld.hist);
        plot(x_theta,y_theta)
        xlim([0 720])
        set(hax,'XTick',[0 90 180 270 360 450 540 630 720])
        ylabel('Number of Spikes')
        xlabel('Theta Degress º')
        title(['mvl: ', num2str(pld.mvl)])
        
        
        hax = subplot(4,4,7)
%         bz_plotISI(spikes,i,'timestamps',timestamps{1},'isi_maxlag',0.5)
        x = [pld_gam.phasebins pld_gam.phasebins+pld_gam.phasebins(end)];
        y = [pld_gam.hist pld_gam.hist];
        bar(x,y)
        hold on;
        x_theta = [0:720];
        y_theta = sin(pi/180*[0:720])*0.05*max(pld.hist)+0.85*max(pld.hist);
        plot(x_theta,y_theta)
        xlim([0 720])
        set(hax,'XTick',[0 90 180 270 360 450 540 630 720])
        ylabel('Number of Spikes')
        xlabel('Gamma Degress º')
        title(['mvl: ', num2str(pld_gam.mvl)])

        subplot(4,4,8)
        plot(tracking.position.x,tracking.position.y,'LineWidth',1,'Color',[0.5 0.5 0.5])
        view(0,-90)
        set(gca,'DataAspectRatio',[1 1 1]);
        
        subplot(4,4,9)
        [n,bin] = histc(spikes.times{i},tracking.timestamps);
        plot(tracking.position.x,tracking.position.y,'LineWidth',1,'Color',[0.5 0.5 0.5])
        hold on
        view(0,-90)
        plot(tracking.position.x(bin(bin>0)),tracking.position.y(bin(bin>0)),'.','MarkerEdgeColor','r','MarkerSize',3);
        hold off
        set(gca,'DataAspectRatio',[1 1 1]);  
        title(['Max Fr: ', num2str(firingMaps.stats{i}{1}.peak)])

        subplot(4,4,10)
        % Occupancy
        imagesc(firingMaps.occupancy{i}{1});
        colormap(jet(15)),colorbar, shading flat
%         view(0,-90)
        title('occupancy')
        set(gca,'DataAspectRatio',[1 1 1]);  

        subplot(4,4,11)
        %count
        imagesc(firingMaps.countMaps{i}{1});
        colormap(jet(15)), colorbar, shading flat
        title('count')
%         view(0,-90)
        set(gca,'DataAspectRatio',[1 1 1]);  

        subplot(4,4,12)
        % rateMap
        imagesc(firingMaps.rateMaps{i}{1});
        colormap(jet(15)), colorbar, shading flat
        title('ratemap')
%         view(0,-90)
        set(gca,'DataAspectRatio',[1 1 1]);  
        xlabel(['Skaggs: ', num2str(firingMaps.stats{i}{1}.skaggs.bitsPerSpike)], 'color', color_skaggs )
        ylabel(['r: ', num2str(firingMaps.stats{i}{1}.spatialCorr.sc)], 'color', color_sc)
      
        subplot(4,4,13)
        % Periodic analysis
        r = bz_autoCorrPearson(firingMaps,i);
        imagesc(r), colormap(jet(15)), colorbar, shading flat
        title('Spatial Autocorrelogram')
        set(gca,'DataAspectRatio',[1 1 1]);
        
        

%         saveas(gcf,['placeFields\placeFields_Cell',num2str(i),'.png']);
    end

    % figure
    % set(gcf,'Position',[100 -100 2500 1200])
    % for i=1:length(firingMaps.rateMaps)    
    %     subplot(7,ceil(size(firingMaps.rateMaps,1)/7),i); % autocorrelogram
    %     h = pcolor(firingMaps.rateMaps{i}{1}); 
    %     colormap(jet(15)), colorbar, shading flat  
    %     set(gca,'DataAspectRatio',[1 1 1]);
    % end
    saveas(gcf,'SummaryFigures\placeFields.png');
else
    for i=1:length(firingMaps.rateMaps)
        for j=1:conditions
            figure,
            set(gcf,'Position',get(0,'ScreenSize'))
            subplot(2,4,1)
            plot(spikes.filtWaveform{i})
            title(['Cell:' , num2str(i)])
            
            subplot(2,4,2)
            %ACG
            area(cell_metrics.acg.narrow(:,i),'LineStyle','none')
            title(['Cell type: ', num2str(cell_metrics.putativeCellType{i})])
            

            subplot(2,4,3)
            % ISI
            area(cell_metrics.isi.log10(:,i),'LineStyle','none')
            plot(cell_metrics.isi.log10(:,i))
            title(['Burst Index: ', num2str(spkTrain{j}{i}.BurstIndex)])
            
            subplot(2,4,4)
            % Trilateralization
            plot(cell_metrics.general.chanCoords.x,cell_metrics.general.chanCoords.y,'.k'), hold on
            plot(cell_metrics.trilat_x(i),cell_metrics.trilat_y(i),'ob'), xlabel('x position (µm)'), ylabel('y position (µm)')
            title([ ' Shank:' ,num2str(spikes.shankID(i)),])
            
            subplot(2,4,5)
            [n,bin] = histc(spikes.times{i}(spikes.times{i} > tracking.events.subSessions(j,1) & spikes.times{i} < tracking.events.subSessions(j,2)),tracking.timestamps(tracking.events.subSessionsMask == j));
            plot(tracking.position.x(tracking.events.subSessionsMask == j),tracking.position.y(tracking.events.subSessionsMask == j),'Color',[0.5 0.5 0.5])
            hold on
            view(0,-90)
            xToPlot = tracking.position.x(tracking.events.subSessionsMask == j);
            yToPlot = tracking.position.y(tracking.events.subSessionsMask == j);
            plot(xToPlot(bin(bin>0)),yToPlot(bin(bin>0)),'.','MarkerEdgeColor','r','MarkerSize',15);
            hold off
            set(gca,'DataAspectRatio',[1 1 1]);
            title(['Mean Firing Rate: ' , num2str(firingMaps.stats{i}{j}.meanFr)])
            
            subplot(2,4,6)
            % Occupancy
            if strcmpi(firingMaps.params.analysis,'tint')
                imagesc(firingMaps.occupancy{i}{j});
                colormap(jet(15)),colorbar, shading flat
                title('occupancy')
                set(gca,'DataAspectRatio',[1 1 1]);
            elseif strcmpi(firingMaps.params.analysis,'buzcode')
                imagesc(firingMaps.occupancy{i}{j});
                colormap(jet(15)),colorbar, shading flat
                view(0,-90)
                title('occupancy')
                set(gca,'DataAspectRatio',[1 1 1]);
            end
            
            
            subplot(2,4,7)
            %count
            if strcmpi(firingMaps.params.analysis,'tint')
                imagesc(firingMaps.countMaps{i}{j});
                colormap(jet(15)), colorbar, shading flat
                title('count')
                set(gca,'DataAspectRatio',[1 1 1]); 
            elseif strcpi(firingMaps.params.analysis,'buzcode')
                imagesc(firingMaps.countMaps{i}{j});
                colormap(jet(15)), colorbar, shading flat
                title('count')
                view(0,-90)
                set(gca,'DataAspectRatio',[1 1 1]);
            end
            
            subplot(2,4,8)
            % rateMap
            if strcmpi(firingMaps.params.analysis,'tint')
                imagesc(firingMaps.rateMaps{i}{j});
                colormap(jet(15)), colorbar, shading flat
                title('ratemap')
                set(gca,'DataAspectRatio',[1 1 1]); 
            elseif strcmpi(firingMaps.params.analysis,'buzcode')
                imagesc(firingMaps.rateMaps{i}{j});
                colormap(jet(15)), colorbar, shading flat
                title('ratemap')
                view(0,-90)
                set(gca,'DataAspectRatio',[1 1 1]); 
            end
            
            saveas(gcf,['placeFields\placeFields_Cell',num2str(i),'folder',tracking.folders{j},'.png']);
        
        
        end
    end
    
    
    
    
end

end

