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


mkdir('placeFields')
conditions = size(tracking.events.subSessions,1);


if conditions == 1

    figure
    set(gcf,'Position',[100 -100 2500 1200])
    for i=1:length(firingMaps.rateMaps)
        subplot(2,4,1)
        plot(spikes.filtWaveform{i})
        title(['Cell:' , num2str(i), ' Shank:' ,num2str(spikes.shankID(i))])

        subplot(2,4,2)
        %ACG

        subplot(2,4,3)
        % ISI

        subplot(2,4,4)
        % Trilateralization

        subplot(2,4,5)
        [n,bin] = histc(spikes.times{i},tracking.timestamps);
        plot(tracking.position.x,tracking.position.y,'Color',[0.5 0.5 0.5])
        hold on
        view(0,-90)
        plot(tracking.position.x(bin(bin>0)),tracking.position.y(bin(bin>0)),'.','MarkerEdgeColor','r','MarkerSize',15);
        hold off
        set(gca,'DataAspectRatio',[1 1 1]);  

        subplot(2,4,6)
        % Occupancy
        imagesc(firingMaps.occupancy{i}{1});
        colormap(jet(15)),colorbar, shading flat
        view(0,-90)
        title('occupancy')
        set(gca,'DataAspectRatio',[1 1 1]);  
    %     h = pcolor(firingMaps.occupancy{i}{1});
    %     colormap(jet(15)), colorbar, shading flat
    %     title('occupancy')
    %     view(0,-90)

        subplot(2,4,7)
        %count
        imagesc(firingMaps.countMaps{i}{1});
        colormap(jet(15)), colorbar, shading flat
        title('count')
        view(0,-90)
        set(gca,'DataAspectRatio',[1 1 1]);  

    %     h = pcolor(firingMaps.countMaps{i}{1});
    %     colormap(jet(15)), colorbar, shading flat  
    %     title('count')

        subplot(2,4,8)
        % rateMap
        imagesc(firingMaps.rateMaps{i}{1});
        colormap(jet(15)), colorbar, shading flat
        title('ratemap')
        view(0,-90)
        set(gca,'DataAspectRatio',[1 1 1]);  
    %     h = pcolor(firingMaps.rateMaps{i}{1}); 
    %     colormap(jet(15)), colorbar, shading flat  
    %     title('rateMaps')

        saveas(gcf,['placeFields\placeFields_Cell',num2str(i),'.png']);
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
            title(['Cell:' , num2str(i), '  MeanFr: ' , num2str(firingMaps.stats{i}{j}.meanFr)])
            
            subplot(2,4,2)
            %ACG
            area(cell_metrics.acg.narrow(:,i),'LineStyle','none')
            title(['Cell type: ' ,num2str(cell_metrics.putativeCellType{i})])

            subplot(2,4,3)
            % ISI
            area(cell_metrics.isi.log10(:,i),'LineStyle','none')
            plot(cell_metrics.isi.log10(:,i))
            
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

