function plot_placeFields(firingMaps,spikes,tracking,varargin)

p = inputParser;
addParameter(p,'basepath',pwd,@isstr);

parse(p,varargin{:})

basepath = p.Results.basepath;

mkdir('placeFields')

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

end

