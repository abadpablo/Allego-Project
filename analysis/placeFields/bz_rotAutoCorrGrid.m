function [grid] =bz_rotAutoCorrGrid(map,varargin)

% Computes Grid Index
% Adapted from rotAutoCorr_Grid2 by Pablo Abad 23/07 to resemble buzCode
%
% USAGE
%   [grid] = bz_rotAutoCorrGrid(map,<options>);
%
% INPUT 
%   map             map obtained using <a href="matlab:help Map">Map</a>
%                    
%   options         optional list of property-value pairs(see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'    By default pwd
%     'smoothing'   smoothed or not smoothed (map) but this has to be organize no to
%                   over smooth the firing maps. (default true)
%    =========================================================================
%
%
%   OUTPUT
%       
%       grid.grid        west
%       grid.gridness        east
%       grid.squareIndex    south
%       grid.r       north

% Pablo Abad 07/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default and Params

p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'smoothing',true,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
smoothing = p.Results.smoothing;

if smoothing
    Map = map.z;
else
    Map = map.zUnSmooth;
end


nan_values = isnan(Map);
[a1,a2] = find(nan_values == 1);
Map(a1,a2) = 0;


[nrClassX nrClassY] = size(Map);


if  prod(nrClassX,nrClassY) > 1000
    fprintf('Please note that grid analysis with classes above 100 can be very timeconsuming')
end

% tetnr = 1;
% cellnr = 1;

MapN = (max(max(Map))*Map)/max(max(Map));
%% Prepare empty Matrices
[amountRowOrig, amountColOrig] = size(Map);
amountColCorr = amountColOrig;
amountRowCorr = amountRowOrig;

%r is the Matrix including the autocorrelation
r = zeros(amountRowCorr*2,amountColCorr*2);

count = 0;
%% Autocorrelation

Matrix_anal = 1;

if Matrix_anal==1
% wh = waitbar(0, (['Calculating Autocorrelogram of tetnr ' num2str(tetnr) ' cellnr ' num2str(cellnr)]));

warning off
for xx=1:1:amountRowCorr
%     waitbar(xx/amountRowCorr,wh);
    for yy=1:1:amountColCorr
        count=count+1;

        AA = Map(1:xx,1:yy);
        BB = Map(amountRowOrig-xx+1:amountRowOrig,amountColOrig-yy+1:amountColOrig);
        AAx= reshape(AA,1,numel(AA));
        BBy= reshape(BB,1,numel(BB));
        temp = corrcoef(AAx,BBy);
        
        if prod(size(temp))>= 4 
            r(xx,yy) = temp(2,1);
        else
            r(xx,yy) = 0;
        end
    
        in1_A=xx:amountRowOrig;
        in1_B=yy:amountColOrig;
        AA = Map(xx:amountRowOrig,yy:amountColOrig);
        BB = Map(1:amountRowOrig-xx+1,1:amountColOrig-yy+1);
        AAx= reshape(AA,1,numel(AA));
        BBy= reshape(BB,1,numel(BB));
        temp = corrcoef(AAx,BBy);

        if prod(size(temp))>=4 
            r(amountRowOrig+xx,amountColOrig+yy) = temp(2,1);
        else
            r(amountRowOrig+xx,amountColOrig+yy) = 0;
        end
    
        AA =Map(1:xx,amountColOrig-yy+1:amountColOrig);
        BB = Map(amountRowOrig-xx+1:amountRowOrig,1:yy);
        AAx= reshape(AA,1,numel(AA));
        BBy= reshape(BB,1,numel(BB));
        temp = corrcoef(AAx,BBy);
    
        if prod(size(temp))>=4 
            r(xx,2*amountColOrig+1-yy) = temp(2,1);
        else
            r(xx,2*amountColOrig+1-yy) = 0;
        end
      
    
%         AA = Matrix(xx:amountRowOrig,1:amountRowOrig-yy+1);
        AA = Map(xx:amountRowOrig,1:amountColOrig-yy+1);
        BB = Map(1:amountRowOrig-xx+1,yy:amountColOrig);
        AAx= reshape(AA,1,numel(AA));
        BBy= reshape(BB,1,numel(BB));
        temp = corrcoef(AAx,BBy);
    
        if prod(size(temp))>=4 
            r(xx+amountRowOrig,amountColOrig+1-yy) = temp(2,1);
        else
            r(xx+amountRowOrig,amountColOrig+1-yy) = 0;
        end

        r(find(isnan(r)))=0; %eliminate NaN
%   xx
%   yy
% PlotGrid=1;
% if PlotGrid==1
% if count==1
%     figure;
%     imagesc(r);hold
% elseif count>1
%     imagesc(r)
% end
% end
%  pause
    end
    end
end



%% Plotting and saving

% arenaSize=[0.5 0.5];
classes = size(Matrix); 
classes = classes(1);

if RotCal==1
figure;
[correlation1,angle1,gridness1,squareIndex1,circumference1] = rotAutoCorrJ(r,arenaSize,classes);
% [correlation1,angle1,gridness1,squareIndex1,circumference1] = rotAutoCorr(r);

AutoCoSource=[ExperCell,'_',num2str(kl),'_Source'];
makedir=0
if makedir==1
% mkdir(['C:\Users\Jorge\Documents\MATLAB\ERB4\QualityClustering\' ExperCell  '\GridAnalysis\' ]);
mkdir(['G:\Proyects\Subiculum\PeriodicCells\'  ExperCell ,'\Spatial\']);

saveas(gcf,['G:\Proyects\Subiculum\PeriodicCells\'  ExperCell ,'\Spatial\' AutoCoSource,'.eps']);

saveas(gcf,['G:\Proyects\Subiculum\PeriodicCells\'  ExperCell ,'\Spatial\' AutoCoSource,'.tiff']);
% mkdir(['C:\Users\Jorge\Documents\MATLAB\ERB4\QualityClustering\' ExperCell  '\GridAnalysis\' ]);

Max1=[ExperCell,'_',num2str(kl),'_Max'];
hgsave(['G:\Proyects\Subiculum\PeriodicCells\'  ExperCell ,'\Spatial\' Max1]);
saveas(gcf,['G:\Proyects\Subiculum\PeriodicCells\'  ExperCell ,'\Spatial\' Max1,'.tiff']);
saveas(gcf,['G:\Proyects\Subiculum\PeriodicCells\'  ExperCell ,'\Spatial\' Max1,'.eps']);

% mkdir(['C:\Users\Jorge\Documents\MATLAB\ERB4\QualityClustering\' ExperCell ]);
% saveas(gcf,['C:\Users\Jorge\Documents\MATLAB\ERB4\QualityClustering\' ExperCell '\' FigName '.tiff']);

%% Calc grid geometry:
end
threshold=0.05;

regionalMax1 = findRegionalMaxima3DJ( r, arenaSize,classes, circumference1, threshold );

[grid.Orientation1 grid.spacingOfGrid1 grid.gridInfoTable1 grid.sumAngle1 ...
grid.OrientationS1 grid.spacingOfGridS1 grid.gridInfoTableS1 grid.sumAngleS1] ...
           = gridMeasureJ( regionalMax1,arenaSize,size(Matrix) );

fivem =[0 0 1;0 0.6 1; 0 1 0;1 1 0;1 0 0] % generates colormap for firing fields in 20% interval similar to tint.
%   figure;imagesc(r);colormap(jet(15));shading flat% colormap(105);%colormap(fivem)
  figure;pcolor(r);shading flat;colormap(jet(15))
else
    grid=0; gridness1=0; squareIndex1=0;
end
GridPlot=[ExperCell,'_',num2str(kl),'_SAC'];

  p.NotVisited=1;  %visualize unvisited places
 p.EdgeColor=1;
 p.ColorMap='tint';   %load colormap tint
  
%  t=pcolor(map1);
% 
% if p.EdgeColor
%     set(t,'EdgeColor','flat')
% else
%     set(t,'EdgeColor','none')
% end
% %shading('interp')
% %set(gca,'CLim',[0 max(max(tmp_matrix))]);
% set(gca,'DataAspectRatio',[1 1 1]);  
% axis tight
% switch p.ColorMap;
%     case 'tint'
%         colormap(gettintcolormap)
%     case 'jet'
%         colormap('jet');
% end
mice=0
if mice==1
  colorbar

hgsave(['G:\Proyects\Subiculum\PeriodicCells\'  ExperCell ,'\Spatial\' GridPlot]);
saveas(gcf,['G:\Proyects\Subiculum\PeriodicCells\'  ExperCell ,'\Spatial\' GridPlot,'.tiff']);
end
close 


