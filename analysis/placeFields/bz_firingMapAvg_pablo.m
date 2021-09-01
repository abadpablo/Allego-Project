function [firingMaps] = bz_firingMapAvg_pablo(positions,spikes,varargin)

% USAGE
% [firingMaps] = bz_firingMapAvg(positions,spikes,varargin)
% Calculates averaged firing map for a set of linear postions 
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%      or
%   behavior  - buzcode format behavior struct - NOT YET IMPLEMENTED
%   <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'			smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'			number of bins (default = 50)
%     'speedThresh'		speed threshold to compute firing rate
%     'minTime'			minimum time spent in each bin (in s, default = 0)
%     'mode'			interpolate' to interpolate missing points (< minTime),
%                   	or 'discard' to discard them (default)
%     'maxDistance'		maximal distance for interpolation (default = 5)
%     'maxGap'			z values recorded during time gaps between successive (x,y)
%                   	samples exceeding this threshold (e.g. undetects) will not
%                	    be interpolated; also, such long gaps in (x,y) sampling
%                 	    will be clipped to 'maxGap' to compute the occupancy map
%                 	    (default = 0.100 s)
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'saveMat'   		- logical (default: false) that saves firingMaps file
%     'CellInspector'  	- logical (default: false) that creates an otuput
%                   	compatible with CellInspector

%
%
% OUTPUT
%
%   firingMaps - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%
% Antonio FR, 10/2019

%% parse inputs
p=inputParser;
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'speedThresh',0.1,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'CellInspector',false,@islogical);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'bodypart','body',@isstr); % Which tracking variable you want to perform the place field on ( 'body', 'head','tail')
addParameter(p,'paradigm','Open Field',@isstr);


parse(p,varargin{:});
smooth = p.Results.smooth;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
saveMat = p.Results.saveMat;
CellInspector = p.Results.CellInspector;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
order = p.Results.orderKalmanVel;
bodypart = p.Results.bodypart;
paradigm = p.Results.paradigm;

  %%% TODO: conditions label
  apparatus = [];
  
  tracking = getSessionTracking();
  if isfield(tracking,'apparatus')
    apparatus = tracking.apparatus{1};
    xLim = apparatus.boundingbox.xmax - apparatus.boundingbox.xmin;
    yLim = apparatus.boundingbox.ymax - apparatus.boundingbox.ymin;
  end
  if isfield(tracking,'pixelsmetre')
    pixels_metre = tracking.pixelsmetre{1};
  end
  
  % Size of bins 2.5 pixels / 1 cm
  pixelsPerCm = 2.5;
  if strcmp(apparatus.name,'Open Field')
          nBins = xLim / pixelsPerCm;
          nBins = round(nBins);
  end
  
  % number of conditions
  if iscell(positions)
     conditions = length(positions); 
  elseif isstruct(positions)
      if strcmp(paradigm,'Open Field')
        conditions = size(positions,2);
        positions = positions;
      elseif strcmp(paradigm, 'Linear Maze')
          positions = positions.maps;
      end
  else
     conditions = 1;
     temp{1} = positions;
     positions = temp{1};
  end
 
%% Calculate
% Erase positions below speed threshold

for iCond=1:size(positions,2)
    post = positions(iCond).timestamps;
    if positions(iCond).description{1} == 'linearMaze' 
       posx = positions(iCond).position.lin;
       [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
       % Absolute speed
       v = sqrt(vx.^2+vy.^2);
       % Compute timestamps where speed is under threshold
    %   positions(iCond)(v<speedThresh,:) = [];
        positions(iCond).timestamps(v<speedThresh,:) = [];
        positions(iCond).position.lin(v<speedThresh,:) = [];
                 
            elseif strcmp(positions(iCond).description{1},'Open Field')
                if strcmp(bodypart,'body')
                    posx = positions(iCond).position.x;
                    posy = positions(iCond).position.y;
                    [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
                    % Abosulte speed
                    v = sqrt(vx.^2+vy.^2);
                    % Compute timestamps where speed is under theshold
%                     positions(iCond).timestamps(v<speedThresh,:) = [];
%                     positions(iCond).position.x(v<speedThresh,:) = [];
%                     positions(iCond).position.y(v<speedThresh,:) = [];
                    
                elseif strcmp(bodypart,'head')
                    posx = positions(iCond).headposition.x;
                    posy = positions(iCond).headposition.y;
                    [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
                    % Abosulte speed
                    v = sqrt(vx.^2+vy.^2);
                    % Compute timestamps where speed is under theshold
%                     positions(iCond).timestamps(v<speedThresh,:) = [];
%                     positions(iCond).headposition.x(v<speedThresh,:) = [];
%                     positions(iCond).headposition.y(v<speedThresh,:) = [];
                    
                elseif strcmp(bodypart,'tail')
                    posx = positions(iCond).tailposition.x;
                    posy = positions(iCond).tailposition.y;
                    [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
                    % Abosulte speed
                    v = sqrt(vx.^2+vy.^2);
                    % Compute timestamps where speed is under theshold
%                     positions(iCond).timestamps(v<speedThresh,:) = [];
%                     positions(iCond).tailposition.x(v<speedThresh,:) = [];
%                     positions(iCond).tailposition.y(v<speedThresh,:) = [];
                    
                end
            else
                warning('This is not a linear nor a 2D space!');   
            end
           

            
        end
end


%% get firign rate maps
if iscell(positions)
    for unit = 1:length(spikes.times)
        for c = 1:conditions
            map{unit}{c} = Map(positions{c},spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
        end
    end
elseif isstruct(positions)
        for unit=1:length(spikes.times)
            for c=1:conditions
                if strcmp(positions(c).description{1},'linearMaze')
                    map{unit}{c} = Map([positions(c).timestamps positions(c).position.lin],spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                    'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                elseif strcmp(positions(c).description{1},'Open Field')
                    if strcmp(bodypart,'body')
                        map{unit}{c} = Map([positions(c).timestamps positions(c).position.x positions(c).position.y],spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                        'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                    elseif strcmp(bodypart,'head')
                        map{unit}{c} = Map([positions(c).timestamps positions(c).headposition.x positions(c).headposition.y],spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                        'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                    elseif strcmp(bodypart,'tail')
                        map{unit}{c} = Map([positions(c).timestamps positions(c).tailposition.x positions(c).tailposition.y],spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                        'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                    end
                end
            end
        end
end

%%% TODO: pass rest of inputs to Map

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
% firingMaps.sessionName = spikes.sessionName;
firingMaps.sessionName = spikes.basename;

try
firingMaps.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end

firingMaps.params.smooth = smooth;
firingMaps.params.minTime = minTime;
firingMaps.params.nBins = nBins;
firingMaps.params.maxGap = maxGap;
firingMaps.params.mode = mode;
firingMaps.params.maxDistance = maxDistance;

firingMaps.bodypart = bodypart;

if iscell(positions)
    for unit = 1:length(spikes.times)
        for c = 1:conditions
        firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
        firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
        firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
        %Get the x bins back in units of meters...
        firingMaps.xbins{unit,1}{c} = map{unit}{c}.x .* ...
            (max(positions{c}(:,2))-min(positions{c}(:,2))) + min(positions{c}(:,2));
        end
    end
elseif isstruct(positions)
    for unit = 1:length(spikes.times)
        for c = 1:conditions
            if strcmp(positions(c).description{1},'linearMaze')
                firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
                firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
                firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
                %Get the x bins back in units of meters...
                firingMaps.xbins{unit,1}{c} = map{unit}{c}.x .* ...
                    (max(positions(c).position.lin)-min(positions(c).position.lin)) + min(positions(c).position.lin);
            elseif strcmp(positions(c).description{1},'Open Field')
                firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
                firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
                firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
%                 if strcmp(bodypart,'body')
%                     %Get the x bins back in units of meters...
%                     firingMaps.xbins{unit,1}{c} = map{unit}{c}.x .* ...
%                         (max(positions(c).position.x)-min(positions(c).position.x)) + min(positions(c).position.x);
%                     firingMaps.ybins{unit,1}{c} = map{unit}{c}.y .*...
%                         (max(positions(c).position.y)-min(positions(c).position.y)) + min(positions(c).position.y);
%                 elseif strcmp(bodypart,'head')
%                     %Get the x bins back in units of meters...
%                     firingMaps.xbins{unit,1}{c} = map{unit}{c}.x .* ...
%                         (max(positions(c).headposition.x)-min(positions(c).headposition.x)) + min(positions(c).headposition.x);
%                     firingMaps.ybins{unit,1}{c} = map{unit}{c}.y .*...
%                         (max(positions(c).headposition.y)-min(positions(c).headposition.y)) + min(positions(c).headposition.y);
%                 elseif strcmp(bodypart,'tail')
%                     %Get the x bins back in units of meters...
%                     firingMaps.xbins{unit,1}{c} = map{unit}{c}.x .* ...
%                         (max(positions(c).tailposition.x)-min(positions(c).tailposition.x)) + min(positions(c).tailposition.x);
%                     firingMaps.ybins{unit,1}{c} = map{unit}{c}.y .*...
%                         (max(positions(c).tailposition.y)-min(positions(c).tailposition.y)) + min(positions(c).tailposition.y);
%                 end
                firingMaps.xbins{unit,1}{c} = map{unit}{c}.x .*...
                    (apparatus.boundingbox.xmax - apparatus.boundingbox.xmin) + apparatus.boundingbox.xmin;
                firingMaps.ybins{unit,1}{c} = map{unit}{c}.y .*...
                    (apparatus.boundingbox.ymax - apparatus.boundingbox.ymin) + apparatus.boundingbox.ymin;
                
        end
    end
    
    
    
end



if saveMat
   save([firingMaps.sessionName '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 
end

end
