function [firingMaps] = bz_firingMapAvg(positions,spikes,varargin)

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
%   behavior  - buzcode format behavior struct - 
%   <options>      optional list of property-value pairs (see table below)
% ===================================================
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
%                .cmBin                 cm/bins ratio
%
% Antonio FR, 10/2019

%% parse inputs
p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
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
addParameter(p,'speedFilter',false,@islogical);
addParameter(p,'cmBin',2.5,@isnumeric);
addParameter(p,'periodicAnalysis',false,@islogical);
addParameter(p,'numRand',10,@isnumeric);
addParameter(p,'spikeShuffling',true,@islogical);

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
speedFilter = p.Results.speedFilter;
cmBin = p.Results.cmBin;
periodicAnalysis = p.Results.periodicAnalysis;
numRand = p.Results.numRand;
basepath = p.Results.basepath;
spikeShuffling = p.Results.spikeShuffling;

%% In case firingMapsAvg already exists
if ~isempty(dir([basepath filesep '*firingMapsAvg.cellinfo.mat']))
    disp('Firing Maps already detected! Loading file.');
    file = dir([basepath filesep '*firingMapsAvg.cellinfo.mat']);
    load(file.name);
    return
end
% We need to obtain the number of bins
tracking = getSessionTracking();
if isfield(tracking,'apparatus')
    if size(tracking.apparatus,2) == 1
        numApparatus = 1;
        apparatus = tracking.apparatus{1};
        apparatus_name = apparatus.name;
        xLim = apparatus.boundingbox.xmax - apparatus.boundingbox.xmin;
        yLim = apparatus.boundingbox.ymax - apparatus.boundingbox.ymin;  
    else
        numApparatus = size(tracking.apparatus,2);
        for i=1:numApparatus
            apparatus{i} = tracking.apparatus{i};
            apparatus_name{i} = apparatus{i}.name;
            xLim{i} = apparatus{i}.boundingbox.xmax - apparatus{i}.boundingbox.xmin;
            yLim{i} = apparatus{i}.boundingbox.ymax - apparatus{i}.boundingbox.ymin;
        end
    end
end

if numApparatus == 1
    nBins{1} = round(xLim/cmBin);
else
    clear nBins
    for i = 1:numApparatus
        nBins{i} = round(xLim{i}/cmBin);
    end
end

if isstruct(positions)
    positions = positions.maps;
end

% number of conditions
  if iscell(positions)
     conditions = length(positions); 
  elseif isvector(positions)
     conditions = 1;
  end
  %%% TODO: conditions label
%   if conditions == 3
%       nBins(1) = numberBins{1};
%       nBins(2) = numberBins{1};
%       nBins(3) = numberBins{2};
%   end
  

for i=1:length(positions)   
    positions{i} = positions{i}(tracking.events.subSessionsMask == i,:);  
end
%% Calculate
% Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    post = positions{iCond}(:,1);
    % - 1D 
    if size(positions{iCond},2)==2
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
    elseif size(positions{iCond},2)==3
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
    else
        warning('This is not a linear nor a 2D space!');
    end
    % Absolute speed
    v = sqrt(vx.^2+vy.^2);    
    if length(v) < length(posx)
        dif = length(posx) - length(v);
        s = v;
        sx = vx;
        sy = vy;
        clear v
        clear vx
        clear vy
        
        v = zeros(1,length(posx));
        vx = zeros(1,length(posx));
        vy = zeros(1,length(posx));
        
        v = [nan(dif,1); s];
        vx = [nan(dif,1); sx];
        vy = [nan(dif,1); sy];        
    end
    tracking.speed.v = v/100; % cm/s
    tracking.speed.vx = vx/100; % cm/s
    tracking.speed.vy = vy/100; % cm/s
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); sessionInfo.rates.lfp = 1250;  
    basepath = sessionInfo.session.path;
    save([basepath filesep sessionInfo.FileName '.Tracking.Behavior.mat'],'tracking');
    % Compute timestamps where speed is under threshold
    if speedFilter == 1
        positions{iCond}(tracking.speed.v<speedThresh,:) = [];
    end
end


%% get firign rate maps & map stats
for unit = 1:length(spikes.times)
    for c = 1:conditions
            map{unit}{c} = Map(positions{c},spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                'nBins',nBins{c},'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
            stats{unit}{c} = MapStats(map{unit}{c},spikes.times{unit},'nBins',nBins{c},'verbose','on');
            
            % shuffling of spikes times
            if spikeShuffling
                shuffling{unit}{c} = bz_SpikeShuffling(positions{c},spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                    'nBins',nBins{c},'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance,'numRand',numRand);
            else
                shuffling = [];
            end
            
            % Periodic Firing
            if periodicAnalysis
                periodic{unit}{c} = bz_PeriodicPower(map{unit}{c},shuffling{unit}{c},'random',true,'nBins',nBins{c});
            else
                periodic{unit}{c} = [];
            end        
    end
end
% cmBin = (max(positions{1}(:,2))-min(positions{1}(:,2)))/nBins;
%%% TODO: pass rest of inputs to Map

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
try firingMaps.sessionName = spikes.sessionName;
catch
    firingMaps.sessionName = spikes.basename;
end
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
firingMaps.params.cmBin = cmBin;
firingMaps.params.numRand = numRand;

for unit = 1:length(spikes.times)
    for c = 1:conditions
    firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
    firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
    firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
    end
end

% Save stats
firingMaps.stats = stats;

% Save shuffling
firingMaps.shuffling = shuffling;

% Save periodic
if periodicAnalysis
    firingMaps.periodic = periodic;
end


for unit = 1:length(spikes.times)
    for c = 1:conditions
        firingMaps.rateMapsUnSmooth{unit,1}{c} =map{unit}{c}.zUnSmooth;
        firingMaps.countMapsUnSmooth{unit,1}{c} = map{unit}{c}.countUnSmooth;
        firingMaps.occupancyUnSmooth{unit,1}{c} = map{unit}{c}.timeUnSmooth;
    end
end


if saveMat
   save([firingMaps.sessionName '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 
end


end
