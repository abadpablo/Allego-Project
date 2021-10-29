 function [shuffling,shuffling_maps]=bz_SpikeShuffling(positions,spikes,varargin)
 
% Computes randomization of the map
% Adapted from SpikeShuffling_spk2 by Pablo Abad 23/07 to resemble buzCode
%
% USAGE
%   [shuffling] = bz_SpikeShuffling(map,<options>);
%
% INPUT 
%   positions       [t x y]
%
%   map             map obtained using <a href="matlab:help Map">Map</a>
%   
%   spikes          spikes structure
%
%   options         optional list of property-value pairs(see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'    By default pwd
%     'smoothing'   smoothed or not smoothed (map) but this has to be organize no to
%                   over smooth the firing maps. (default true)
%     'smooth'      smoothing size in bins ( 0 = no smoothing, default = 2)
%     'nBins'       
%    =========================================================================
%
%
%   OUTPUT
%       
%       shuffling


% Pablo Abad 07/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Default and Params

p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'numRand',1000,@isnumeric);
addParameter(p,'sr',30000,@isnumeric);
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'type','lll',@isstr);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'nBins',20,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'showFig',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
numRand = p.Results.numRand;
sr = p.Results.sr;
smooth = p.Results.smooth;
type = p.Results.type;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
minTime = p.Results.minTime;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
showFig = p.Results.showFig;

t = positions(:,1);
minA = 20; % Min value 20 secs ( as in Kupric)
maxA = max(t)- min(t);
ValRand = abs(ceil((maxA-minA).*rand(numRand,1)-20 ));%
numSpikes = length(spikes);
ValRand(1) = 20;

%% ORIGINAL FIRING MAP
map = Map(positions,spikes,'smooth',smooth,'minTime',minTime,...
            'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
        
stats_original = MapStats(map,spikes,'verbose','off');
%%

spikes_aux = spikes(spikes > t(1) & spikes < t(end));


for i=1:numRand
    spikes_rand = spikes_aux+ValRand(i);
    indexunder = find(spikes_rand < t(end));
    indexover = find(spikes_rand > t(end));
    overlim = isempty(indexover);
    if overlim
        spikes_rand3 = spikes_rand;
    else
        spikes_rand1 = sort(spikes_rand(indexunder));
        spikes_rand2 = abs(spikes_rand(indexover)-t(end) + t(1));
        spikes_rand3 = sort([spikes_rand1; spikes_rand2]);
    end
    
    shuffling{i} = Map(positions,spikes_rand3,'smooth',smooth,'minTime',minTime,...
                        'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                    
    stats{i} = MapStats(shuffling{i},spikes_rand3,'verbose','off');
    
    % Spatial Stability
    [r p] = corrcoef(map.z,shuffling{i}.z);
    stability{i} = r(1,2);
    stats{i}.stability = stability{i};
    
    
    if showFig
        figure,
        imagesc(shuffling{i}.z), colormap(jet)
    end
end


%% restructure into cell info data type

shuff = shuffling;
clear shuffling;

shuffling_maps.maps = shuff;
shuffling_maps.stats = stats;


%% Statistics

specificity_shuffling = zeros(1,length(stats));
sc_shuffling = zeros(1,length(stats));
bitsPerSpike_shuffling = zeros(1,length(stats));
stability_shuffling = zeros(1,length(stats));
west_shuffling = zeros(1,length(stats));
east_shuffling = zeros(1,length(stats));
north_shuffling = zeros(1,length(stats));
south_shuffling = zeros(1,length(stats));


for i=1:length(stats)    
    sc_shuffling(i) = stats{i}.spatialCorr.sc;
    bitsPerSpike_shuffling(i) = stats{i}.skaggs.bitsPerSpike;
    specificity_shuffling(i) = stats{i}.specificity;
%     stability_shuffling(i) = stats{i}.stability;
    west_shuffling(i) = stats{i}.borderIndex.west;
    east_shuffling(i) = stats{i}.borderIndex.east;
    north_shuffling(i) = stats{i}.borderIndex.north;
    south_shuffling(i) = stats{i}.borderIndex.south;    
end


sc_shuffling_mean = mean(sc_shuffling);
sc_shuffling_std = std(sc_shuffling);

bitsPerSpike_shuffling_mean = mean(bitsPerSpike_shuffling);
bitsPerSpike_shuffling_std = std(bitsPerSpike_shuffling);

specificity_shuffling_mean = mean(specificity_shuffling);
specificity_shuffling_std = std(specificity_shuffling);

% stability_shuffling_mean = mean(stability_shuffling);
% stability_shuffling_std = std(stability_shuffling);

west_shuffling_mean = mean(west_shuffling);
west_shuffling_std = std(west_shuffling);

east_shuffling_mean = mean(east_shuffling);
east_shuffling_std = std(east_shuffling);

north_shuffling_mean = mean(north_shuffling);
north_shuffling_std = std(north_shuffling);

south_shuffling_mean = mean(south_shuffling);
south_shuffling_std = std(south_shuffling);

   
    
% Finding specific percentil 99% and 95% for cutting value
sc_shuffling99 = prctile(sc_shuffling,99);
sc_shuffling95 = prctile(sc_shuffling,95);
   
bitsPerSpike_shuffling99 = prctile(bitsPerSpike_shuffling,99);
bitsPerSpike_shuffling95 = prctile(bitsPerSpike_shuffling,95);

specificity_shuffling99 = prctile(specificity_shuffling,99);
specificity_shuffling95 = prctile(specificity_shuffling,95);

% stability_shuffling99 = prctile(stability_shuffling,99);
% stability_shuffling95 = prctile(stability_shuffling,95);

west_shuffling99 = prctile(west_shuffling,99);
west_shuffling95 = prctile(west_shuffling,95);

east_shuffling99 = prctile(east_shuffling,99);
east_shuffling95 = prctile(east_shuffling,95);

north_shuffling99 = prctile(north_shuffling,99);
north_shuffling95 = prctile(north_shuffling,95);

south_shuffling99 = prctile(south_shuffling,99);
south_shuffling95 = prctile(south_shuffling,95);


shuffling.mean.sc = sc_shuffling_mean;
shuffling.mean.bitsPerSpike = bitsPerSpike_shuffling_mean;
shuffling.mean.specificity = specificity_shuffling_mean;
% shuffling.mean.stability = stability_shuffling_mean;
shuffling.mean.west = west_shuffling_mean;
shuffling.mean.east = east_shuffling_mean;
shuffling.mean.north = north_shuffling_mean;
shuffling.mean.south = south_shuffling_mean;

shuffling.std.sc = sc_shuffling_std;
shuffling.std.bitsPerSpike = bitsPerSpike_shuffling_std;
shuffling.std.specificity = specificity_shuffling_std;
% shuffling.std.stability = stability_shuffling_std;
shuffling.std.west = west_shuffling_std;
shuffling.std.east = east_shuffling_std;
shuffling.std.north = north_shuffling_std;
shuffling.std.south = south_shuffling_std;


shuffling.prctile99.sc = sc_shuffling99;
shuffling.prctile99.bitsPerSpike = bitsPerSpike_shuffling99;
shuffling.prctile99.specificity = specificity_shuffling99;
% shuffling.prctile99.stability = stability_shuffling99;
shuffling.prctile99.west = west_shuffling99;
shuffling.prctile99.east = east_shuffling99;
shuffling.prctile99.north = north_shuffling99;
shuffling.prctile99.south = south_shuffling99;


shuffling.prctile95.sc = sc_shuffling95;
shuffling.prctile95.bitsPerSpike = bitsPerSpike_shuffling95;
shuffling.prctile95.specificity = specificity_shuffling95;
% shuffling.prctile95.stability = stability_shuffling95;
shuffling.prctile95.west = west_shuffling95;
shuffling.prctile95.east = east_shuffling95;
shuffling.prctile95.north = north_shuffling95;
shuffling.prctile95.south = south_shuffling95;   
    
   
end


