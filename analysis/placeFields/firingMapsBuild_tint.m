function [map] = firingMapsBuild_tint(position,times,varargin)
% Map - Map z on (x,y) where x, y and z are time-varying variables (samples).
% Function based on TINT analysis for place cells.
%
%  Compute a continuous map, where one time-varying variable z is represented
%  as a function of one or two time-varying variables x and y. The variable z
%  can either be a point process (typically, a list of spike timestamps) or a
%  continuous measure (e.g. the instantaneous velocity of the animal, the
%  spectral power of an LFP channel in a given frequency band, the coherence
%  between two oscillating LFP channels, etc.) Typical examples of x and y
%  include spatial coordinates and angular directions.
%
%  An occupancy map is also computed.
%
%  USAGE
%
%    map = firingMapsBuild(,,<options>)
%    position       [t1 x y]
%    times          [t2 z]
%    t1             timestamps for x and y
%    x              x values in [0,1]
%    y              optional y values in [0,1]
%    t2             timestamps for z
%    z              optional z values
%    <options>      optional list of property-value pairs (see table below)
%
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'       number of horizontal and vertical bins (default = [50 50])
%     'minTime'     minimum time spent in each bin (in s, default = 0)
%     'mode'        'interpolate' to interpolate missing points (< minTime),
%                   or 'discard' to discard them (default)
%     'maxDistance' maximal distance for interpolation (default = 5)
%     'maxGap'      z values recorded during time gaps between successive (x,y)
%                   samples exceeding this threshold (e.g. undetects) will not
%                   be interpolated; also, such long gaps in (x,y) sampling
%                   will be clipped to 'maxGap' to compute the occupancy map
%                   (default = 0.100 s)
%     'type'        three letters (one for X, one for Y and one for Z) indi-
%                   cating which coordinates are linear ('l') and which are
%                   circular ('c') - for 1D data, only two letters are used
%                   (default 'lll')
%    =========================================================================

%% Parse Inputs

p = inputParser;
addParameter(p,'smooth',7,@isnumeric);
addParameter(p,'nBins',20,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'type','lll',@isstr);
addParameter(p,'boundingbox',[],@isstruct);
addParameter(p,'var2binby','position',@isstr);
addParameter(p,'pixels_metre',[],@isnumeric);
addParameter(p,'tracking',[],@isstruct);

parse(p,varargin{:});
smooth = p.Results.smooth;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
maxDistance = p.Results.maxDistance;
mode = p.Results.mode;
type = p.Results.type;
boundingbox = p.Results.boundingbox;
var2binby = p.Results.var2binby;
pixels_metre = p.Results.pixels_metre;
tracking = p.Results.tracking;

%% Default Output
map.count = [];
map.time = [];
map.z = [];
map.countUnSmooth = [];
map.timeUnSmooth = [];
map.zUnSmooth = [];


if isempty(position) || size(position,1) < 2
    return
end

pointProcess = (isempty(times) | size(times,2) == 1);
t = position(:,1);
x = position(:,2);
if size(position,2) >= 3
    y = position(:,3);
else
    y = [];
end

% Number of bins for x and y
nBinsX = nBins(1);
if length(nBins) == 1
    nBinsY = nBinsX;
    nBins(2) = nBins;
else
    nBinsY = nBins(2);
end

xy = [x y];

% bin the position data
nBins = pixels_metre*2.5 / 100;
[pos_binned_array] = bz_bin_pos_data(var2binby,nBins,xy,1:numel(xy)/2,'boundingbox',boundingbox,'pixels_metre',pixels_metre);

% bin the spike data
[clust_binned_array] = bz_bin_pos_data(var2binby,nBins,xy,times,'boundingbox',boundingbox,'pixels_metre',pixels_metre);

% unsmoothed_rate = (clust_binned_array./pos_binned_array)*tracking.samplingRate(1);
unsmoothed_rate = (clust_binned_array./pos_binned_array);

% smooth the data
[smoothed_pos,smoothed_spikes,smoothed_rate] = smooth_field_plot(pos_binned_array,clust_binned_array,smooth,'boxcar');

map.count = smoothed_spikes;
map.time = smoothed_pos;
map.z = smoothed_rate*tracking.samplingRate(1);

map.countUnSmooth = clust_binned_array;
map.timeUnSmooth = pos_binned_array;
map.zUnSmooth = unsmoothed_rate*tracking.samplingRate(1);;

end

