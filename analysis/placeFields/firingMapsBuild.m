function [map] = firingMapsBuild(position,times,varargin)
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
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'nBins',20,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'type','lll',@isstr);

parse(p,varargin{:});
smooth = p.Results.smooth;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
maxDistance = p.Results.maxDistance;
mode = p.Results.mode;
type = p.Results.type;


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


% Make sure x and y are normalized
if max(x) > 1 || min(x) < 0
    x = ZeroToOne(x);
end
if ~isempty(y)
    if max(y) > 1 || min(y) < 0
        y = ZeroToOne(y);
    end
end

% Number of bins for x and y
nBinsX = nBins(1);
if length(nBins) == 1
    nBinsY = nBinsX;
    nBins(2) = nBins;
else
    nBinsY = nBins(2);
end

% Bin x and y
x = Bin(x,[0 1],nBinsX);
if ~isempty(y)
    y = Bin(y,[0 1],nBinsY);
end

% Duration for each (X,Y) sample (clipped to maxGap)
dt = diff(t);
dt(end+1)=dt(end);
dt(dt>maxGap) = maxGap;

if pointProcess,
	% Count occurrences for each (x,y) timestamp
	n = CountInIntervals(times,[t t+dt]);
else
	% Interpolate z at (x,y) timestamps
	[z,discarded] = Interpolate(z,t,'maxGap',maxGap);
	if isempty(z), return; end
	if strcmp(type(end),'c'),
		range = isradians(z(:,2));
		z(:,2) = exp(j*z(:,2));
	end
	n = 1;
end


% Computations
if isempty(y)
	% 1D (only x)
	map.x = linspace(0,1,nBinsX);
	map.count = Accumulate(x,n,nBinsX);
	map.time = Accumulate(x,dt,nBinsX);
	valid = map.time > minTime;
	map.count = Smooth(Interpolate1(map.x,map.count,valid,mode,maxDistance),smooth,'type',type(1))';
	map.time = Smooth(Interpolate1(map.x,map.time,valid,mode,maxDistance),smooth,'type',type(1))';
	if pointProcess,
		map.z = map.count./(map.time+eps);
	else
		map.z = Accumulate(x,z(:,2),nBinsX);
		map.z = Smooth(Interpolate1(map.x,map.z,valid,mode,maxDistance),smooth,'type',type(1))';
		map.z = map.z./(map.count+eps);
	end
else
	% 2D (x and y)
	map.x = linspace(0,1,nBinsX);
	map.y = linspace(0,1,nBinsY);
	map.count = Accumulate([x y],n,nBins);
	map.time = Accumulate([x y],dt,nBins);
    map.count = rot90(map.count);
    map.time = rot90(map.time);
    
	valid = map.time > minTime;
    count = map.count; %unsmooth
    time = map.time; % unsmooth
    
    [map.time,map.count,map.z] = smooth_field_plot(time,count,7,'boxcar');
    
    map.countUnSmooth = count;
    map.timeUnSmooth = time;
    
	if pointProcess,
        map.zUnSmooth = map.countUnSmooth./(map.timeUnSmooth);
	else
		map.z = Accumulate([x y],z(:,2),nBins);
		map.z = Smooth(Interpolate2(map.x,map.y,map.z,valid,mode,maxDistance),smooth,'type',type(1:2)).';
		map.z = map.z./(map.count+eps);
	end
end



end

