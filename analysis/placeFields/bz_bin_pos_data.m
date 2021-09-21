function [binned_array,grid_values] = bz_bin_pos_data(var2binby,binsize,posdata,cluster_xy_index,varargin)
% Inputs:
% var2binby = 'position', 'direction', 'speed' or 'pxd' - the variables
%               to be binned by, in order (the values at the centre of each bin to be put
%               in grid_values): units of cm for position, degrees for direction, cm/s for speed
%               - depend on pixels_per_metre in .pos file header *with y increasing upwards*.
% binsize = [8], [8], [8 8] - the sizes of the bins in the each grid, in order matching
%               var2binby. Binsize units for 'position' are camera pixels, degrees for 'direction'
%               Square bins are assumed for 'position', and the range of position is
%               the area tracked by the camera (using window_min_x etc in .pos header).
%               For 'direction' binsize should be a factor of 360.
% posdata - position data in the format of global: TintStructure(index).data{i}.pos
% pos2use - list of samples to be binned (e.g. the position samples corresponding to spikes being fired
%               by a given cell, or the entire list of positions, or those facing North etc etc). Can include
%               repeats of a given position sample.
% Outputs:
% binned_array = requested data binned as required - ij format (y, x) for binning by 'position',
%       column vector for 'direction', 'speed'.
% grid_values = matching shaped array of values corresponding to the centre of each bin
%               e.g. for {'direction' 'position'} you get grid_values(thetai, yi,
%               xi).theta, .x and .y
%
%
%

%% Parse Inputs

p = inputParser;
addParameter(p,'boundingbox',[],@isstruct);
addParameter(p,'pixels_metre',[],@isnumeric);

parse(p,varargin{:});
boundingbox = p.Results.boundingbox;
pixels_metre = p.Results.pixels_metre;


win_max_x = boundingbox.xmax;
win_min_x = boundingbox.xmin;
win_max_y = boundingbox.ymax;
win_min_y = boundingbox.ymin;

extent_x = win_max_x - win_min_x;
extent_y = win_max_y - win_min_y;

if length(binsize) > 1
    binsize = binsize(1);
end

%% We need to convert from cm to pixels for this function
win_max_x = win_max_x*pixels_metre / 100;
win_min_x = win_min_x*pixels_metre / 100;
win_max_y = win_max_y*pixels_metre / 100;
win_min_y = win_min_y*pixels_metre / 100;

posdata = posdata*pixels_metre / 100;



switch var2binby
    case 'position'
        [MX MY] = meshgrid((win_min_x:binsize:win_max_x) + (binsize/2), (win_min_y:binsize:win_max_y) + (binsize/2));        
        if isempty(find(MX>extent_x)) == 0 | isempty(find(MY >= extent_y)) == 0
            MX(:,size(MX,2)) = [];
            MX(size(MX,1),:) = [];
            MY(size(MY,1),:) = [];
            MY(:,size(MY,2)) = [];
        end        
        grid_values = {MX MY};
        
        if isempty(cluster_xy_index) == 1
            binned_array = zeros(size(MX));
        else
            binned_array = hist_nd(posdata(cluster_xy_index,:),grid_values);
        end
end

