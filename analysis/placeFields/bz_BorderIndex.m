function [borderIndex] = bz_BorderIndex(map,varargin)
% Computes Border Index
% Adapted from borderindex by Pablo Abad 23/07 to resemble buzCode
%
% USAGE
%   [borderIndex] = bz_BorderIndex(map,<options>);
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
%       borderIndex.west        west
%       borderIndex.east        east
%       borderIndex.south       south
%       borderIndex.north       north

% Pablo Abad 07/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'minFSize',12,@isnumeric);
addParameter(p,'smoothing',true,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
minFSize = p.Results.minFSize;
smoothing = p.Results.smoothing;

if smoothing
    Map = map.z;
else
    Map = map.zUnSmooth;
end

disp('Nan changed by 0')

indnan=(find(isnan(Map)));
Map(indnan)=0;

[bn jn] = size(Map); % to determine the number of bins to consider a border firing.


% MINPFSIZE = 12;

WALLDIST = round(jn/3);
BIASVALUE = 0;
MAPDATA = Map;

mapdata = MAPDATA((1+BIASVALUE):(end-BIASVALUE),(1+BIASVALUE):(end-BIASVALUE));
MaxF1=max(max(mapdata));
Mf1 = nanmean(nanmean(mapdata));
mapV = mapdata(:);
Sd1 = nanstd(mapV);
Se1 = Sd1/(sqrt(bn*jn));
Lim1 = Mf1+Se1;
thFF1 = (Lim1*100)/MaxF1;
tresFF = Lim1/MaxF1;
zzz = max(mapdata(:));


% if zzz>0 
if zzz>1 %max  frec treshold
    mapdata(mapdata<(tresFF*zzz)) = 0;
    rmap = mapdata;
    mapdata(mapdata>0) = 1;
end

% We need to convert Nan values to 0
[nan_values] = isnan(mapdata);
[a1,a2] = find(nan_values == 1);
mapdata(a1,a2) = 0;

rmap(a1,a2) = 0;


vvv = bwlabel(squeeze(mapdata));
PLACEFIELDS = [];
% figure;pcolor(vvv);shading flat
ii = 0;
for jj=1:max(vvv(:))
    ttt = sum(vvv(:)==jj);
    if(ttt>=minFSize)
        ii = ii + 1;
        PLACEFIELDS.size{ii} = ttt;
        PLACEFIELDS.data{ii} = mapdata;
        PLACEFIELDS.data{ii}(vvv~=jj) = 0;
    else
        mapdata(vvv==jj) = 0;
        rmap(vvv==jj) = 0;
    end
end

cM = 0;
dM = 0;

% if ii>0 
if ii>0 && zzz>1
    disty = WALLDIST;
    minsize = min(size(mapdata))/2;
    xsize = size(mapdata,1);
    ysize = size(mapdata,2);
    rmap = rmap/(sum(rmap(:)));
    
    cM1 = 0;
    cM2 = 0;
    cM3 = 0;
    cM4 = 0;
    
    for jj=1:ii
        aa = sum(PLACEFIELDS.data{jj}(1:end,1:disty),2); % west
        aa(aa>0) = 1;
        cM1 = max([cM1 ...
            sum(aa)/xsize]);
        aa = sum(PLACEFIELDS.data{jj}(1:end,(end-disty+1):end),2); % east
        aa(aa>0) = 1;
        cM2 = max([cM2 ...
            sum(aa)/xsize]);
        aa = sum(PLACEFIELDS.data{jj}(1:disty,1:end),1); % north
        aa(aa>0) = 1;
        cM3 = max([cM3 ...
            sum(aa)/xsize]);
        aa = sum(PLACEFIELDS.data{jj}((end-disty+1):end,1:end),1);% south
        aa(aa>0) = 1;
        cM4 = max([cM4 ...
            sum(aa)/xsize]);
        
        
%         
%         cM2 = max([cM3 ...
%             sum(sum(PLACEFIELDS.data{ii}(1:end,(end-disty+1):end)))/(xsize*disty)]);
%         cM3 = max([cM3 ...
%             sum(sum(PLACEFIELDS.data{ii}(1:disty,1:end)))/(ysize*disty)]);
%         cM4 = max([cM4 ...
%             sum(sum(PLACEFIELDS.data{ii}((end-disty+1):end,1:end)))/(ysize*disty)]);
        
    end
    
    cM = [cM1 cM2 cM3 cM4];
    
    distmap = rmap;
    
    for xx = 1:xsize
        for yy = 1:ysize
            distmap(xx,yy) = rmap(xx,yy) * min([xx-1 yy-1 xsize-xx+1 ysize-yy+1]);
        end
    end
    
    dM = sum(distmap(:))/minsize;
    
    
    BORDERSCORE = (cM-dM)./(cM+dM);
    
    
else
    BORDERSCORE = -1;
end

%% output

% CCC you will get the BORDERVALUE for the 4 walls ([south north  west east])
% % 
% % DDD are the CM values
% % 
% % and EEE is the DM value
ss=size(BORDERSCORE,2);
if ss>1 %w e n s..
west=BORDERSCORE(1);
east=BORDERSCORE(2);
north=BORDERSCORE(3);% change order 6/6/2018
south=BORDERSCORE(4);% change order 6/6/2018
elseif BORDERSCORE==-1
west=BORDERSCORE ;
east=BORDERSCORE;
south=BORDERSCORE;
north=BORDERSCORE;
end
maxBor1=max(BORDERSCORE);

borderIndex = [];
borderIndex.west = west;
borderIndex.east = east;
borderIndex.north = north;
borderIndex.south = south;
borderIndex.maxBorder = maxBor1;
