function [firingField]=bz_FiringFieldSize(map,varargin)

% Computes Firing Field Size
% Adapted from FiringFieldSize by Pablo Abad 23/07 to resemble buzCode
%
% USAGE
%   [firingField] = bz_FiringFieldSize(map,<options>);
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
%       firingField.numFF      number of firing fields
%       firingField.FFArea     area of the firing field
%       firingField.FFAreatot  % of the total area       
%       firingField.patchs
%       firingField.patchsArea
%       firingField.patchsAreatot
%       firingField.FFArevspatchAr
%       firingField.TotSizeRat
%       firingField.ponx
%       firingField.pony
%       firingField.MaxFfre
%       firingField.thFF1
%
% Pablo Abad 07/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'smoothing',true,@islogical);
addParameter(p,'minFFSize',8,@isnumeric); 
addParameter(p,'specificity',true,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
smoothing = p.Results.smoothing;
minFFSize = p.Results.minFFSize;
specificity = p.Results.specificity;

%%

if smoothing
    Mapf = map.z;
    Map = map.z;
else
    Mapf = map.zUnSmooth;
    Map = map.zUnSmooth;
end

indnan=(find(isnan(Mapf)));
Mapf(indnan)=0;
firmat=Mapf;

MaxF1=max(max(Mapf)); % max firing rate
% we need to take into account that some values are NaNs
Mf1 = nanmean(nanmean(Map));
MapfV = Mapf(:);
Sd1 = nanstd(MapfV); % standard deviation of the firing matrix %

[orr abb] = size(Mapf);

totalbins = abb*orr;
Se1 = Sd1;%/(sqrt(orr));
Lim1 = Mf1+Se1;%;Sdthe
thFF1 = (Lim1*100)/MaxF1;
tresFF = Lim1/MaxF1;
zzz = max(Mapf(:));

if zzz>1 %max  frec treshold
    Mapf(Mapf<(tresFF*zzz)) = 0;
    rMapf = Mapf;
    Mapf(Mapf>0) = 1;
end

[mm pp]=size(Mapf);

totBins = mm*pp;
minFFSize = 8; %(4*totBins)/100;
% disp('Warning % OF arena active as treshold')

% MINPFSIZE=15
vvv = bwlabel(squeeze(Mapf));
firingField = [];

% figure;pcolor(vvv);shading flat
%%
ii = 0;
for jj=1:max(vvv(:))
    ttt = sum(vvv(:)==jj);
    indV = find(vvv~=jj);
    firmod = firmat ;
    firmod(indV)=0;
    MaxF = max( max(firmod) );
    if(ttt>=minFFSize) & MaxF>1
        ii = ii + 1;
        firingField.size{ii} = ttt;
        firingField.sizeperc{ii} = (ttt*100)/totalbins;
        firingField.sizearea{ii} =  (6.25)*ttt;
        firingField.data{ii} = vvv;
        firingField.data{ii}(vvv~=jj) = 0;
        
        [posy posx ] =find(firmat==MaxF);
        firingField.positionx{ii}=posx(1); 
        firingField.positiony{ii}=posy(1);
        firingField.MaxF{ii}=MaxF;
        firingField.Nosize{ii}=0;
        firingField.Nosizeperc{ii} = 0;
        firingField.Nosizearea{ii} = 0;
%         figure; imagesc(firmod)
        
    else
         ii = ii + 1;
        firingField.Nosize{ii}=ttt;
        firingField.Nosizeperc{ii} = (ttt*100)/totalbins;
        firingField.Nosizearea{ii} =  (6.25)*ttt;
        % No size firing field parameters
        firingField.size{ii} = 0;
        firingField.sizeperc{ii} = 0;
        firingField.sizearea{ii} =  0;
        firingField.data{ii} = vvv;
        firingField.data{ii}(vvv~=jj) = 0;
        
        [posy posx]=find(firmat==0);
        firingField.positionx{ii}=0; 
        firingField.positiony{ii}=0;
        firingField.MaxF{ii}=0;
        firingField.Nosize{ii}=ttt;

        
        if specificity
            T = sum(sum(map.time));
            if T == 0
                stats.specificity = 0;
            else
                occupancy = map.time/(T+eps);
                m = sum(sum(map.count))/(sum(sum(map.time)+eps));
                if m == 0
                    stats.specificity = 0;
                else
                    logArg = map.count/m;
                    logArg(logArg <= 1) = 1;
                    stats.specificity = sum(sum(map.count.*log2(logArg).*occupancy))/m;
                    firingField.specificity{ii} = stats.specificity;
                end
            end
        end



    end
end


%%
if  ~isempty(firingField)
    iNumFF = find(cell2mat(firingField.size));
    NumFF = length(iNumFF);
    FFArea = mean(cell2mat(firingField.sizeperc(iNumFF)));
    FFAreatot = sum(cell2mat(firingField.sizeperc(iNumFF)));
    ipatchs = find((cell2mat(firingField.size))==0) ;
    patchs = length(ipatchs);
    patchsArea = mean(cell2mat(firingField.Nosizeperc(ipatchs)));
    if isnan(patchsArea)==1
        patchsArea=0;
    end
    patchsAreatot=sum(cell2mat(firingField.Nosizeperc(ipatchs)));
    [MaxFfre posmax]=max(cell2mat(firingField.MaxF));
    ponx=cell2mat(firingField.positionx(posmax));
    pony=cell2mat(firingField.positiony(posmax));
    Ffvspatch=NumFF/patchs;
    if patchsArea>0
        FFArevspatchAr=(FFArea/NumFF)/(patchsArea/patchs);
        TotSizeRat=(FFAreatot/patchsAreatot)/FFAreatot;
    else
        FFArevspatchAr=(FFArea/NumFF)/1;
        TotSizeRat=(FFAreatot/1)/FFAreatot;
    end
else
    
    iNumFF=0;
    NumFF=0;
    FFArea=0;
    FFAreatot=0;
    ipatchs=0 ;
    patchs=0;
    patchsArea=0;    
    patchsAreatot=0;
    FFArevspatchAr=0;
    TotSizeRat=0;
    ponx=0;
    pony=0;
    MaxFfre=0; 
end
%%
firingField.numFF = NumFF;
firingField.areaFF = FFArea;
firingField.areaTotFF = FFAreatot;
firingField.patchs = patchs;
firingField.patchsArea = patchsArea;
firingField.patchsAreaTot = patchsAreatot;
firingField.FFArevspatchAr = FFArevspatchAr;
firingField.TotSizeRat = TotSizeRat;
firingField.ponx = ponx;
firingField.pony = pony;
firingField.MaxFfre = MaxFfre;



