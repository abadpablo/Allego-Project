function [periodic] = bz_PeriodicPower(map,shuffling,varargin)
% Power spectrum in the 2D firing map 
% calculates the power spectrum of a firing matrix or in general 2D data matrix
% (c)2016 Jorge R Brotons-Mas and ValElsa Corps
% Adapted by Pablo Abad 27/07 to resemble buzCode 
% INPUT: 
%   map          map obtained using <a href="matlab:help Map">Map</a>
%
%   options      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'    By default pwd    
%    =========================================================================
%
%
%  OUTPUT:
%       
%
%
%
%
%
%
%
%
%
% Pablo Abad 27/07

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'nBins',20,@isnumeric);
addParameter(p,'showFig',false,@islogical);
addParameter(p,'saveFig',false,@islogical);
addParameter(p,'random',true,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
nBins = p.Results.nBins;
showFig = p.Results.showFig;
saveFig = p.Results.saveFig;
random = p.Results.random;

mapa = map.zUnSmooth;
mapa (isnan(mapa)) = 0;
mapa_sm = mapa-mean(mapa(:));
B = padarray(mapa_sm,[round((256-nBins)/2) round((256-nBins)/2)]); 
TF = fftshift(fft2(B));  
TFP=abs(TF.^2);
%%%%%%%%% Frequency %% NEED TO CHECK
max_num = max(TFP(:)) ;
[YY XX]=find(TFP == max_num); 

% Finds the max value in FFT and then calculate the distance to this point. The frequency is then 
% converted in pixeles and cm, so expressed in cycles by cm. 

evalua12 = [(256/2)-1 (256/2)-1;XX(1) YY(1)];%
dd12 = pdist(evalua12);
frec = (dd12/256)/0.025   % conversion to cycles per cm

thre= prctile(TFP(:),99) ; %thresholding for power 
thre50= prctile(TFP(:),50);

% Rotation of the original map to optain the power polar plot
gradStep=1;

TF2= TFP.*(TFP>thre);%>thre) ; CHECK this 121016
figure, imagesc(TF2);  
cleanmat=1;

while cleanmat==0
    disp('Eliminate local peaks ')
    BW = roipoly(TF2) ;
    e=double(BW);
    e=abs(e-1);
    TF2=e.*TF2;
    figure, imagesc(TF2);  figure, imagesc(TF2);  
    cleanmat=input('Clean matrix ok?');
end
 clear cleanmat

iii=1;
for grad=0:gradStep:360-1

    B = imrotate(abs(TF2) ,grad); %rotate
    rotate=B;
    [l1 l2]=size(rotate);
    BC(iii)=sum(rotate(round(l1/2),round(l2/2):end )/l2);
    iii=iii+1 ;
end

%  plot(BC)
%%
theta=pi*(0:gradStep:360-1)/180;
 
 %%%% Plotting and saving the figures por the FFT analysis
if showFig

    figure,imagesc(mapa);title(['Raw Firing Map Cell '])% original firing map
    figure, imagesc(abs(TF2));title(['Raw Power FFT^2 Cell '])%colormap gray
    figure;imagesc(abs(TF2).*(TFP>thre));title(['Thres 99% Power FFT^2 Cell '])
    figure;polar(theta,smooth(BC)');  title(['Polar FFT thrs 99% Cell '])
end
 
[maxPolar posPolar] = max(BC);
Orient = rad2deg(theta(posPolar));


%% Randomization

ratemapR2 = cell2mat(shuffling.maps);

if random
    [hh,jj] = size(ratemapR2(1).z);
    kk = size(ratemapR2,2);
else
    kk = 10;
end

BCR_out = zeros(360,kk);

for hj=1:kk
    
    mapR = ratemapR2(hj).zUnSmooth;
    mapa_smR = mapR- mean( mapR(:));   
    BR = padarray( mapa_smR  ,[round((256-nBins)/2)  round((256-nBins)/2) ] ); % padding to increase resolution
    TFR = fftshift(fft2(BR));% Power spectrum 
    TFRP= abs(TFR.^2);
    threE = prctile(TFRP(:),99) ; %thresholding for power 
    TFR2= TFRP.*(TFRP>threE);


%     figure, subplot(3,1,1), imagesc(mapa), title('Original Map');
%     subplot(3,1,2), imagesc (mapa_smR), title('Random Map');
%     subplot(3,1,3),imagesc(TFR2),title('Random Map FFT2');


    TFR_out(:,:,hj) = TFR;
    
    meanTFR_out=mean(TFR_out,3);
    meanTFR_out2=mean(meanTFR_out(:));
    threR_out(hj,1) =prctile(meanTFR_out2(:),95) ;
    iii=1;
    for grad=0:gradStep:360-1
        %Matrix B represents the rotated Matrix
        BR = imrotate(TFR2,grad); %rotate
        rotateR=BR;
        [l1 l2]=size(rotateR);
        BCR(iii)=sum(rotateR(round(l1/2),round(l2/2):end )/l2);
         iii=iii+1 ;
    end
    
    BCR_out(:,hj)= BCR';    
end


figure;imagesc(abs(TF2).*(TFP>thre));title(['Thres 99% Power FFT^2'])

% Polar plots 
threR992= prctile(BCR ,99); % value for interval. 
threR952= prctile(BCR ,95); % value for interval. 
threR992_out= prctile(BCR_out,99,2)

 
[maxPolar posPolar]=max(BC);
Orient=rad2deg(theta(posPolar));
ind=find(BC>threR992_out')
if isempty(ind)==0 
    Periodic=1
else
    Periodic=1
end


%% Reestructure into cell info data type

periodic.Periodic = Periodic;
periodic.TFP = TFP;
periodic.BC = BC;
periodic.meanTFR = meanTFR_out;
periodic.threR992 = threR992_out;
periodic.frec = frec;
periodic.maxPolar = maxPolar;
periodic.Orient = Orient;


end


