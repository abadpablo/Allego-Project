function [r] = bz_autoCorrPearson(firingMaps,id,varargin)

% Autocorrelation based on Pearson's product moment correlation coefficient
% with corrections for edge effects and unvisited locations.
%   z(x,y) is denoting the average fire rate of a cell at location(x,y). 
%   n is the summation of all pixels
%   tauXX and tauYY are the spacial lags
%
%
%
%
%
%
%
%
%
% Pablo Abad 2021, based on autoCorrPearson from Stefan Schaffelhofer 

% Default and Params

p = inputParser;

addParameter(p,'basepath',pwd,@isdir)

parse(p,varargin{:})

basepath = p.Results.basepath;

Matrix = firingMaps.rateMaps{id}{1};
Matrix(isnan(Matrix)) = 0;

%% Prepare empty Matrices
% amount of Rows of the Spike Density Map:
amountRowOrig = length(Matrix(:,1));
% amount of Columns of the choosen Spike Density Map:
amountColOrig = length(Matrix(1,:));

amountColCorr = amountRowOrig;
amountRowCorr = amountColOrig;

%r is the Matrix including the autocorrelation
r=zeros(amountColCorr*2,amountRowCorr*2);

%% Autocorrelation
for xx=1:1:amountRowCorr
%     waitbar(xx/amountRowCorr,wh);
    for yy=1:1:amountColCorr
    
    AA = Matrix(1:xx,1:yy);
    BB = Matrix(amountRowOrig-xx+1:amountRowOrig,amountColOrig-yy+1:amountColOrig);
    AAx= reshape(AA,1,numel(AA));
    BBy= reshape(BB,1,numel(BB));
    temp = corrcoef(AAx,BBy);
    if prod(size(temp))>=4 
        r(xx,yy) = temp(2,1);
    else
        r(xx,yy) = 0;
    end
    
    AA = Matrix(xx:amountRowOrig,yy:amountColOrig);
    BB = Matrix(1:amountRowOrig-xx+1,1:amountColOrig-yy+1);
    AAx= reshape(AA,1,numel(AA));
    BBy= reshape(BB,1,numel(BB));
    temp = corrcoef(AAx,BBy);
    
    if prod(size(temp))>=4 
        r(amountRowOrig+xx,amountColOrig+yy) = temp(2,1);
    else
        r(amountRowOrig+xx,amountColOrig+yy) = 0;
    end
    
    AA = Matrix(1:xx,amountColOrig-yy+1:amountColOrig);
    BB = Matrix(amountRowOrig-xx+1:amountRowOrig,1:yy);
    AAx= reshape(AA,1,numel(AA));
    BBy= reshape(BB,1,numel(BB));
    temp = corrcoef(AAx,BBy);
    
    if prod(size(temp))>=4 
        r(xx,2*amountColOrig+1-yy) = temp(2,1);
    else
        r(xx,2*amountColOrig+1-yy) = 0;
    end
      
    
    AA = Matrix(xx:amountRowOrig,1:amountRowOrig-yy+1);
    BB = Matrix(1:amountRowOrig-xx+1,yy:amountColOrig);
    AAx= reshape(AA,1,numel(AA));
    BBy= reshape(BB,1,numel(BB));
    temp = corrcoef(AAx,BBy);
    
    if prod(size(temp))>=4 
        r(xx+amountRowOrig,amountColOrig+1-yy) = temp(2,1);
    else
        r(xx+amountRowOrig,amountColOrig+1-yy) = 0;
    end

r(find(isnan(r)))=0; %eliminate NaN
    
    end
end




end

