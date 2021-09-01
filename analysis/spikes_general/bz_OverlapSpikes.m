function [overlap] = bz_OverlapSpikes(spikes,varargin)

% USAGE 
% [overlap] = bz_OverlapSpikes(spikes,varargin)
% Calculates the overlapping in Gaussian distribution of the different
% spikes

% INPUTS
%
%   spikes      - buzcode format .cellinfo. strct with the following fields
%                   .times
%
% <options>     - optinal list of property-value pairs
% 
% 
% OUTPUT
% 
% overlap - cellinfo 
% 
%
% Pablo Abad Pérez, 07/2021

p = inputParser;
addParameter(p,'show',true,@islogical);

parse(p,varargin{:});
show = p.Results.show;


numCells = spikes.numcells;

% overlapFP = zeros(numCells,numCells);
% overlapFN = zeros(numCells,numCells);

addpath('C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\OpenGSpike\QualityClustering\quality_measures')



sessionInfo = bz_getSessionInfo();

for i=1:length(sessionInfo.AnatGrps)
    
    spikesChan = find(spikes.shankID == i);
    disp(['Computing overlaptool for shank ', (num2str(i))])
    disp(['Total Neurons in Shank ', num2str(i), ': ' , num2str(length(spikesChan))])

    for ii=1:length(spikesChan)-1
        
        figure,
        
        for jj=1:length(spikesChan)-1

            spikes.UID(spikesChan(ii))
            spikes.UID(spikesChan(jj+1))

            if spikes.total(spikesChan(ii)) > 150 && spikes.total(spikesChan(jj+1)) > 150

                [confusion,x1,x2,w] = gaussian_overlap(spikes.waveforms.filt{spikesChan(ii)},spikes.waveforms.filt{spikesChan(jj+1)});
%                 [x1,x2,w] = testPlotClusterOverlapV3(spikes.waveforms.filt{spikesChan(ii)}, spikes.waveforms.filt{spikesChan(jj+1)},1);
            end
            
            ComFP = 1-((1-confusion(1,1))*(1-confusion(2,1)));  
            overlapFP{i}(ii,jj) = ComFP;
            ComFN = 1-((1-confusion(1,2))*(1-confusion(2,2)));
            overlapFN{i}(ii,jj) = ComFN;
            
            maIn = 1:1:(length(spikesChan)-1)*(length(spikesChan)-1);
            B = reshape(maIn,length(spikesChan)-1,(length(spikesChan)-1))';
            
            subplot((length(spikesChan)-1), length(spikesChan)-1, B(ii,(jj)));
            OverlapPLot(x1,x2,w,spikes.waveforms.filt{spikesChan(ii)},spikes.waveforms.filt{spikesChan(jj+1)});
        end



    end

end


end

