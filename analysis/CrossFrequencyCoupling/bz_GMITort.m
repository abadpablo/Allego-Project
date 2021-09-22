function [MI] = bz_GMITort(lfp,phaserange,amprange,varargin)
% [GMI_Tort] = bz_GMITort(lfp,varargin)
%
%
% This function computes the Gamma Modulation Index as in Tort¿?* between a phase
% filtered signal and an amplitude filtered signal
%
%   INPUT
%
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                                                   lfp.channels
%
%   <options>       optional list of property-value pairs (see table below)
%
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     phaseCh      channel to compute phase. If empty take first channel
%     ampChans     channels to compute amplitude. If empty take first channel
%     method       ['wavelet'(default)|'hilbert']. Method to extract power
%        of specific bands, using wavelet (default) or hilbert
%
%     makePlot      default true
%     filterType    default 'fir1'. Method of filtering for the phase
%        bands, in case of method = 'hilbert' it also defines the filter of
%        bands in the amplitude range.
%     filterOrder   default 4. Order of the filter used for the phase
%        bands, in case of method = 'hilbert' it also defines the filter order of
%        bands in the amplitude range.
%     numBins       default 50. Number of bins to calculate the
%        Amplitude/Phase distribution.
%     perm_test     default false. To whether calculate or not a surrogate
%     test for CFC values.
%     alpha         default 0.05. Alpha for the surrogate test
%     num_inter     default 200.  Number of permutations
%     units         default MI. alternative: zscore. Whether the units are 
%     in MI or zscore after the permutation test
%    =========================================================================
%
%
%   OUTPUT
%
%
%
%
% Dependencies
%   bz_Filter
%   bz_WaveSpec
%
%   Pablo Abad Pérez 2021

%% Default params
p = inputParser();
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'sr',1250,@isnumeric);
addParameter(p,'phaseCh',lfp.channels(1),@isnumeric);
addParameter(p,'ampCh',lfp.channels(1),@isnumeric);
addParameter(p,'filterType','fir1',@isstring);
addParameter(p,'filterOrder',4,@isnumeric);
addParameter(p,'numBins',50,@isnumeric);
addParameter(p,'method','wavelet',@ischar);

addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'PhaseFreqVector',[1:2:20],@isnumeric);
addParameter(p,'AmpFreqVector',[10:5:200],@isnumeric);
addParameter(p,'PhaseFreq_BandWidth',4,@isnmuerical);
addParameter(p,'AmpFreq_BandWidth',20,@isnmuerical);
addParameter(p,'foldername',[],@isstr);
addParameter(p,'showFig',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'exist_file',false,@islogical);
addParameter(p,'speed_filter',false,@islogical); % to filter the signal for speed (HomeCage recordings do not have tracking information)¿?¿?¿?¿
addParameter(p,'noise_filter',false,@islogical);
addParameter(p,'epoch_analysis',false,@islogical); % Equal to previous parameter. If including in the analysis the filtering of the signals to remove TD ratio noise
% params for Chronux functions

parse(p,varargin{:})
basepath = p.Results.basepath;
sr = p.Results.sr;
phaseCh = p.Results.phaseCh;
ampCh = p.Results.ampCh;
filterType = p.Results.filterType;
filterOrder = p.Results.filterOrder;
numBins = p.Results.numBins;
method = p.Results.method;

nBins = p.Results.nBins;
PhaseFreqVector = p.Results.PhaseFreqVector;
AmpFreqVector = p.Results.AmpFreqVector;
PhaseFreq_BandWidth = p.Results.PhaseFreq_BandWidth;
AmpFreq_BandWidth = p.Results.AmpFreq_BandWidth;
foldername = p.Results.foldername;
showFig = p.Results.showFig;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;
exist_file = p.Results.exist_file;
speed_filter = p.Results.speed_filter;
noise_filter = p.Results.noise_filter;
epoch_analysis = p.Results.epoch_analysis;

%% Load sessionInfo
sessionInfo = bz_getSessionInfo(basepath);
shanks = sessionInfo.AnatGrps;
%% Input Channels
if isequal(length(phaseCh),length(ampCh)) && length(phaseCh) ~= 1
    disp('Same number of phaseChannels and amplitudeChannels. More than 1 channel')
elseif length(phaseCh) == 1 && ~isequal(length(phaseCh),length(ampCh))
    disp('One phaseChannel (ref) and more than 1 channel for amplitude')
elseif length(phaseCh) == 1 && length(ampCh) == 1
    disp('Computing GMI for only 1 channel')
end


%% Filtering the signals
filtered_phase = bz_Filter(lfp,'passband',phaserange,'filter',filterType,'order',filterOrder,'channels',phaseCh);
phase = filtered_phase.phase;

filtered_amp = bz_Filter(lfp,'passband',amprange,'filter',filterType,'order',filterOrder,'channels',ampCh);
amp = filtered_amp.amp;

%% ModIndex_V1 Tort
% Define phase bins
position = zeros(1,nBins); % This variable will get the beginning (not center) of each phase bin ( in rads)
winsize = 2*pi/nBins;

for jj=1:nBins
    position(jj) = -pi+(jj-1)*winsize;
end

% Now we search for a Phase-Amp relation between these frequencies by
% calculating the mean amplitude of the AmpFreq in each phase bin of the
% PhaseFreq

% Computing the mean amplitude in each phase

for ii=1:length(filtered_amp.channels)
    MeanAmp{ii} = zeros(1,nBins);
    for jj=1:nBins
        I = find(phase(:,ii) < position(jj)+winsize & phase(:,ii) >= position(jj));
        MeanAmp{ii}(jj) = mean(amp(I,ii));
    end
    % the center of each bin ( for plotting purposes) is position+winsize/2

    % quantifying the amount of amp modulation by means of a normalized entropy
    % index (Tort et al PNAS 2008):

    MI{ii} = (log(nBins)-(-sum((MeanAmp{ii}/sum(MeanAmp{ii})).*log((MeanAmp{ii}/sum(MeanAmp{ii}))))))/log(nBins);
end

if showFig
    if length(phaseCh) == 1 && length(ampCh) == 1
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        bar(linspace(1,720,nBins*2),[MeanAmp{1},MeanAmp{1}]/sum(MeanAmp{1}),'k')
        xlim([0 720])
        set(gca,'xtick',0:360:720)
        xlabel('Phase (Deg)')
        ylabel('Amplitude')
        title(['MI = ' num2str(MI{1}), ' Ch:', num2str(phaseCh)])

        if saveFig
            if ~isempty(foldername)
                saveas(gcf,['lfpAnalysisFigures\GMI_Tort_v1.',foldername,'_',num2str(amprange(1)), num2str(amprange(end)),'.png'])
            else
                saveas(gcf,['lfpAnalysisFigures\GMI_Tort_v1','_',num2str(amprange(1)), num2str(amprange(end)),'.png'])
            end
        end
    else
        meanamp = cell2mat(MeanAmp);
        minamp = min(meanamp/sum(meanamp));
        maxpamp = max(meanamp/sum(meanamp));
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj=1:size(shanks,2)
            count = jj;
            for ii=1:length(shanks(jj).Channels)
                subplot(sessionInfo.nChannels/sessionInfo.nElecGps,sessionInfo.nChannels/sessionInfo.nElecGps,count);
                bar(linspace(1,720,nBins*2),[MeanAmp{shanks(jj).Channels(ii)+1} , MeanAmp{shanks(jj).Channels(ii)+1}] / sum(MeanAmp{shanks(jj).Channels(ii)+1}),'k')
                xlim([0 720])
                ylim([0 0.025])
                if jj == 1 && ii ~= length(shanks(jj).Channels)
                    ylabel('Amplitude');
                    set(gca,'XTIck',[]);
                elseif jj == 1 && ii == length(shanks(jj).Channels)
                    xlabel('Phase (Deg)');
                    set(gca,'xtick',0:360:720)
                    ylabel('Amplitude');
                elseif jj~= 1 && ii == length(shanks(jj).Channels)
                    xlabel('Phase (Deg)')
                    set(gca,'xtick',0:360:720)
                    set(gca,'YTick',[])
                else
                    set(gca,'XTick',[])
                    set(gca,'YTIck',[])
                end
                title(['MI = ' num2str(MI{shanks(jj).Channels(ii)+1}), ' Ch:', num2str(shanks(jj).Channels(ii)+1)])
                count = count + size(shanks,2);
            end
        end
        
        if saveFig
            if ~isempty(foldername)
                saveas(gcf,['lfpAnalysisFigures\GMI_Tort_v1_allCh.',foldername,'_',num2str(amprange(1)),num2str(amprange(end)),'.png'])
            else
                saveas(gcf,['lfpAnalysisFigures\GMI_Tort_v1_allCh','_',num2str(amprange(1)),num2str(amprange(end)),'.png'])
            end
        end
    end
end


mi = MI;
GMI = [];
GMI.MI = mi;
GMI.MeanAmp = MeanAmp;

clear MI
clear MeanAmp

%% Compute MI and Comodulogram ModIndex_v2

PhaseFreqVector = phaserange(1):2:phaserange(2);
AmpFreqVector = amprange(1):5:amprange(2);
PhaseFreq_BandWidth = 2;
AmpFreq_BandWidth = 5;


% Fitering and Hilbert Transform

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj);
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq = bz_Filter(lfp,'passband',[Pf1 Pf2],'filter',filterType,'order',filterOrder,'channels',phaseCh);
    PhaseFreqTransformed{jj} = PhaseFreq.phase; % getting the phase time series
end

for ii=1:length(AmpFreqVector)
   Af1 = AmpFreqVector(ii); 
   Af2 = Af1+AmpFreq_BandWidth; 
   AmpFreq = bz_Filter(lfp,'passband',[Af1 Af2],'filter',filterType,'order',filterOrder,'channels',phaseCh);
   AmpFreqTransformed{ii} = AmpFreq.amp;
end

%% Compute MI and Comodulogram


for ch = 1:length(PhaseFreq.channels)
    counter1=0;
    for ii=1:length(PhaseFreqVector)        
    counter1=counter1+1;

    Pf1 = PhaseFreqVector(ii);
    Pf2 = Pf1+PhaseFreq_BandWidth;

    counter2=0;
        for jj=1:length(AmpFreqVector)
        counter2=counter2+1;
            Af1 = AmpFreqVector(jj);
            Af2 = Af1+AmpFreq_BandWidth;

            [MI{ch}{jj},MeanAmp{ch}{jj}]=ModIndex_v2(PhaseFreqTransformed{ii}(:,ch), AmpFreqTransformed{jj}(:,ch), position);
            Comodulogram{ch}(counter1,counter2)=MI{ch}{jj};
        end
    end
end
    
%% Plotting Comodulogram

if showFig
    if length(phaseCh) == 1 && length(ampCh) == 1
        
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        clf
        contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram{1}',30,'lines','none')
        set(gca,'fontsize',14)
        ylabel('Amplitude Frequency (Hz)')
        xlabel('Phase Frequency (Hz)')
        colormap(jet)
        colorbar
        if saveFig
            if ~isempty(foldername)
                saveas(gcf,['lfpAnalysisFigures\Comodulogram_Tort.',foldername,'_',num2str(amprange(1)), num2str(amprange(end)),'.png'])
            else
                saveas(gcf,['lfpAnalysisFigures\Comodulogram_Tort','_',num2str(amprange(1)), num2str(amprange(end)),'.png'])
            end
        end
    else
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj=1:size(shanks,2)
            count = jj;
            for ii=1:length(shanks(jj).Channels)
                subplot(sessionInfo.nChannels/sessionInfo.nElecGps,sessionInfo.nChannels/sessionInfo.nElecGps,count)
                contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram{shanks(jj).Channels(ii)+1}',30,'lines','none')
%                 set(gca,'fontsize',14)
                colormap(jet)
                colorbar
                if jj==1 && ii ~= length(shanks(jj).Channels)
                    ylabel('Amplitude Frequency (Hz)')
                    set(gca,'XTick',[])
                elseif jj == 1 && ii == length(shanks(jj).Channels)
                    xlabel('Phase Frequency (Hz)')
                    ylabel('Amplitude Frequency (Hz)')
                elseif jj~=1 && ii == length(shanks(jj).Channels)
                    xlabel('Phase Frequency (Hz)')
                    set(gca,'YTick',[])
                else
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                end
                title(['Comodulogram Ch:', num2str(shanks(jj).Channels(ii)+1)])
                count = count + size(shanks,2);
            end
        end             
    end
end

%% Save matlab structure

if ~isempty(foldername)
    GMI.foldername = foldername;
end
GMI.Channels.phaseCh = phaseCh;
GMI.Channels.ampCh = ampCh;

GMI.Comodulogram = Comodulogram;
GMI.PhaseFreqVector = PhaseFreqVector;
GMI.PhaseFreq_BandWidth = PhaseFreq_BandWidth;
GMI.AmpFreqVector = AmpFreqVector;
GMI.AmpFreq_BandWidth = AmpFreq_BandWidth;


if saveMat
    if ~isempty(folder)
        try
            save([basepath filesep sessionInfo.FileName, '.', foldername, '.GMI_Tort_', num2str(amprange(1)), num2str(amprange(end)),'.lfp.mat'],'GMI')
        catch
            save([basepath filesep sessionInfo.FileName, '.', foldername, '.GMI_Tort_', num2str(amprange(1)), num2str(amprange(end)),'.lfp.mat'],'GMI','-v7.3')
        end
    else
        try
           save([basepath filesep sessionInfo.FileName,'.GMI_Tort_', num2str(amprange(1)), num2str(amprange(end)),'.lfp.mat'],'GMI') 
        catch
            save([basepath filesep sessionInfo.FileName,'.GMI_Tort_', num2str(amprange(1)), num2str(amprange(end)),'.lfp.mat'],'GMI','-v7.3')
        end
    end
end

end

