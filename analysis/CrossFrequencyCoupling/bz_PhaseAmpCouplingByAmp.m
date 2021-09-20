function [comod ] = bz_PhaseAmpCouplingByAmp(lfp,sig1band,sig2band,varargin)
% [comod] = bz_PhaseAmpCouplingByAmp(lfp,sig1band,sig2band,varargin)
%calculates phase amplitude coupling between the phase of signal 1 and the
%amplitude of signal 2, with respect to the amplitude of signal 1.
%
%
%%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                                                   lfp.channels
%                   
%    sig1band [min max] - band of the reference signal (Hz)
%    sig2band [min max] - band of the tested signal (Hz)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     phaseCh      channel to compute phase. If empty takes first channel
%     ampCh     channels to compute amplitude. If empty takes first channel
%     method       ['hilbert'(default)|'wavelet']. Method to extract power
%        of both signals, using wavelet or hilbert (default)
%     
%     makePlot      default true
%     filterType    default 'fir1'. Method of filtering for the phase
%        bands, in case of method = 'hilbert' it also defines the filter of
%        bands in the amplitude range.
%     filterOrder   default 4. Order of the filter used for the phase
%        bands, in case of method = 'hilbert' it also defines the filter order of
%        bands in the amplitude range.
%     ampNumBins    default 50. Number of bins to calculate the
%        amplitude distribution.
%     phaseNumBins  default 50. Number of bins to calculate the phase
%     distribution.
%    =========================================================================
%
% OUTPUT
% comod.ampbins = ampbins;
% comod.phasebins = phasebins;
% comod.sig2powerskew = sig2powerskew;
% comod.sig2prefangle = sig2prefangle;
% comod.phaseamphist = phaseamphist;
% comod.params        Collection of parameters used to calculate this comod
%
%Dependencies
%   bz_Filter
%   bz_WaveSpec
%
%
%DLevenstein Fall 2016
%Improved documentation and I/O. Added extraction of signal features - Eliezyer de Oliveira 2019
%
%%


if ~bz_isLFP(lfp)
    error('Not following the lfp structure required, see documentation')
end

p = inputParser;
addParameter(p,'ampCh',lfp.channels(1),@isnumeric);
addParameter(p,'phaseCh',lfp.channels(1),@isnumeric);
addParameter(p,'filterType','fir1',@ischar);
addParameter(p,'filterOrder',4,@isnumeric);
addParameter(p,'ampNumBins',50,@isnumeric);
addParameter(p,'phaseNumBins',50,@isnumeric);
addParameter(p,'intervals',[0 inf],@isnumeric);
addParameter(p,'makePlot',true,@islogical);
addParameter(p,'method','hilbert',@ischar);
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'foldername',[],@isstr);
addParameter(p,'exist_file',false,@islogical);


parse(p,varargin{:});
ampCh = p.Results.ampCh;
phaseCh = p.Results.phaseCh;
makePlot = p.Results.makePlot;
filterType = p.Results.filterType;
filterOrder = p.Results.filterOrder;
ampNumBins = p.Results.ampNumBins;
phaseNumBins = p.Results.phaseNumBins;
intervals = p.Results.intervals;
method = p.Results.method;
basepath = p.Results.basepath;
saveFig = p.Results.saveFig;
saveMat = p.Results.saveMat;
foldername = p.Results.foldername;
exist_file = p.Results.exist_file;

% Load sessionInfo
sessionInfo = bz_getSessionInfo(basepath);

% Loading .mat file if exists
if ~isempty(foldername)
    if ~isempty(dir([basepath filesep sessionInfo.FileName, '.',foldername, '.*PhaseAmpCouplingByAmp_',num2str(sig2band(1)),num2str(sig2band(end)),'.lfp.mat']))
        disp(['PhaseAmpCouplingByAmp_', num2str(sig2band(1)),num2str(sig2band(end)),' ' , foldername, ' already detected ! Loading file.'])
        file = dir([basepath filesep sessionInfo.FileName ,'.',foldername, '.*PhaseAmpCouplingByAmp_',num2str(sig2band(1)),num2str(sig2band(end)),'.lfp.mat'])
        load(file.name);
        exist_file = true;
        %return
    end
else
    if ~isempty(dir([basepath filesep 'PhaseAmpCouplingByAmp_',num2str(sig2band(1)),num2str(sig2band(end)),'.lfp.mat']))
        disp('Coherogram already detected! Loading file.')
        file = dir([basepath filesep 'PhaseAmpCouplingByAmp_',num2str(sig2band(1)),num2str(sig2band(end)),'.lfp.mat']);
        load(file.name);
        exist_file = true;
        %return
    end
end
    


%%
if ~exist_file
    phasebins = linspace(-pi,pi,phaseNumBins+1);
    phasebins=phasebins(1:end-1)+diff(phasebins(1:2));
    ampbins = linspace(-2.5,2.5,ampNumBins);
    % ampNumBins = length(ampbins);

    %%
    sig1 = bz_Filter(lfp,'passband',sig1band,'filter',filterType,'order',filterOrder,'channels',phaseCh);
    sig1phase = sig1.phase;
    switch method
        case 'wavelet'
            sig1 = bz_Wavelet(lfp,'frange',sig1band,'nfreqs',1,'chanID',phaseCh);
            sig1amp = abs(sig1.data);

            sig2 = bz_Wavelet(lfp,'frange',sig1band,'nfreqs',1,'chanID',ampCh);
            sig2amp = abs(sig2.data);
        case 'hilbert'
            sig1 = bz_Filter(lfp,'passband',sig1band,'filter',filterType,'order',filterOrder,'channels',phaseCh);
            sig1amp = sig1.amp;

            sig2 = bz_Filter(lfp,'passband',sig2band,'filter',filterType,'order',filterOrder,'channels',ampCh);
            sig2amp = sig2.amp;

    end

    sig1amp = zscore(sig1amp);
    sig2amp = zscore(sig2amp);

    %%
    sig1binpower = interp1(ampbins,ampbins,sig1amp,'nearest');
    sig1binphase = interp1(phasebins,phasebins,sig1phase,'nearest');


    powerhist = zeros(ampNumBins,ampNumBins);
    sig2prefangle = zeros(ampNumBins,1);
    sig2powerskew = zeros(ampNumBins,1);
    for bb = 1:ampNumBins
        for bbb = 1:phaseNumBins
            ampwintimes = sig2amp(sig1binpower==ampbins(bb) & sig1binphase==phasebins(bbb));
            phaseamphist(bb,bbb) = mean(ampwintimes);
        end
        sig2powerskew(bb) = mean(sig2amp(sig1binpower==ampbins(bb)).*exp(1i.*sig1phase(sig1binpower==ampbins(bb))));
        sig2prefangle(bb) = angle(sig2powerskew(bb));
        sig2powerskew(bb) = abs(sig2powerskew(bb));
    end


    comod.ampbins = ampbins;
    comod.phasebins = phasebins;
    comod.sig2powerskew = sig2powerskew;
    comod.sig2prefangle = sig2prefangle;
    comod.phaseamphist = phaseamphist;
    comod.params.sig1band = sig1band;
    comod.params.sig2band = sig2band;
    comod.params.method = method;
    comod.params.filterType = filterType;
    comod.params.filterOrder = filterOrder;
    comod.sig1amp = sig1amp;
    comod.sig2amp = sig2amp;
    comod.sig1phase = sig1phase;
    % Adding channels for analysis
    comod.phaseCh = phaseCh;
    comod.ampCh = ampCh;
    if ~isempty(foldername)
        comod.foldername = foldername;
    end
    %% Figure


    if makePlot
        rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
        plotx = linspace(-pi,3*pi,100);
        if ~isempty(foldername)
            figure('Name',foldername)
        else
            figure('Name',basepath)
        end
        subplot(2,2,1)
        hold on
        imagesc(phasebins,ampbins,phaseamphist)
        imagesc(phasebins+2*pi,ampbins,phaseamphist)
        plot(sig2prefangle,ampbins,'.k')
        plot(sig2prefangle+2*pi,ampbins,'.k')
        plot(plotx,cos(plotx),'k')
        colormap(gca,rwbcolormap)
        axis xy
        axis tight
        ColorbarWithAxis([-0.5 0.5],['Mean Amp.'])
        caxis([-0.5 0.5])
        %  xlim([-pi 3*pi]);ylim(ampbins([1 end]))
        xlabel('Signal 1 Phase');ylabel('Signal 1 Amp (Z)')
        subplot(4,2,2)
        plot(ampbins,sig2powerskew,'k','LineWidth',1)
        xlabel('Signal 1 Amp. (Z)');
        ylabel('Phase-Amp. Modulation (mrl)')
        axis tight
        subplot(4,2,6)
        histogram(sig1amp,ampbins)
        xlabel('Signal 1 Amp. (Z)');
        ylabel('Occupancy')
        axis tight
        title('Signal1/2 Amp. Distributions')

        subplot(4,2,8)
        histogram(sig2amp,ampbins)
        xlabel('Signal 2 Amp. (Z)');
        ylabel('Occupancy')
        axis tight

        subplot(2,2,3)
        gasphist = hist3([sig1amp,sig2amp],{ampbins,ampbins});
        imagesc(ampbins,ampbins,gasphist)
        axis xy
        xlabel('Signal 1 Power');ylabel('Signal 2 Power')


        if saveFig
            if ~isempty(foldername)
               saveas(gcf,['lfpAnalysisFigures\PhaseAmpCouplingByAmp.',foldername,'_',num2str(sig2band(1)),num2str(sig2band(end)),'.png'])
            else
                saveas(gcf,['lfpAnalysisFigures\PhaseAmpCouplingByAmp_',num2str(sig2band(1)),num2str(sig2band(end)),'.png'])
            end
        end

    end

    if saveMat
        if ~isempty(foldername)
            try
                save([basepath filesep sessionInfo.FileName,'.',foldername,'.PhaseAmpCouplingByAmp','_',num2str(sig2band(1)),num2str(sig2band(end)),'.SubSession.lfp.mat'],'comod')
            catch
                save([basepath filesep sessionInfo.FileName,'.',foldername,'.PhaseAmpCouplingByAmp','_',num2str(sig2band(1)),num2str(sig2band(end)),'.SubSession.lfp.mat'],'comod','-v7.3')
            end
        else
            try
                save([basepath filesep sessionInfo.FileName,'.','PhaseAmpCouplingByAmp','_',num2str(sig2band(1)),num2str(sig2band(end)),'.lfp.mat'],'comod')
            catch
                save([basepath filesep sessionInfo.FileName,'.','PhaseAmpCouplingByAmp','_',num2str(sig2band(1)),num2str(sig2band(end)),'.lfp.mat'],'comod')
            end
        end


    end
else
    
    if makePlot
        rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
        plotx = linspace(-pi,3*pi,100);
        if ~isempty(foldername)
            figure('Name',foldername)
        else
            figure('Name',basepath)
        end
        subplot(2,2,1)
        hold on
        imagesc(comod.phasebins,comod.ampbins,comod.phaseamphist)
        imagesc(comod.phasebins+2*pi,comod.ampbins,comod.phaseamphist)
        plot(comod.sig2prefangle,comod.ampbins,'.k')
        plot(comod.sig2prefangle+2*pi,comod.ampbins,'.k')
        plot(plotx,cos(plotx),'k')
        colormap(gca,rwbcolormap)
        axis xy
        axis tight
        ColorbarWithAxis([-0.5 0.5],['Mean Amp.'])
        caxis([-0.5 0.5])
        %  xlim([-pi 3*pi]);ylim(ampbins([1 end]))
        xlabel('Signal 1 Phase');ylabel('Signal 1 Amp (Z)')
        subplot(4,2,2)
        plot(comod.ampbins,comod.sig2powerskew,'k','LineWidth',1)
        xlabel('Signal 1 Amp. (Z)');
        ylabel('Phase-Amp. Modulation (mrl)')
        axis tight
        subplot(4,2,6)
        histogram(comod.sig1amp,comod.ampbins)
        xlabel('Signal 1 Amp. (Z)');
        ylabel('Occupancy')
        axis tight
        title('Signal1/2 Amp. Distributions')

        subplot(4,2,8)
        histogram(comod.sig2amp,comod.ampbins)
        xlabel('Signal 2 Amp. (Z)');
        ylabel('Occupancy')
        axis tight

        subplot(2,2,3)
        gasphist = hist3([comod.sig1amp,comod.sig2amp],{comod.ampbins,comod.ampbins});
        imagesc(comod.ampbins,comod.ampbins,gasphist)
        axis xy
        xlabel('Signal 1 Power');ylabel('Signal 2 Power')


    end    
end

end

