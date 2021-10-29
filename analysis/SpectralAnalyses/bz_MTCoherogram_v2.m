
function [coherogram] = bz_MTCoherogram_v2(lfp1,lfp2,varargin)

%MTCoherogram - Compute LFP coherogram by multi-taper estimation.
%
%  USAGE
%
%    [coherogram,phase,t,f] = MTCoherogram(lfp1,lfp2,<options>)
%
%    lfp1,lfp2      wide-band LFPs (one channel each).
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = from timestamps if
%                   available, otherwise 1250Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'window'      duration (in s) of the time window (default = 5)
%     'overlap'     overlap between successive windows (default = window/2)
%     'step'        step between successive windows (default = window/2)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help cohgramc">cohgramc</a>) (default = 0)
%     'show'        plot results (default = 'off')
%     'cutoffs'     cutoff values for color plot (default = [0 1])
%    =========================================================================
%
%  NOTES
%
%    The LFP can be provided either as a time stamped matrix (list of time-voltage
%    pairs), or as a voltage vector - in which case the frequency must be specified.
%
%    The time displacement between successive short time coherences can be supplied
%    either as a 'step' (explicit time difference) or as an 'overlap' (between
%    successive time windows).
%
%  OUTPUT
%
%    coherogram     coherogram magnitude
%    phase          coherogram phase
%    t              time bins
%    f              frequency bins
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.
%
%  SEE
%
%    See also MTCoherence, MTSpectrum, MTSpectrogram, PlotColorMap.

% Copyright (C) 2010-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% Function modified by Pablo Abad to incluse inputParser
% 
% Modified by Pablo Abad 2021 to include filtering by TD ratio and speed 
% (if tracking exists)

% Make sure chronux is installed and functional
CheckChronux('cohgramc');

% Default Params
p = inputParser();
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'f',1250,@isnumeric);
addParameter(p,'frequency',[],@isnumeric);
addParameter(p,'window',5,@isnumeric);
addParameter(p,'range',[],@isnumeric);
addParameter(p,'overlap',[],@isnumeric);
addParameter(p,'step',[],@isnumeric);
addParameter(p,'show','on',@isstr);
addParameter(p,'tapers',[3 5], @isnumeric);
addParameter(p,'pad',0,@isnumeric);
addParameter(p,'cutoffs',[0 1],@isnumeric);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'foldername',[],@isstr);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'exist_file',false,@islogical);
addParameter(p,'filtering',true,@islogical);


parse(p,varargin{:})

basepath = p.Results.basepath;
f = p.Results.f;
frequency = p.Results.frequency;
window = p.Results.window;
range = p.Results.range;
overlap = p.Results.overlap;
step = p.Results.step;
show = p.Results.show;
tapers = p.Results.tapers;
pad = p.Results.pad;
cutoffs = p.Results.cutoffs;
saveFig = p.Results.saveFig;
foldername = p.Results.foldername;
saveMat = p.Results.saveMat;
exist_file = p.Results.exist_file;
filtering = p.Results.filtering;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
end

% Check parameter sizes
if size(lfp1.data,2) ~= 1 && size(lfp1.data,2) ~= 2,
	error('Parameter ''lfp1'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
end
if size(lfp2.data,2) ~= 1 && size(lfp2.data,2) ~= 2,
	error('Parameter ''lfp2'' is not a vector or a Nx2 matrix (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
end

% Load SessionInfo
sessionInfo = bz_getSessionInfo(basepath,'noPrompts',true);

% Loading .mat file if exists
if ~isempty(foldername)
    if ~isempty(dir([basepath filesep sessionInfo.FileName, '.',foldername, '.*Coherogram.lfp.mat']))
        disp(['Coherogram ', foldername, ' already detected ! Loading file.'])
        file = dir([basepath filesep sessionInfo.FileName ,'.',foldername, '.*Coherogram.lfp.mat'])
        load(file.name);
        exist_file = true;
        %return
    end
else
    if ~isempty(dir([basepath filesep '.Coherogram.lfp.mat']))
        disp('Coherogram already detected! Loading file.')
        file = dir([basepath filesep '.Coherogram.lfp.mat']);
        load(file.name);
        exist_file = true;
        return
    end
end



% Determine LFP frequency
if isempty(frequency)
	if size(lfp1,2) == 2
		frequency = 1/median(diff(lfp1(:,1)));
	else
		frequency = f;
	end
end

% Determine step/overlap
if isempty(step)
	if isempty(overlap)
		overlap = window/2;
	end
else
	if isempty(overlap),
		overlap = window-step;
	elseif overlap ~= window-step,
		error('Incompatible ''step'' and ''overlap'' parameters (type ''help <a href="matlab:help MTCoherogram">MTCoherogram</a>'' for details).');
	end
end


% Compute and plot coherogram
parameters.Fs = frequency;
if ~isempty(range)
    parameters.fpass = range; 
end
parameters.tapers = tapers;
parameters.pad = pad;
% [coherogram,phase,~,S1,S2,t,f] = cohgramc(double(lfp1.data),double(lfp2.data),[window window-overlap],parameters);
[coherogram,phase,~,S1,S2,t,f] = cohgramc(double(lfp1.data),double(lfp2.data),[1 0.01],parameters);

% t = t'+lfp1(1,1);
f = f';
coherogram = coherogram';
% coherogram = permute(coherogram,[2 1 3]);  % Previous code by Gabrielle Girardeau, keep it around just in case
phase = phase';


if ~isempty('*.MergePoints.events.mat')
    file = dir('*.MergePoints.events.mat');
    load(file.name)
end

if ~isempty('*.Tracking.Behavior.mat')
    file = dir('*.Tracking.Behavior.mat');
    load(file.name);
end
    
if strcmpi(show,'on')
    disp('Showing Figure for Whole Recording')
    % Figure for whole recording
    minS1 = min(10*log10(mean(S1)));
    minS2 = min(10*log10(mean(S2)));
    minS = min(minS1,minS2);
    maxS1 = max(10*log10(mean(S1)));
    maxS2 = max(10*log10(mean(S2)));
    maxS = max(maxS1,maxS2);
    % Coherence and Power
    figure,
    set(gcf,'Position',get(0,'ScreenSize'))
    subplot(1,4,1), plot(f,mean(coherogram,2),'LineWidth',1), ylabel('Coherence(r)'), xlabel('Frequency(f)'), title(['Coherence Ch: ', num2str(lfp1.channels), ' and Ch: ', num2str(lfp2.channels)]), axis([0 200 0 1])
    subplot(1,4,2), plot(f,mean(phase,2),'LineWidth',1), ylabel('Coherence in frequency (r)'), xlabel('Frequency(f)'), title(['Coherence in Frequency Ch: ', num2str(lfp1.channels), ' and Ch: ', num2str(lfp2.channels)]), axis([0 200 -pi pi])
    subplot(1,4,3), plot(f,10*log10(mean(S1)), 'LineWidth', 1), ylabel('10*log10(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(lfp1.channels)]), ylim([minS-1 maxS+1])
    subplot(1,4,4), plot(f,10*log10(mean(S2)), 'LineWidth', 1), ylabel('10*log10(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(lfp2.channels)]), ylim([minS-1 maxS+1])
    
    if saveFig
      if ~isempty(foldername)
        saveas(gcf,['lfpAnalysisFigures\Power_WholeRecording_',foldername,'.png']);
      else    
        saveas(gcf,'lfpAnalysisFigures\Power_WholeRecording.png');
      end
    end

    figure, 
    set(gcf,'Position',get(0,'ScreenSize'))
    subplot(5,1,1),imagesc(t,f,10*log10(S1*10e12)'),view(0,-90),colorbar,title(['Power Ch: ', num2str(lfp1.channels)]), colormap(jet)
    subplot(5,1,2);imagesc(t,f,10*log10(S2*10e12)'), view(0,-90),colorbar,title(['Power Ch: ',num2str(lfp2.channels)]), colormap(jet)
    subplot(5,1,3);imagesc(t,f,coherogram),view(0,-90),colorbar, title(['Coherence Ch:  ', num2str(lfp1.channels), 'and Ch: ', num2str(lfp2.channels)])
    subplot(5,1,4);imagesc(t,f,phase),view(0,-90),colorbar, title(['Phase Coherence:  ', num2str(lfp1.channels), 'and Ch: ', num2str(lfp2.channels)])
    if ismember(foldername,tracking.folders)
        cd(foldername)
        if ~isempty('*Tracking.Behavior.mat')
            file = dir('*Tracking.Behavior.mat');
            load(file.name);
        end
        subplot(5,1,5);plot(tracking.timestamps,tracking.position.speed),title('Speed vs time'),colorbar,axis([0 max(tracking.timestamps) 0  max(tracking.position.speed)+0.01]); 
        cd(basepath)
    end
    
    if saveFig
      if ~isempty(foldername)
        saveas(gcf,['lfpAnalysisFigures\Periodogram_WholeRecording_',foldername,'.png']);
      else    
        saveas(gcf,'lfpAnalysisFigures\Periodogram_WholeRecording.png');
      end
    end
end


if filtering
    S1_TD_thres = 3;
    S2_TD_thres = 3;
    
    if ~isempty('*.MergePoints.events.mat')
        file = dir('*MergePoints.events.mat');
        load(file.name)
    end

    if ~isempty('*Tracking.Behavior.mat')
        file = dir('*Tracking.Behavior.mat');
        load(file.name);
    end

        
    FreqBands = [1 3; 4 12];
    nBands = size(FreqBands,1);
    
    % Let's see if exists tracking variable
    
    if ismember(foldername,tracking.folders)
        cd(foldername)
        if ~isempty('*Tracking.Behavior.mat')
            file = dir('*Tracking.Behavior.mat');
            load(file.name)
        end
        v = tracking.position.speed;
        Vi = 0.05;
        Vs = 0.50;
        indS1 = find(v<Vi);
        indS2 = find(v>Vs);
        mat = v;
        mat=mat';
        mat2 = mat;
        mat2(indS1) = 0;
        mat2(indS2) = 0;
        mat2(isnan(mat2)) = 0;
        
        % Finding epochs
        vvv = bwlabel(squeeze(mat2));
        MINPFSIZE = tracking.samplingRate * 2;
        [n2,bin2] = histc(tracking.timestamps,t);
        
        for jj=1:max(vvv(:))
            ttt = sum(vvv(:) == jj);
            [rb cb] = find(bwlabel(vvv) == jj);
            
            if ttt > MINPFSIZE
                C_aux = coherogram';
                phi_aux = phase';
                S1_aux = S1;
                S2_aux = S2;
                t_aux = t;
                cb_aux = cb;
                
                C_epoch = C_aux(bin2(cb(1)):bin2(cb(end)),:);
                phi_epoch = phi_aux(bin2(cb(1)):bin2(cb(end)),:);
                S1_epoch = S1_aux(bin2(cb(1)):bin2(cb(end)),:);
                S2_epoch = S2_aux(bin2(cb(1)):bin2(cb(end)),:);
                t_epoch = t_aux(bin2(cb(1)):bin2(cb(end)));
                cb_epoch = cb_aux;
                
                coherency_epochs{jj} = C_epoch;
                phi_epochs{jj} = phi_epoch;
                S1_epochs{jj} = S1_epoch;
                S2_epochs{jj} = S2_epoch;
                t_epochs{jj} = t_epoch;
                cb_epochs{jj} = cb_epoch;
            end
        end
        
        C_spd = cell2mat(coherency_epochs');
        phi_spd = cell2mat(phi_epochs');
        S1_spd = cell2mat(S1_epochs');
        S2_spd = cell2mat(S2_epochs');
        t_spd = cell2mat(t_epochs);
        cb_spd = cell2mat(cb_epochs);
        
        
        matPow_S1_spd = S1_spd';
        matPow_S2_spd = S2_spd';
        matPow_S1_spd_mean = mean(matPow_S1_spd,2);
        matPow_S2_spd_mean = mean(matPow_S2_spd,2);
        
        coherence_spd = C_spd';
        coherence_spd_mean = mean(coherence_spd,2);
        
        phi_spd = phi_spd';
        phi_spd_mean = mean(phi_spd,2);
        
        CV = 1;
        
        for it=1:length(t_spd)
            for k=1:nBands
                iDx = f >= FreqBands(k,1) & f < FreqBands(k,2);
                
                t_Pow_S1_spd = matPow_S1_spd(:,it);
                t_Pow_S2_spd = matPow_S2_spd(:,it);
                c_spd = coherence_spd(:,it);
                phase_spd = phi_spd(:,it);
                
                Eband_S1_spd(k,it) = CV*trapz(f(iDx),t_Pow_S1_spd(iDx));
                Eband_S2_spd(k,it) = CV*trapz(f(iDx),t_Pow_S2_spd(iDx));
                
                CoherenceBand_spd(k,it) = mean(c_spd(iDx));
                CoherenceFrequencyBand_spd(k,it) = mean(phase_spd(iDx));
                
            end
        end
               
        MeanPowerSpecS1_spd = mean(Eband_S1_spd,2)';
        MeanPowerSpecS2_spd = mean(Eband_S2_spd,2)';
        
        MeanCoherenceSpec_spd = mean(CoherenceBand_spd,2)';
        MeanCoherenceFrequency_spd = mean(CoherenceFrequencyBand_spd,2)';
        
        RatTD_S1_spd = Eband_S1_spd(2,1:end) ./ Eband_S1_spd(1,1:end);
        ind_td_S1_spd = find(RatTD_S1_spd >= S1_TD_thres);
        RatTD_S2_spd = Eband_S2_spd(2,1:end) ./ Eband_S2_spd(1,1:end);
        ind_td_S2_spd = find(RatTD_S2_spd >= S2_TD_thres);
        Crit_spd = [RatTD_S1_spd;RatTD_S2_spd];
        
        indx=find(Crit_spd(1,:)>=S1_TD_thres & Crit_spd(2,:)>=S2_TD_thres);
        [r_spd p_spd]=corrcoef(RatTD_S1_spd,RatTD_S2_spd);
        MeanPowerSpecRatS1_spd = mean(Eband_S1_spd(:,indx),2)';
        MeanPowerSpecRatS2_spd = mean(Eband_S2_spd(:,indx),2)';
        MeanCoherenceRat_spd = mean(CoherenceBand_spd(:,indx),2)';
        MeanCoherenceFrequencySpecRat_spd = mean(CoherenceFrequencyBand_spd(:,indx),2)';
        
        matPow_S1_spd_rat = matPow_S1_spd(:,indx);
        matPow_S1_spd_rat_mean = mean(matPow_S1_spd_rat,2);
        matPow_S2_spd_rat = matPow_S2_spd(:,indx);
        matPow_S2_spd_rat_mean = mean(matPow_S2_spd_rat,2);
        
        S1_spd_rat = S1_spd(indx,:);
        S2_spd_rat = S2_spd(indx,:);
        C_spd_rat = C_spd(indx,:);
        phi_spd = phi_spd';
        phi_spd_rat = phi_spd(indx,:);
        
        t_spd_rat = t_spd(indx);
        f_spd_rat = f;
        
        cd(basepath)
        if strcmpi(show,'on')
            figure,
            set(gcf,'Position',get(0,'ScreenSize'))
            subplot(1,4,1), plot(f,mean(C_spd_rat),'LineWidth',1),ylabel('coherence (r)'),xlabel ('Frequency (f)'), title(['Coherence Ch: ', num2str(lfp1.channels), 'and Ch: ', num2str(lfp2.channels),' filtered']),axis([0 200 0 1])
            subplot(1,4,2),plot(f,mean(phi_spd_rat),'LineWidth',1), ylabel(' Coherence in frequency(r)'), xlabel('Frequency(f)'), title(['Coherence Ch: ', num2str(lfp1.channels), 'and Ch: ', num2str(lfp2.channels), 'filtered']), axis([0 200 -pi pi])
            subplot(1,4,3),plot(f,10*log10(mean(S1_spd_rat)),'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(lfp1.channels), 'filtered'])
            subplot(1,4,4),plot(f,10*log10(mean(S2_spd_rat)),'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(lfp2.channels), 'filtered'])
            
            if saveFig
                if ~isempty(foldername)
                    saveas(gcf,['lfpAnalysisFigures\Power_',foldername,'_filtered','.png']);
                else    
                    saveas(gcf,'lfpAnalysisFigures\Power_filtered.png');
                end
            end
            
            
            figure,
            set(gcf,'Position',get(0,'ScreenSize'))
            subplot(5,1,1),imagesc(t_spd_rat,f,10*log10(S1_spd_rat*10e12)'),view(0,-90),colorbar,title(['Power Ch: ', num2str(lfp1.channels)]), colormap(jet)
            subplot(5,1,2);imagesc(t_spd_rat,f,10*log10(S2_spd_rat*10e12)'), view(0,-90),colorbar,title(['Power Ch: ', num2str(lfp1.channels)]), colormap(jet)
            subplot(5,1,3);imagesc(t_spd_rat,f,C_spd_rat'),view(0,-90),colorbar, title(['Coherence Ch:  ', num2str(lfp1.channels), 'and Ch: ', num2str(lfp2.channels)])
            subplot(5,1,4);imagesc(t_spd_rat,f,phi_spd_rat'),view(0,-90),colorbar, title(['Phase Coherence Ch:  ', num2str(lfp1.channels), 'and Ch: ', num2str(lfp2.channels)])
            subplot(5,1,5);plot(tracking.timestamps,v),title('Speed vs time'),colorbar, axis([0 max(tracking.timestamps) 0  max(v)+0.01]); hold on;
            subplot(5,1,5); plot(tracking.timestamps(cb_spd),v(cb_spd),'r'), title('Speed vs time epoch'),colorbar, axis([0 max(tracking.timestamps) 0  max(v)+0.01]);
            
            if saveFig
                if ~isempty(foldername)
                    saveas(gcf,['lfpAnalysisFigures\Periodogram_',foldername,'_filtered','.png']);
                else    
                    saveas(gcf,'lfpAnalysisFigures\Periodogram_filtered.png');
                end
            end
        end
    end
    cd(basepath)
end
%% Saving matlab struct

if ~isempty('*Tracking.Behavior.mat')
    file = dir('*Tracking.Behavior.mat');
    load(file.name)
end

if ismember(foldername,tracking.folders)
    c = C_spd_rat';
    coherogram = [];
    coherogram.coherogram = c;
    coherogram.phase = phi_spd_rat';
    coherogram.t = t_spd_rat;
    coherogram.f = f;
    coherogram.S1 = S1_spd_rat;
    coherogram.S2 = S2_spd_rat;
    coherogram.ch1 = lfp1.channels;
    coherogram.ch2 = lfp2.channels;
    if ~isempty(foldername)
        coherogram.foldername = foldername;
    end
else
    c = coherogram;
    coherogram = [];
    coherogram.coherogram = c;
    coherogram.phase = phase;
    coherogram.t = t;
    coherogram.f = f;
    coherogram.S1 = S1;
    coherogram.S2 = S2;
    coherogram.ch1 = lfp1.channels;
    coherogram.ch2 = lfp2.channels;
    if ~isempty(foldername)
        coherogram.foldername = foldername;
    end
end


if saveMat
    if ~isempty(foldername)
        try
            save([basepath filesep sessionInfo.FileName, '.', foldername, '.Coherogram.lfp.mat'], 'coherogram');
        catch
            save([basepath filesep sessionInfo.FileName, '.', foldername, '.Coherogram.lfp.mat'], 'coherogram', '-v7.3');
        end
    else
        try
            save([basepath filesep sessionInfo.FileName '.Coherogram.lfp.mat'], 'coherogram');
        catch
            save([basepath filesep sessionInfo.FileName '.Coherogram.lfp.mat'], 'coherogram','-v7.3');
        end
    end
end
   
end


