function [coherogram] = bz_MTCoherogram(lfp1,lfp2,varargin)

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
%         return
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

if ~exist_file

    % Compute and plot coherogram
    parameters.Fs = frequency;
    if ~isempty(range)
        parameters.fpass = range; 
    end
    parameters.tapers = tapers;
    parameters.pad = pad;
    [coherogram,phase,~,S1,S2,t,f] = cohgramc(double(lfp1.data),double(lfp2.data),[window window-overlap],parameters);
%     [coherogram,phase,~,S1,S2,t,f] = cohgramc(double(lfp1.data),double(lfp2.data),[0.1 0.01],parameters);

    % t = t'+lfp1(1,1);
    f = f';
    coherogram = coherogram';
    % coherogram = permute(coherogram,[2 1 3]);  % Previous code by Gabrielle Girardeau, keep it around just in case
    phase = phase';


    minS1 = min(10*log10(mean(S1)));
    minS2 = min(10*log10(mean(S2)));
    minS = min(minS1,minS2);
    maxS1 = max(10*log10(mean(S1)));
    maxS2 = max(10*log10(mean(S2)));
    maxS = max(maxS1,maxS2);

    if strcmpi(show,'on')
      figure
      set(gcf,'Position',get(0,'ScreenSize'))
      hold on;
      subplot(2,1,1);
      PlotColorMap(coherogram,'x',t,'y',f,'cutoffs',cutoffs,'newfig','off');
      xlabel('Time (s)');
      ylabel('Frequency (Hz)');
      title(['Coherogram Amplitude Ch ', num2str(lfp1.channels),' vs Ch ',num2str(lfp2.channels)]);
      colorbar
      subplot(2,1,2);
      PlotColorMap(phase,'x',t,'y',f,'cutoffs',[-pi pi],'newfig','off');
      xlabel('Time (s)');
      ylabel('Frequency (Hz)');
      title(['Coherogram Phase Ch ', num2str(lfp1.channels), ' vs Ch ', num2str(lfp2.channels)]);
      colorbar

      if saveFig
          if ~isempty(foldername)
            saveas(gcf,['lfpAnalysisFigures\Coherogram_',foldername,'.png']);
          else    
            saveas(gcf,'lfpAnalysisFigures\Coherogram.png');
          end
      end

      figure,
      set(gcf,'Position',get(0,'ScreenSize'))
      subplot(1,4,1)
      plot(f,mean(coherogram'))
      xlabel('Frequency (f)');
      ylabel('Coherence (r)')
      title(['Coherence Ch ', num2str(lfp1.channels), ' vs Ch ', num2str(lfp2.channels)])
      xlim(range)
      ylim([min(mean(coherogram')) max(mean(coherogram'))])
      subplot(1,4,2)
      plot(f,mean(phase'))
      xlabel('Frequency (f)');
      ylabel('Coherence Phase')
      title(['Coherence Phase Ch ', num2str(lfp1.channels), ' vs Ch ', num2str(lfp2.channels)])
      xlim(range)
      ylim([min(mean(phase')) max(mean(phase'))])
      subplot(1,4,3)
      plot(f,10*log10(mean(S1)));
      xlabel('Frequency (f) ');
      ylabel('10*log10');
      title(['Power Spectrum Ch ', num2str(lfp1.channels)]);
      xlim([range])
      ylim([minS maxS])
      subplot(1,4,4)
      plot(f,10*log10(mean(S2)));
      xlabel('Frequency (f) ');
      ylabel('10*log10');
      title(['Power Spectrum Ch ', num2str(lfp2.channels)]);
      xlim([range])
      ylim([minS maxS])

      if saveFig
          if ~isempty(foldername)
            saveas(gcf,['lfpAnalysisFigures\Coherence_',foldername,'.png']);
          else    
            saveas(gcf,'lfpAnalysisFigures\Coherence.png');
          end
      end

    end

    %% Saving matlab struct
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
else
    minS1 = min(10*log10(mean(coherogram.S1)));
    minS2 = min(10*log10(mean(coherogram.S2)));
    minS = min(minS1,minS2);
    maxS1 = max(10*log10(mean(coherogram.S1)));
    maxS2 = max(10*log10(mean(coherogram.S2)));
    maxS = max(maxS1,maxS2);
    
    if strcmpi(show,'on')
      if ~isempty(foldername)
        figure('Name',coherogram.foldername),
      else
        figure('Name',basepath)
      end
      set(gcf,'Position',get(0,'ScreenSize'))
      hold on;
      subplot(2,1,1);
      PlotColorMap(coherogram.coherogram,'x',coherogram.t,'y',coherogram.f,'cutoffs',cutoffs,'newfig','off');
      xlabel('Time (s)');
      ylabel('Frequency (Hz)');
      title(['Coherogram Amplitude Ch ', num2str(lfp1.channels),' vs Ch ',num2str(lfp2.channels)]);
      colorbar
      subplot(2,1,2);
      PlotColorMap(coherogram.phase,'x',coherogram.t,'y',coherogram.f,'cutoffs',[-pi pi],'newfig','off');
      xlabel('Time (s)');
      ylabel('Frequency (Hz)');
      title(['Coherogram Phase Ch ', num2str(lfp1.channels), ' vs Ch ', num2str(lfp2.channels)]);
      colorbar

      if ~isempty(foldername)
        figure('Name',coherogram.foldername),
      else
        figure('Name',basepath)
      end
      set(gcf,'Position',get(0,'ScreenSize'))
      subplot(1,4,1)
      plot(coherogram.f,mean(coherogram.coherogram'))
      xlabel('Frequency (f)');
      ylabel('Coherence (r)')
      title(['Coherence Ch ', num2str(lfp1.channels), ' vs Ch ', num2str(lfp2.channels)])
      xlim(range)
      ylim([min(mean(coherogram.coherogram')) max(mean(coherogram.coherogram'))])
      subplot(1,4,2)
      plot(coherogram.f,mean(coherogram.phase'))
      xlabel('Frequency (f)');
      ylabel('Coherence Phase')
      title(['Coherence Phase Ch ', num2str(lfp1.channels), ' vs Ch ', num2str(lfp2.channels)])
      xlim(range)
      ylim([min(mean(coherogram.phase')) max(mean(coherogram.phase'))])
      subplot(1,4,3)
      plot(coherogram.f,10*log10(mean(coherogram.S1)));
      xlabel('Frequency (f) ');
      ylabel('10*log10');
      title(['Power Spectrum Ch ', num2str(lfp1.channels)]);
      xlim([range])
      ylim([minS maxS])
      subplot(1,4,4)
      plot(coherogram.f,10*log10(mean(coherogram.S2)));
      xlabel('Frequency (f) ');
      ylabel('10*log10');
      title(['Power Spectrum Ch ', num2str(lfp2.channels)]);
      xlim([range])
      ylim([minS maxS])

    end
    
    
    
end

