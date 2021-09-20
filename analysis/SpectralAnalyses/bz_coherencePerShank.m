function [coherence_Shanks] = bz_coherencePerShank(varargin)
% bz_coherencePerShank - Compute coherence between shanks (one channel per
%   shank)
%
%   USAGE
%   
%   [] = bz_coherencePerShank(varargin)
%   
%   INPUTS
%   <options> optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       
%      'lfp'       lfp struct from buzcode
%     
%      PARAMS
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
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.

%% Make sure Chronuz is installed
CheckChronux('cohgramc')

%% Default params

p = inputParser();
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'lfp',[],@bz_islfp);
addParameter(p,'timestampsSubSession',[],@isnumeric);
addParameter(p,'discardShanks',[],@isnumeric);
addParameter(p,'frequency',1250,@isnumeric);
addParameter(p,'range',[0 200],@isnumeric);
addParameter(p,'window',5,@isnumeric);
addParameter(p,'overlap',2.5,@isnumeric);
addParameter(p,'step',2.5,@isnumeric);
addParameter(p,'tapers',[3 5], @isnumeric);
addParameter(p,'pad',0,@isnumeric);
addParameter(p,'show','on',@isstr);
addParameter(p,'cutoffs',[0 1], @isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'exist_file',false,@islogical);
addParameter(p,'foldername',[],@isstr);

parse(p,varargin{:})
lfp = p.Results.lfp;
timestampsSubSession = p.Results.timestampsSubSession;
discardShanks = p.Results.discardShanks;
frequency = p.Results.frequency;
range = p.Results.range;
window = p.Results.window;
overlap = p.Results.overlap;
step = p.Results.step;
tapers = p.Results.tapers;
pad = p.Results.pad;
show = p.Results.show;
cutoffs = p.Results.cutoffs;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;
basepath = p.Results.basepath;
exist_file = p.Results.exist_file;
foldername = p.Results.foldername;

[sessionInfo] = bz_getSessionInfo(basepath,'noPrompts',true);
session = loadSession;
excludeChannels = [];

shanks = sessionInfo.AnatGrps;
shanks(discardShanks) = [];

if ~isempty(foldername)
    if ~isempty(dir([basepath filesep sessionInfo.FileName, '.',foldername, '.*Coherence_Shanks.lfp.mat']))
        disp(['Coherence per Shank', foldername, 'already detected ! Loading file.'])
        file = dir([basepath filesep sessionInfo.FileName ,'.',foldername, '.*Coherence_Shanks.lfp.mat'])
        load(file.name);
        exist_file = true;
        %return
    end
else
    if ~isempty(dir([basepath filesep '.Coherence_Shanks.lfp.mat']))
        disp('Coherence per Shank already detected! Loading file.')
        file = dir([basepath filesep '.Coherence_Shanks.lfp.mat']);
        load(file.name);
        exist_file = true;
        %return
    end
end

for ii=1:length(discardShanks)
    excludeChannels = session.extracellular.electrodeGroups.channels(excludeShanks(ii));
end

if exist('chanMap.mat','file')
    load('chanMap.mat','xcoords','ycoords','kcoords');
    xcoords = xcoords - min(xcoords); xcoords = xcoords/max(xcoords);
    ycoords = ycoords - min(ycoords); ycoords = ycoords/max(ycoords); 
    xcoords(excludeChannels) = [];
    ycoords(excludeChannels) = [];
    kcoords(excludeChannels) = [];
            
else
    xcoords = NaN;
    ycoords = NaN;
    kcoords = NaN;
end
% Let's pick the first channel for each shank
channel = 1;
for jj=1:size(shanks,2)
   ch(jj) = shanks(jj).Channels(channel);     
end

if ~exist_file
    % Define params for chronux functions
    params.Fs = frequency;
    if ~isempty(range)
        params.fpass = range;
    end
    params.tapers = tapers;
    params.pad = pad;
    % Chronux function
    for jj=1:length(ch)
        lfp1 = bz_GetLFP(ch(jj),'restrict',timestampsSubSession);
        for ii=1:length(ch)
            lfp2 = bz_GetLFP(ch(ii),'restrict',timestampsSubSession);
            [coherogram{jj}{ii},phase{jj}{ii},S12{jj}{ii},S1{jj}{ii},S2{jj}{ii},t,f] = cohgramc(double(lfp1.data),double(lfp2.data),[window window-overlap],params);
        end
    end

    % Plot figures
    count = 1;
    if strcmp(show,'on')

        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj=1:length(ch)
            for ii=1:length(ch)
                subplot(size(shanks,2),size(shanks,2),count);
                PlotColorMap(coherogram{jj}{ii}','x',t,'y',f','cutoffs',cutoffs,'newfig','off');
                count = count+1;
                if jj ~= length(ch) && ii == 1
                    ylabel('Trial');
                    set(gca,'XTick',[]);
                elseif jj == length(ch) && ii == 1
                    xlabel('Time (s)')
                    ylabel('Trial')
                elseif jj == length(ch) && ii ~= 1
                    xlabel('Time(s)')
                else
                    set(gca,'YTick',[],'XTick',[]);
                end

            end
        end

        if saveFig
           saveas(gcf,'lfpAnalysisFigures\Coherogram_per_Shanks.png') 
        end

        count = 1;
        figure,
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj=1:length(ch)
            for ii=1:length(ch)
                subplot(size(shanks,2),size(shanks,2),count);
                plot(f,mean(coherogram{jj}{ii}))
                count = count+1;
                 if jj ~= length(ch) && ii == 1
                    ylabel('Coherence (r)');
                    set(gca,'XTick',[]);
                elseif jj == length(ch) && ii == 1
                    xlabel('Frequency (f)')
                    ylabel('Coherence (r) ')
                elseif jj == length(ch) && ii ~= 1
                    xlabel('Frequency (f)')
                else
                    set(gca,'YTick',[],'XTick',[]);
                 end
                 ylim([0.4 1])

            end

        end

        if saveFig
            saveas(gcf,'lfpAnalysisFigures\Coherency_per_Shanks.png') 
        end
    end

    %% Save results in matlab struct
    coherence_Shanks = [];
    coherence_Shanks.coherogram = coherogram;
    coherence_Shanks.phase = phase;
    coherence_Shanks.S12 = S12;
    coherence_Shanks.S1 = S1;
    coherence_Shanks.S2 = S2;
    coherence_Shanks.t = t;
    coherence_Shanks.f = f;
    coherence_Shanks.ch = ch;
    
    if ~isempty(foldername)
        coherence_Shanks.foldername = foldername;
    end

    if saveMat
        if ~isempty(foldername)
           save([basepath filesep sessionInfo.FileName, '.', foldername, '.Coherence_Shanks.lfp.mat'],'coherence_Shanks','-v7.3');
        else
           save([basepath filesep sessionInfo.FileName, '.Coherence_Shanks.lfp.mat'],'coherence_Shanks','-v7.3');  
        end
    end
else
    % File already exist and it is loaded, so we only plot the data
% Plot figures
    count = 1;
    if strcmp(show,'on')
        if ~isempty(foldername)
            figure('Name',coherence_Shanks.foldername),
        else
            figure('Name',basepath)
        end
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj=1:length(ch)
            for ii=1:length(ch)
                subplot(size(shanks,2),size(shanks,2),count);
                PlotColorMap(coherence_Shanks.coherogram{jj}{ii}','x',coherence_Shanks.t,'y',coherence_Shanks.f','cutoffs',cutoffs,'newfig','off');
                count = count+1;
                if jj ~= length(ch) && ii == 1
                    ylabel('Trial');
                    set(gca,'XTick',[]);
                elseif jj == length(ch) && ii == 1
                    xlabel('Time (s)')
                    ylabel('Trial')
                elseif jj == length(ch) && ii ~= 1
                    xlabel('Time(s)')
                else
                    set(gca,'YTick',[],'XTick',[]);
                end

            end
        end
        
        count = 1;
        if ~isempty(foldername)
            figure('Name',coherence_Shanks.foldername),
        else
            figure('Name',basepath)
        end
        set(gcf,'Position',get(0,'ScreenSize'))
        for jj=1:length(ch)
            for ii=1:length(ch)
                subplot(size(shanks,2),size(shanks,2),count);
                plot(coherence_Shanks.f,mean(coherence_Shanks.coherogram{jj}{ii}))
                count = count+1;
                 if jj ~= length(ch) && ii == 1
                    ylabel('Coherence (r)');
                    set(gca,'XTick',[]);
                elseif jj == length(ch) && ii == 1
                    xlabel('Frequency (f)')
                    ylabel('Coherence (r) ')
                elseif jj == length(ch) && ii ~= 1
                    xlabel('Frequency (f)')
                else
                    set(gca,'YTick',[],'XTick',[]);
                 end
                 ylim([0.4 1])

            end

        end
            
            
    end
end



end

