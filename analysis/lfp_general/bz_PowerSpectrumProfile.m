function [powerProfile] = bz_PowerSpectrumProfile(frange,varargin)
% power distribution per channels
%
%   INPUTS
%   frange      (Hz, If empty take by default 6 12, theta)
%
%   (optional)
%    lfp        a buzcode-formatted lfp structure (use bz_GetLFP)
%               needs fields: lfp.data, lfp.timestamps, lfp.samplingRate.
%               If empty or no exist, look for lfp in basePath folder
%    winsize    size of the silding time window (s, default 2)
%    dt         sliding time interval (s, default 1)
%    channels   subset of channels to calculate PowerSpectrumSlope
%               (default: all)
%    showfig    true/false - show a summary figure of the results
%               (default:false)
%    saveMat    put your basePath here to save/load
%               baseName.PowerSpectrumProfile_'frange'.lfp.mat  (default: true)
%   forceDetect (default false)
%
%   OUTPUTS
%   powerProfile
%       .mean       log10-transformed mean amplitude of the spectrogram
%       .std
%       .median
%       .channels
%       .channels_shank
%       .frange
%       .f
%       .ic95
%
% MV-BuzsakiLab 2019
% Edited by Peter Petersen

% TODO
% Handle bad channels
% Remove sessionInfo dependencies


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Parsing inputs
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
p = inputParser;
addParameter(p,'winSize',4,@isscalar)
addParameter(p,'dt',2, @isscalar)
addParameter(p,'channels','all')
addParameter(p,'showfig',true,@islogical)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'lfp',[])
addParameter(p,'forceDetect',false,@islogical)
addParameter(p,'useParfor',true,@islogical)
addParameter(p,'timestampsSubSession',[],@isnumeric);
addParameter(p,'foldername',[],@isstr);
addParameter(p,'force',false,@islogical);
addParameter(p,'tapers',[3 5], @isnumeric);
addParameter(p,'pad',1,@isnumeric);

parse(p,varargin{:})
showfig = p.Results.showfig;
saveMat = p.Results.saveMat;
channels = p.Results.channels;
winSize = p.Results.winSize;
lfp = p.Results.lfp;
dt = p.Results.dt;
forceDetect = p.Results.forceDetect;
useParfor = p.Results.useParfor;
timestampsSubSession = p.Results.timestampsSubSession;
foldername = p.Results.foldername;
force = p.Results.force;
tapers = p.Results.tapers;
pad = p.Results.pad;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Resolving inputs   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
try [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
    
    if ~isempty(foldername)
        filename = dir(['*.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.SubSession.channelinfo.mat']);
        if ~isempty(filename) && ~forceDetect
            if length(filename) >1
                for i=1:length(filename)
                    name{i} = filename(i).name;
                    name{i} = strsplit(name{i},'.');
                    name{i} = name{i}{2};
                end
                if any(strcmpi(name,foldername))
                        disp(['Power Spectrum Profile_',num2str(frange(2)),'_',num2str(frange(2)),' already detected ! Loading file: ', filename(strcmpi(name,foldername)).name])
                        load(filename(strcmpi(name,foldername)).name)
                        return
                end
            else
                name = filename.name;
                name = strsplit(name,'.');
                name = name{2};
                
                if strcmpi(name, foldername)
                    disp(['Power Spectrum Profile_',num2str(frange(2)),'_',num2str(frange(2)),' already detected ! Loading file: ', filename.name])
                    load(filename.name);
                    return
                end
            end
        end
        
    elseif (exist([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.channelinfo.mat'],'file') || ...
            exist([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.lfp.mat'],'file')) ...
            && ~forceDetect
        disp(['Power spectrum profile already calculated for ', sessionInfo.FileName, ' in the range: ', num2str(frange(1)),' - ',num2str(frange(2)), 'Hz. Loading file.']);
        try load([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.channelinfo.mat']);
        catch
            load([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.lfp.mat']);
        end
        return
    end
    
    if ischar('channels') && strcmpi(channels,'all')
        channels = [1:sessionInfo.nChannels]-1;
    end
    
catch disp('No session info file!!');
    if ischar('channels') && strcmpi(channels,'all')
        channels = lfp.channels;
    end
    keyboard;
    sessionInfo.rates.lfp = lfp.samplingRate;
    useParfor = false;
    sessionInfo.FileName = date;
    sessionInfo.AnatGrps(1).Channels = lfp.channels; 
end
    
if ~exist('frange') || isempty(frange)
    frange = [6 12];
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Calculate spectrogram per channel
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
params.Fs = sessionInfo.rates.lfp; params.fpass = frange; 
params.tapers = tapers; params.pad = pad;
tic
powerProfileMean = [];
powerProfileStd = [];
powerProfileIc95 = [];
powerProfileMedian = [];
powerProfileChannels = [];
disp('Calculating spectrograms channelwise')
if useParfor
    parfor (ii = 1:length(channels),18)
        fprintf('Channel %3.i/%3.i, ',ii, length(channels));
        if ~isempty(timestampsSubSession)
            lfp = bz_GetLFP(channels(ii),'restrict',timestampsSubSession,'noPrompts', true);
            [S,t,f] = mtspecgramc_fast(single(lfp.data),[4 2],params);
            S = 10 * log10(S);
            powerProfileMean(ii) = mean(mean(S,2));
            powerProfileStd(ii) = std(mean(S,2));
            powerProfileIc95(ii) = 1.96 * std(mean(S,2))/sqrt(length(mean(S,2)));
            powerProfileMedian(ii) = median(median(S,2));
            powerProfileChannels(ii) = channels(ii);
        else
            lfp = bz_GetLFP(channels(ii),'noPrompts', true);
            [S,t,f] = mtspecgramc_fast(single(lfp.data),[4 2],params);
            S = 10 * log10(S);
            powerProfileMean(ii) = mean(mean(S,2));
            powerProfileStd(ii) = std(mean(S,2));
            powerProfileIc95(ii) = 1.96 * std(mean(S,2))/sqrt(length(mean(S,2)));
            powerProfileMedian(ii) = median(median(S,2));
            powerProfileChannels(ii) = channels(ii);
        end
    end
else
    
    if ~isempty(timestampsSubSession)
        
        if isempty(lfp)
            lfp = bz_GetLFP('all','restrict',timestampsSubSession,'noPrompts',true);
        end
        for ii=1:length(channels)
            fprintf('Channel %3.i/%3.i ,\n',ii, length(channels));
            [S,t,f] = mtspecgramc_fast(single(lfp.data(:,channels(ii)+1)),[4 2],params);
            S = 10 * log10(S);
            powerProfileMean(ii) = mean(mean(S,2));
            powerProfileStd(ii) = std(mean(S,2));
            powerProfileIc95(ii) = 1.96 * std(mean(S,2))/sqrt(length(mean(S,2)));
            powerProfileMedian(ii) = median(median(S,2));
            powerProfileChannels(ii) = channels(ii);
        end
    else
        if isempty(lfp)
            lfp = bz_GetLFP('all','noPrompts', true);
        end
        for ii = 1:length(channels)
            fprintf('Channel %3.i/%3.i ,\n',ii, length(channels));
            [S,t,f] = mtspecgramc_fast(single(lfp.data(:,channels(ii)+1)),[4 2],params);
            S = 10 * log10(S);
            powerProfileMean(ii) = mean(mean(S,2));
            powerProfileStd(ii) = std(mean(S,2));
            powerProfileIc95(ii) = 1.96 * std(mean(S,2))/sqrt(length(mean(S,2)));
            powerProfileMedian(ii) = median(median(S,2));
            powerProfileChannels(ii) = channels(ii);
        end
    end
end
toc
powerProfile.mean = powerProfileMean;
powerProfile.std = powerProfileStd;
powerProfile.ic95 = powerProfileIc95;
powerProfile.median = powerProfileMedian;
powerProfile.channels = powerProfileChannels;


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Saving the result to basename.PowerSpectrumProfile_frange.channelinfo.mat
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
powerProfile.processinginfo.function = 'bz_PowerSpectrumProfile';
powerProfile.processinginfo.date = now;
powerProfile.processinginfo.params.winSize = winSize;
powerProfile.processinginfo.params.dt = dt;
powerProfile.processinginfo.params.frange = frange;
if saveMat
    if ~isempty(timestampsSubSession)
        save([strcat(sessionInfo.FileName,'.',foldername), '.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.SubSession.channelinfo.mat'], 'powerProfile');
    else
        save([sessionInfo.FileName,'.PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.channelinfo.mat'],'powerProfile');

    end
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Plotting
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if showfig
    figure('Position', [20 20 1800 900],'Name',sessionInfo.FileName)
    cmap = jet(size(sessionInfo.AnatGrps,2));
    plt1 = [];
    for ii = 1:size(sessionInfo.AnatGrps,2)
        [Lia] = ismember(sessionInfo.AnatGrps(ii).Channels, channels);
        nC = 1:length(sessionInfo.AnatGrps(ii).Channels);
        nC = nC(Lia);
        hold on
        dev1 = powerProfile.mean(sessionInfo.AnatGrps(ii).Channels(Lia)+1) - powerProfile.ic95(sessionInfo.AnatGrps(ii).Channels(Lia)+1);
        dev2 = powerProfile.mean(sessionInfo.AnatGrps(ii).Channels(Lia)+1) + powerProfile.ic95(sessionInfo.AnatGrps(ii).Channels(Lia)+1);
        
        hold on
        fill([dev1 flip(dev2)],[nC flip(nC)],cmap(ii,:),'FaceAlpha',.2,'EdgeColor','none')
        try plt1(ii) = plot(powerProfile.mean(sessionInfo.AnatGrps(ii).Channels(Lia)+1),nC(Lia),'color',cmap(ii,:)); end
    end
    
    ax=axis; axis tight; xlim(ax(1:2));
    ylabel('Channels'); xlabel('Power'); title(strcat('Freq range:',num2str(frange),'Hz'),'FontWeight','normal');
    set(gca,'YDir','reverse');
    legend(plt1,{num2str([1:size(sessionInfo.AnatGrps,2)]')})
    if ~exist('SummaryFigures','dir')
        mkdir('SummaryFigures')
    end
    
    if ~isempty(timestampsSubSession)
        saveas(gcf,['SummaryFigures\PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'_',foldername,'.png']);
    else
        saveas(gcf,['SummaryFigures\PowerSpectrumProfile_',num2str(frange(1)),'_',num2str(frange(2)),'.png']);
    end
end

end