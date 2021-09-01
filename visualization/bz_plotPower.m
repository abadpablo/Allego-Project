function  bz_plotPower(varargin)
% bz_plotPower(basepath) plots the power spectrum and spectogram per shank
%
%
%
%
%
%
%
%
%

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'ignoreShanks',[],@isnumeric);
addParameter(p,'fpass',[1 200],@isnumeric);
addParameter(p,'tapers',[3 5], @isnumeric);
addParameter(p,'pad',1,@isnumeric);
addParameter(p,'savefig',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
ignoreShanks = p.Results.ignoreShanks;
fpass = p.Results.fpass;
tapers = p.Results.tapers;
pad = p.Results.pad;
savefig = p.Results.savefig;


%%
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
shanks = sessionInfo.AnatGrps;
shanks(ignoreShanks) = [];
lfpSampleRate = sessionInfo.lfpSampleRate;

params.Fs = lfpSampleRate;
params.fpass = fpass;
params.tapers = tapers; 
params.pad = pad;

for jj=1:length(shanks)
    lfp = bz_GetLFP(shanks(jj).Channels,'noPrompts',true);
    [S,t,f] = mtspecgramc(single(lfp.data),[2 1], params);
    S = log10(S); % in Db
    count1=1;
    count2=2;
    figure,
    set(gcf,'Position',[100 100 1400 600])
    
    for ii=1:length(lfp.channels)
        
        s = S(:,:,ii);
        s_det= bsxfun(@minus,s,polyval(polyfit(f,mean(s,1),2),f)); % detrending 
        subplot(2,8,count1)
        imagesc(t,f,s_det',[-1.5 1.5]);
        set(gca,'XTick',[]); ylabel('Freqs');
        title(['Channel : ' num2str(lfp.channels(ii))])
        count1=count1+2;
        subplot(2,8,count2)
        plot(mean(s,1),f,'k');
        set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
        count2=count2+2;
        
        if savefig
            saveas(gcf,['SummaryFigures\spectrogramAllSession_shank',num2str(ii),'.png']);
        end
    end
    
end

