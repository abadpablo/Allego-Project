function [pld,pld_gamma] = bz_SpikeLFP(varargin)

% Defaults
p = inputParser;
addRequired(p,'spikes',@bz_isCellInfo);
addRequired(p,'lfp',@bz_isLFP);
addRequired(p,'passband',@isnumeric);
addParameter(p,'intervals',[0 inf],@isnumeric)
addParameter(p,'samplingRate',1250,@isnumeric)
addParameter(p,'method','hilbert',@isstr)
addParameter(p,'plotting',true,@islogical)
addParameter(p,'numBins',180,@isnumeric)
addParameter(p,'powerThresh',2,@isnumeric)
addParameter(p,'saveMat',false,@islogical)
addParameter(p,'force',false,@islogical);
addParameter(p,'unit',[],@isnumeric);
addParameter(p,'nBins',19,@isnumeric);

parse(p,varargin{:})

spikes = p.Results.spikes;
lfp = p.Results.lfp;
passband = p.Results.passband;

intervals = p.Results.intervals; % interval(s) over which to calculate
samplingRate = p.Results.samplingRate; % sampling rate of continuous signal (LFP)
method = p.Results.method;
plotting = p.Results.plotting;
numBins = p.Results.numBins;
powerThresh = p.Results.powerThresh;
saveMat = p.Results.saveMat;
force = p.Results.force;
unit = p.Results.unit;
nBins = p.Results.nBins;

if ~isempty(unit)
    spikeTimes{1} = spikes.times{unit};
    spikes = rmfield(spikes,'times');
    spikes.times = spikeTimes;
end

BCC = fir1(400,[passband(1)/(samplingRate/2) passband(2)/(samplingRate/2)],'bandpass'); % theta filter convolution bands as in Monyer paper Neuron
th_all = filtfilt(BCC,[1],lfp.data) ;%theta band filtered
thetaEnv = abs(hilbert(th_all));  % envelop
Threstheta2=mean(thetaEnv)-std( thetaEnv ) ;   % threshold for peak amp detection

k=length(lfp.data);
EEG_ts_step=1/samplingRate;
tstep = lfp.timestamps;

count1=0;% for phase calculation
count2=0; % for epochs of 2 secs meeting theta/Delta ratio
count3=0; % cicles meeting criteria for analysis
count4=0; % counter for theta avergage cicle
count7=0; % count to use in Anal1

for kl=2:2:(length(tstep))/samplingRate-2
    t1=kl;                                             % Time 1 of epoch 
    t2=kl+2;                                            % Time 2 of epoch
    y = lfp.data(samplingRate*t1:samplingRate*t2);                               % Signal epoch 
    Eptime = tstep(samplingRate*t1:samplingRate*t2);                        % Epoch timed vector
    EptInt = EEG_ts_step*(1:length(Eptime)) + intervals(1);              % helps to find interval for degress CM position
    y = locdetrend(y',samplingRate,[0.5 0.2])';                % Data detrended  
    Args.Fs = samplingRate;
    [Eband Eband2 RatTD p  f] =MaKPower2(y,Args,0);     % Power spectra calculation for theta/delta ratio
    %% Theta Filter
    BC = fir1(400,[passband(1)/(samplingRate/2) passband(2)/(samplingRate/2)],'bandpass');         % theta filter convolution bands 
    th=filtfilt(BC,[1],y) ;%theta band filtered
    thetaEnv_ep = abs(hilbert(th));  % envelop
    Threstheta=mean(thetaEnv)+std( thetaEnv )*5 ;   % threshold  for noise detection    
    %% Gamma
    B = fir1(400,[30/(samplingRate/2) 65/(samplingRate/2)],'bandpass');        %gamma filter convolution as Buzsaki 2003
    gam= filtfilt(B,[1],y); %gamma filtered
    gammEnv_ep = abs(hilbert(gam));
 
    if (RatTD >6) &&  (max(thetaEnv_ep)< Threstheta)  % Only Epochs of strong theta are included for analysis eliminating high amplitude artifacts

        count2=count2+1;
    %% theta trhough detection
    [Negmmth Negppth] = findpeaks((-th),'minpeakdistance',round(samplingRate*0.07));
    ThAmp=mean(abs(th)); 
    Stdth=std((abs(th))); % Mean values of amplituded 
    AmpThr=ThAmp-Stdth; % Theta Amplitude Signal thresholded
    
    %% Gamma thorugh detection
    [Negmmth_g Negppth_g] = findpeaks((-gam),'minpeakdistance',round(samplingRate*0.0154));
    GamAmp=mean(abs(gam)); 
    Stdgam=std((abs(gam))); % Mean values of amplituded 
    Ampgamr=GamAmp-Stdgam; % Theta Amplitude Signal thresholded
    %% Theta phase detection and degress epoch detection

    X = hilbert(th) ; % Hilbert transformed
    phi = angle(X);   % angle of the different signal samples
    gradthe=rad2deg(phi)+180; % corrected to 0-360º
    [mmGradth ppGradth] = findpeaks(gradthe,'minpeakdistance',round(samplingRate*0.07));     % finding neg peaks so transitons between cicles
    %% Gamma phase detection and degress epoch detection

    X_g = hilbert(gam) ; % Hilbert transformed
    phi_g = angle(X_g);   % angle of the different signal samples
    gradthe_g=rad2deg(phi_g)+180; % corrected to 0-360º
    [mmGradth_g ppGradth_g] = findpeaks(gradthe_g,'minpeakdistance',round(samplingRate*0.0154));     % finding neg peaks so transitons between cicles

    
    
    %% PLoting data
    plott=0;
    Anal1=1;
        if Anal1==1
            count7=count7+1;
            if plott==1
            [S2,t,f2] = nbt_TF_plot2(y,0.1,1,50,1,0,Args);                 % wavelets
            S2=sqrt(abs(S2));S2=flipud(S2);
            tt2=t;
            figure;                                                              % plot signal and components
            subplot(3,1,1);hold; plot( Eptime,y,'k','linewidth',1.50 );ylabel('Amp (micv)') % plot Raw Signal
            axis([Eptime(1)  Eptime(end) -0.8 0.8]);
            plot( Eptime,gradthe/1000,'c','linewidth',1.50 );                    % plot Degrees Signal 
            text( Eptime(ppGradth),mmGradth/1000,'PeakGr');                      % texting Neg Degrees Peaks  
            plot( Eptime,th,'linewidth',1.25 );                                  % Plot theta filtered signal
            text( Eptime(Negppth),-Negmmth,'NegT');                              % Texting theta valleys
            plot( Eptime,gam,'r','linewidth',1.25 );                             % Plot gamma signal  
            subplot(3,1,2); imagesc(t,f2,S2); colorbar('SouthOutside' );axis xy % Plots waveltes
            subplot(3,1,3);plot( Eptime,y);hold;plot( Eptime,th,'r');plot( Eptime,gradthe/1000,'c' );   %plot Degrees Signal 
            end
            for ci=1:length(spikes.times)
                spikeTimes = spikes.times{ci};
                status = InIntervals(spikeTimes,[Eptime(1) Eptime(end)]);
                spikeTimes = spikeTimes(status);
                [tbin bbin2] =histc(spikeTimes,Eptime); % trying to bin spike times to Eptime 
                intT=find(bbin2>0);
                ytext=zeros(1,length(intT));
                phase_locking = gradthe(bbin2);                                               % Spike Phase Locking                                            
                phase_locking2 = gradthe(Negppth); 
                phase_locking_gamma = gradthe_g(bbin2);
                phase_locking2_gamma = gradthe(Negppth_g);
                
                % Negative peaks to test phase location
%                 figure;rose(deg2rad(gradthe(Negppth)))%-180))                             % To check phase locking of spikes
            if plott==1
                for lo=1:length(bbin2)
                    text(Eptime(bbin2( lo)),ytext(lo),'O' )   
                    text(Eptime(bbin2( lo)),ytext(lo)-0.2,num2str(phase_locking(lo)))
                end
            end
            phase_locking_out{count7,ci} = phase_locking ;
            phase_locking_out2{count7,1} = phase_locking2 ;
            phase_locking_gamma_out{count7,ci} = phase_locking_gamma;
            phase_locking_gamma_out2{count7,1} = phase_locking2_gamma;
            end
        end
%%%%%%% _end of function for Anal1
    end
end
edges = (-pi:(2*pi)/nBins:pi);
v = genvarname(['Cell_ ', num2str(unit)]);
v2 = genvarname(['Cell_',num2str(unit),'_gamma']);

eval([v '=cell2mat(transpose(phase_locking_out(:,unit)))'])
eval([v2 '=cell2mat(transpose(phase_locking_gamma_out(:,unit)))'])

totalUse = length(eval(v));
totalUse_gamma = length(eval(v2));
data2 = deg2rad(eval(v)-180);
data2_gamma = deg2rad(eval(v2)-180);

data3 = deg2rad(eval(v));
data3_gamma = deg2rad(eval(v2));
hist = histc(data2,edges);
hist(end) = [];
hist_gamma = histc(data2_gamma,edges);
hist_gamma(end)=[];
% Rayliegh test
[p,z]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist);
[mu ul ll] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist,2);  
meanAn=rad2deg(mu)+180; 
degCell=eval(v);
kkk=1:360/(nBins):360; % Intervals for deg hist
[oo pp] = histc(degCell,kkk);
InfInter=rad2deg(ul)+180 ;
SupInter=rad2deg(ll)+180 ;
mvl = circ_r((edges(1:end-1)+edges(2:end))/2,hist,[],2);

% Rayliegh test Gamma
[p_gamma,z_gamma]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist_gamma);
[mu_gamma ul_gamma ll_gamma] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist_gamma,2);  
meanAn_gamma=rad2deg(mu_gamma)+180; 
degCell_gamma=eval(v2);
kkk_gamma=1:360/(nBins):360; % Intervals for deg hist
[oo_gamma pp_gamma] = histc(degCell_gamma,kkk_gamma);
InfInter_gamma=rad2deg(ul_gamma)+180 ;
SupInter_gamma=rad2deg(ll_gamma)+180 ;
mvl_gamma = circ_r((edges(1:end-1)+edges(2:end))/2,hist_gamma,[],2);




%% OUTPUT VALUES
pld.totalUse = totalUse;
pld.hist = hist;
pld.p = p;
pld.z = z;
pld.mu = mu;
pld.ul = ul;
pld.ll = ll;
pld.meanAn = meanAn;
pld.InfInter = InfInter;
pld.SupInter = SupInter;
pld.mvl = mvl;
pld.phasebins = kkk;

pld_gamma.totalUse = totalUse_gamma;
pld_gamma.hist = hist_gamma;
pld_gamma.p = p_gamma;
pld_gamma.z = z_gamma;
pld_gamma.mu = mu_gamma;
pld_gamma.ul = ul_gamma;
pld_gamma.ll = ll_gamma;
pld_gamma.meanAn = meanAn_gamma;
pld_gamma.InfInter = InfInter_gamma;
pld_gamma.SupInter = SupInter_gamma;
pld_gamma.mvl = mvl_gamma;
pld_gamma.phasebins = kkk_gamma;
end

