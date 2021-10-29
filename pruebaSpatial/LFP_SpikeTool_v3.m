

function  [totalUse hist p z mu ul ll meanAn InfInter SupInter mvl]=LFP_SpikeTool_v3(CA1,times_unit,data_aux,nrtet,iClust,reg_rec)

% (C) Jorge R. Brotons-Mas (2015)
% This script is and attemp to obtain the relationship between spike
% ocurrence and theta phase. 
% It is based in params3 and different work such as that from Hangya "
% Complementary spatial firing in place cell-interneuron pairs
% Algorithim 1:
% 1. Signal is filtered between 4-12 hz and the hilbert transform obtained to
% determine phase of the ongoin cycle. 
% 2. The signal is broken i 2 second epochs and the theta/ delta ratio
% calculated, if this >6 then the epoch is entered for analysis.
% 3 Timestamps of spikes is obtained by binning this with the LFP clock.
% 4. Data is organized and circular statistics calculated.
% 5. Finally data is save in different files

% Input:
% SK1,CUTO1,duration,mouse,experiment,filenamCon1,nrtet
% OutPut:
% mu ul ll p z  hist NSpikes  totalUse

% Algorithim 2:
% % 1.  Need further testing and would be based in thw work of Csicvari 1999
%  
% %% B. Load Cutting Info
% %... of all 3 conditions
% % Tetrode organizing
% tintimportedold=0
% if tintimportedold==1
% 
% if str2double(nrtet)==1
%     filename=[filenamCon1,'.egf'] ;
% elseif str2double(nrtet)==2
%     ch=5;
% elseif str2double(nrtet)==3
%     ch=9;
% elseif str2double(nrtet)==4
%     ch=13;
% elseif str2double(nrtet)==5
%     ch=17;
% elseif str2double(nrtet)==6
%     ch=21;
% elseif str2double(nrtet)==7
%     ch=25;
% elseif str2double(nrtet)==8
%     ch=29;  
% end
%     if str2double(nrtet)>1
%    filename=[filenamCon1,'.egf',num2str(ch)] 
%     end
% % filename=[filenamCon1,'.egf',num2str(ch)]
% [CA1, Fs] = getegf(filename);
% CA1=decimate(CA1,2)'/8738.1; % This is an aprox value for amp correction obtained by calculating the ratio between Original EGF and the imported data from Nex.
% Args.Fs=Fs/2;
% RunOne=1;
% end
% 
% %%%%%% NEW IMPORTED
% 
% if nrtet==1
%     ch=9
% elseif nrtet==2
%     ch=5
% elseif nrtet==3
%     ch=9
% elseif nrtet==4
%     ch==13
% end
% 
% LFP_analysis=data_aux(ch,1:end);
Args.Fs=2400;
% LFP_analysis=decimate(LFP_analysis,2);
% 
% [LFP_analysis]=Notched(LFP_analysis,Args,50,100); % Notched for 50 and 100 Hz
% LFP_analysis = locdetrend(LFP_analysis',Args.Fs,[1 0.2])';
%   CA1=LFP_analysis;
RunOne=1
if RunOne==1
% filenamCon1='D:\CCK\CCK229\2411CCK-R229\2411CCK-R229_OF3';%[dire,filename]
% ExperCell='2411CCK-R229_OF3';% =[experiment,'_',num2str(nrtet),'_','c',num2str(nrcell)];
% nrtet='3'
% GridTet=nrtet;
% SK1=spikeobj;   %create spikeobject for three conditions
% SK1=importToolSpk(SK1,'AXONA',[filenamCon1 '.' nrtet]); %load the spikes of axona of all 3 Cond.
%  %%%%
% % B. Load Cutting Info
% CUTO1=cutobj; %***
% [CUTO1]=import(CUTO1,'TINT',[filenamCon1 '_' nrtet 'CCD.cut']);

%_____________________ Organizing Spike data _____________________ 

% kk=cell2mat(SK1.spike);
% SpikeWaveform=[kk(:,:,1) kk(:,:,2) kk(:,:,3) kk(:,:,4)];
% CellIdALL=cell2mat(CUTO1.neuronid);
% spiketime=cell2mat(SK1.spiketime); 
% 
% nameVar=zeros(length(spiketime),203);
% nameVar(:,2)=CellIdALL;
% nameVar(:,3)=spiketime(:,1);
% nameVar(:,4:end)=SpikeWaveform;
% 
% % [sett] = readsettingsTSpk([filenamCon1 '.set'] ,str2num(nrtet));% Gets settings from recording
% % duration=str2num(sett.Duration);
% % Cell1=8;
% % CellId1=find(CellIdALL==Cell1);
% % CellId1=CellId1';
% % CellTS1=spiketime(CellId1);

%_____________________  LFP slecting and pre-processing_________________
addpath('C:\Users\Jorge\OneDrive - Fundación Universitaria San Pablo CEU\MATLAB\GlobAnalysis-1.0\LAB\SignalProcessing\PowerMak')
TheL1=4;
TheL2=12;
% ch=5; % Channel to analyse  LFP Vs Spike
% CA1=Args2.LFP(ch,1:end); 
% CA1=Notched(CA1,Args,50,100); % NOTCH FILTER FOR 50 an 100 Hz Arguments (Signal, Frequency 1 , Frequency 2)
% CA1=locdetrend(CA1',Args.Fs,[2 0.25])';
BCC = fir1(400,[TheL1/(Args.Fs/2) TheL2/(Args.Fs/2)],'bandpass');         % theta filter convolution bands as in Monyer paper Neuron
th_all=filtfilt(BCC,[1],CA1) ;%theta band filtered
thetaEnv = abs(hilbert(th_all));  % envelop
Threstheta2=mean(thetaEnv)-std( thetaEnv ) ;   % threshold for peak amp detection
  
end
% figure; pwelch(CA1,[],[],[],2400,'twosided'); % helps to understand the
% filter

%______________________Time vector construction  __________________________

k=length(CA1);
Fs=Args.Fs;
EEG_ts_step=1/Fs;
tstep=EEG_ts_step*(1:k);

%___________________ Counts for organizing data and loop __________________

 count1=0;% for phase calculation
 count2=0; % for epochs of 2 secs meeting theta/Delta ratio
 count3=0; % cicles meeting criteria for analysis
 count4=0; % counter for theta avergage cicle
 count7=0; % count to use in Anal1
 
 %____________________ Start Loop for data selecting _____________________
 for kl=2:2:(length(tstep))/Args.Fs-2
    t1=kl ;                                             % Time 1 of epoch 
    t2=kl+2;                                            % Time 2 of epoch
    y=CA1(1,Fs*t1:Fs*t2);                               % Signal epoch 
    Eptime=tstep(1,Fs*t1:Fs*t2);                        % Epoch timed vector
    EptInt=EEG_ts_step*(1:length(Eptime));              % helps to find interval for degress CM position
    y=locdetrend(y',Args.Fs,[0.5 0.2])';                                     % Data detrended                                      
    [Eband Eband2 RatTD p  f] =MaKPower2(y,Args,0);            % Power spectra calculation for theta/delta ratio
%% Theta Filter
BC = fir1(400,[TheL1/(Fs/2) TheL2/(Fs/2)],'bandpass');         % theta filter convolution bands 
th=filtfilt(BC,[1],y) ;%theta band filtered

% BC = fir1(400,[4/(Args.Fs/2) 12/(Args.Fs/2)],'bandpass');% filter generation
%     thetaF=filtfilt(BC,[1],CA1); % filter signal
thetaEnv_ep = abs(hilbert(th));  % envelop
Threstheta=mean(thetaEnv)+std( thetaEnv )*5 ;   % threshold  for noise detection    
% Threstheta2=mean(thetaEnv)+std( thetaEnv ) ;  % threshold for peak amp detection
%% Gamma
B = fir1(400,[25/(Fs/2) 150/(Fs/2)],'bandpass');        %gamma filter convolution as Buzsaki 2003
gam= filtfilt(B,[1],y); %gamma filtered

if reg_rec==1
tresh_1=3;
elseif reg_rec==2
    tresh_1=1.2;
end


if (RatTD >tresh_1) &&  (max(thetaEnv_ep)< Threstheta)  % Only Epochs of strong theta are included for analysis eliminating high amplitude artifacts

    count2=count2+1;



%% theta trhough detection
[Negmmth Negppth] = findpeaks((-th),'minpeakdistance',round(Fs*0.07));%,'minpeakheight',Threstheta2);  % Using Amp thresholdto detect poor theta cycles    % Theta Valley 
%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 

% [mmth ppth] = findpeaks(th,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.09));                                                % Theta Peaks 
%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 
ThAmp=mean(abs(th)); Stdth=std((abs(th)));                                          % Mean values of amplituded 
AmpThr=ThAmp-Stdth;%-Stdth;                                                         % Theta Amplitude Signal thresholded
% [mmGa ppGa] = findpeaks(gam);                                               % Gamma Peaks detection
%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 


%% Theta phase detection and degress epoch detection

 X = hilbert(th) ;                                                          % Hilbert transformed
 phi = angle(X);                                                             % angle of the different signal samples
%     y = abs(phi);
gradthe=rad2deg(phi)+180;                                                     % corrected to 0-360º
[mmGradth ppGradth] = findpeaks(gradthe,'minpeakdistance',round(Fs*0.07));     % finding neg peaks so transitons between cicles

%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 
%% PLoting data
plott=0;
Anal1=1;
if Anal1==1
    count7=count7+1;
    if plott==1;
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

% subplot(2,1,2); imagesc(tt2,f2,S3); colorbar('SouthOutside' );axis xy % Plots waveltes

subplot(3,1,3);plot( Eptime,y);hold;plot( Eptime,th,'r');plot( Eptime,gradthe/1000,'c' );   %plot Degrees Signal 
% pause% 
    end
   
% for ci=1:max(CellIdALL);
ci=iClust;
% Cell1=ci;
% CellId1=find(CellIdALL==Cell1);
% CellId1=CellId1';

CellTS1=times_unit;%;
VectimeSpike=find(CellTS1>=Eptime(1) & CellTS1<=Eptime(end));
VectimeSpike_2=CellTS1(VectimeSpike);
% CellTS1(CellTS1>=Eptime(1) & CellTS1<=Eptime(end));
[tbin bbin2] =histc(VectimeSpike_2,Eptime); % trying to bin spike times to Eptime 

intT=find(bbin2>0);
% Vec1=Eptime(intT);
ytext=zeros(1,length(intT));
phase_locking=gradthe(bbin2);                                               % Spike Phase Locking                                            
phase_locking2=gradthe(Negppth);                                            % Negative peaks to test phase location
% figure;rose(deg2rad(gradthe(Negppth)))%-180))                             % To check phase locking of spikes
 if plott==1
for lo=1:length(bbin2)
 text(Eptime(bbin2( lo)),ytext(lo),'O' )   
text(Eptime(bbin2( lo)),ytext(lo)-0.2,num2str(phase_locking(lo)))

end
 end
% % figure;  rose(deg2rad(phase_locking))%-180))
% % 
% pause (0.5) 
close all

phase_locking_out{count7,ci}=phase_locking ;
phase_locking_out2{count7,1}=phase_locking2 ;
% end
end
%%%%%%% _end of function for Anal1
 end
 end
% TestingOutputLFP_SpikeTool
%% Cicles analysis
% %  for lmk=4 
% CycleAnal=0;
% % [bbin3] =histc(Eptime(Negppth),t+kl);                     % bining epoch times to negative degress phase
% %                                                         % Also I have proved against  the EEG filtered negative peaks  
% % indx3=find(bbin3==1);        
% if CycleAnal==1;
%     for  lmk= 1:length(Negppth)-1
% count1=count1+1;
% count3=count3+1
%         % Variable generation for each cicle
%         % Signal vectors
% 
% Dgradthe2=(gradthe (Negppth(lmk):(Negppth(lmk+1))));
% y2=(y(Negppth(lmk):(Negppth(lmk+1))));                % Raw epoched signal
% th2=(th(Negppth(lmk):(Negppth(lmk+1))));              % theta filtered epoched signal
% Pole=round(200*length(y2)/Args.Fs);
% % testing for refiltering after cutting 
% pole=round(length(y2)/5);
% BC = fir1(pole,[TheL1/(Args.Fs/2) TheL2/(Args.Fs/2)],'bandpass');         % theta filter convolution bands as in Monyer paper Neuron
% th2_2=filtfilt(BC,[1],y2) ;%theta band filtered
% 
% 
% %% Theta phase detection and degress epoch detection  for short epochs
% 
%  X2 = hilbert(th2) ;                                                   % Hilbert transformed
%  phi2 = angle(X2);                                                     % angle of the different signal samples
% %     y = abs(phi);
% gradthe2=rad2deg(phi2)+180;                                                     % corrected to 0-360º
% %------Inserted for plotting and checking the function _____
% if plott==1
% figure;hold;plot(y2)
% % 
% plot(th2,'r')
% plot(Dgradthe2/1000,'c')
% plot(gradthe2/1000,'r')
% 
% [mmGradth2 ppGradth2] = findpeaks(-gradthe2,'minpeakdistance',round(Fs*0.06));     % finding neg peaks so transitons between cicles
% 
% %,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascen
% %d'); % Peaks Above 20 Hz 
% close all
% end
% % [NegmmthSS ppthS] = findpeaks((-th2)),'minpeakdistance',round(Fs*0.05))%,'minpeakheight',Threstheta2);  % Using Amp thresholdto detect poor theta cycles    % Theta Valley 
% 
% % EmptyPeakTheta=isempty(ppthS);                                       % Cover for error of non peak detected
% % if EmptyPeakTheta==0
%     if (length(th2)>round(Fs*0.03)) && (length(th2)<=round(Fs*0.11))     % Amplitude and too long or too short periods                     
% %   for kjl=1:length(indx3)-1
% % vallim=indx3(end);
% %-----
% for ci=1:max(CellIdALL)
%     count7=count7+1;
% Cell1=ci;
% CellId1=find(CellIdALL==Cell1);
% CellId1=CellId1';
% CellTS1=spiketime(CellId1,1);%;
% VectimeSpike=find(CellTS1>=Eptime(Negppth(lmk)) & CellTS1<=Eptime(Negppth(lmk+1)));
% VectimeSpike_2=CellTS1(VectimeSpike);
% % CellTS1(CellTS1>=Eptime(1) & CellTS1<=Eptime(end));
% [tbin bbin2] =histc(VectimeSpike_2,Eptime); % trying to bin spike times to Eptime 
% 
% intT=find(bbin2>0);
% % Vec1=Eptime(intT);
% ytext=zeros(1,length(intT));
% phase_locking=gradthe(bbin2);
% 
% % figure;rose(deg2rad(gradthe(Negppth)))%-180))  % to check phase locking of spikes
% % figure;  rose(deg2rad(phase_locking))%-180))
% %
% if plott==1
% figure;hold;plot(Eptime(Negppth(lmk):Negppth(lmk+1)),y2)  
% % 
% plot(Eptime(Negppth(lmk):Negppth(lmk+1)), th2,'r') 
% plot(Eptime(Negppth(lmk):Negppth(lmk+1)),Dgradthe2/1000,'c') 
% plot(Eptime(Negppth(lmk):Negppth(lmk+1)),gradthe2/1000,'r') 
% 
% text(Eptime(bbin2),ytext,'O')
% 
% end
% % pause(0.5)
% close all
% 
% phase_locking_out{count7,ci}=phase_locking ;
% end
%% OUTPUT

edges=(-pi:0.3307:pi);  % preparing intervals for phase locking testing
 OutDeg=zeros(1,12);
%  OutSpikes=zeros(max(CellIdALL),2);
% for ii=1:max(CellIdALL)
% 
% end
%  for ci2=1:max(CellIdALL)

%  [Sm1 Sm2 SpikeAmpD WidthD]=SpikeWaveformPlotOG (SK1,CUTO1,ci2,ExperCell);
%  OutSpikes(ii,1:end)=[SpikeAmpD WidthD]
  v = genvarname(['cell_Deg_',num2str(iClust) ]);
  
%   CellId1=find(CellIdALL== ci2);
addpath('C:\Users\Jorge\OneDrive - Fundación Universitaria San Pablo CEU\MATLAB\GlobAnalysis-1.0\LAB\CircStat')
  
  
  eval([v '=cell2mat(transpose(phase_locking_out(:,iClust)));']) 
%   NSpikes=length(time_unit);
%   MeanF= NSpikes/duration;
  totalUse=length(eval(v))
 data2= deg2rad(eval(v) -180);
%  data2= circ_ang2rad(eval(v) -180);%%%%%%%%%%%%%%%%% CHECK FUNCTION  
 data3=deg2rad(eval(v));
% data2= deg2rad(gradthe(Negppth)-180)
 hist=histc(data2,edges); hist(end)=[];
% Rayleight test
 [p,z]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist)   %% 

 [mu ul ll] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist);  
 meanAn=rad2deg(mu)+180  
 
 
%  meanAn=circ_ang2rad(mu)+180 %%%%%%%%%%%%%%%%% CHECK FUNCTION 
%  figure;rose(data3) ; title(['Cell ',num2str(ci2),'Mean º',num2str(meanAn)])
%  pause
degCell=eval(v);
kkk=1:20:360; % Intervals for deg hist
[oo pp]=histc(degCell,kkk);
% figure;bar(oo/length(times_unit));
 InfInter=rad2deg(ul)+180 ;
 SupInter=rad2deg(ll)+180 ;
 mvl=circ_r((edges(1:end-1)+edges(2:end))/2,hist)
