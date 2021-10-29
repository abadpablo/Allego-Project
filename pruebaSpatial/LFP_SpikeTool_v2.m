

function  [OutDeg]=LFP_SpikeTool_v2(SK1,CUTO1,duration,mouse,group,experiment,filenamCon1,nrtet,ExperCell)

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
% 1.  Need further testing and would be based in thw work of Csicvari 1999
 
%% B. Load Cutting Info
%... of all 3 conditions
% Tetrode organizing

if str2double(nrtet)==1
    filename=[filenamCon1,'.egf'] ;
elseif str2double(nrtet)==2
    ch=5;
elseif str2double(nrtet)==3
    ch=9;
elseif str2double(nrtet)==4
    ch=13;
elseif str2double(nrtet)==5
    ch=17;
elseif str2double(nrtet)==6
    ch=21;
elseif str2double(nrtet)==7
    ch=25;
elseif str2double(nrtet)==8
    ch=29;  
end
    if str2double(nrtet)>1
   filename=[filenamCon1,'.egf',num2str(ch)] 
    end
% filename=[filenamCon1,'.egf',num2str(ch)]
[CA1, Fs] = getegf(filename);
CA1=decimate(CA1,2)'/8738.1; % This is an aprox value for amp correction obtained by calculating the ratio between Original EGF and the imported data from Nex.
Args.Fs=Fs/2;
RunOne=1;
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

kk=cell2mat(SK1.spike);
SpikeWaveform=[kk(:,:,1) kk(:,:,2) kk(:,:,3) kk(:,:,4)];
CellIdALL=cell2mat(CUTO1.neuronid);
spiketime=cell2mat(SK1.spiketime); 

nameVar=zeros(length(spiketime),203);
nameVar(:,2)=CellIdALL;
nameVar(:,3)=spiketime(:,1);
nameVar(:,4:end)=SpikeWaveform;

% [sett] = readsettingsTSpk([filenamCon1 '.set'] ,str2num(nrtet));% Gets settings from recording
% duration=str2num(sett.Duration);
% Cell1=8;
% CellId1=find(CellIdALL==Cell1);
% CellId1=CellId1';
% CellTS1=spiketime(CellId1);

%_____________________  LFP slecting and pre-processing_________________

TheL1=4;
TheL2=12;
% ch=5; % Channel to analyse  LFP Vs Spike
% CA1=Args2.LFP(ch,1:end); 
CA1=Notched(CA1,Args,50,100); % NOTCH FILTER FOR 50 an 100 Hz Arguments (Signal, Frequency 1 , Frequency 2)
CA1=locdetrend(CA1',Args.Fs,[2 0.25])';
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
 
if (RatTD >6) &&  (max(thetaEnv_ep)< Threstheta)  % Only Epochs of strong theta are included for analysis eliminating high amplitude artifacts

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
for ci=1:max(CellIdALL);
Cell1=ci;
CellId1=find(CellIdALL==Cell1);
CellId1=CellId1';
CellTS1=spiketime(CellId1,1);%;
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
end
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
 OutDeg=zeros(max(CellIdALL),12);
%  OutSpikes=zeros(max(CellIdALL),2);
% for ii=1:max(CellIdALL)
% 
% end
 for ci2=1:max(CellIdALL)

 [Sm1 Sm2 SpikeAmpD WidthD]=SpikeWaveformPlotOG (SK1,CUTO1,ci2,ExperCell);
%  OutSpikes(ii,1:end)=[SpikeAmpD WidthD]
  v = genvarname(['cell_Deg_',num2str(ci2) ]);
  
  CellId1=find(CellIdALL== ci2);

  
  
  eval([v '=cell2mat(transpose(phase_locking_out(:,ci2)));']) 
  NSpikes=length(CellId1);
  MeanF= NSpikes/duration;
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
% figure;bar(oo/NSpikes);
 InfInter=rad2deg(ul)+180 ;
 SupInter=rad2deg(ll)+180 ;
 mvl=circ_r((edges(1:end-1)+edges(2:end))/2,hist)
 
 if MeanF>0.2 && MeanF<6.9  && WidthD>300 && p<0.05
     
      Type='Pyr';
    ExperCell_Cell=[experiment,'_',nrtet,'_c',num2str(ci2),'_',num2str(group),'_',Type];
    save(['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\OutPutDegress_Pyr\',ExperCell_Cell] ,'meanAn','mu','ul','ll','p','z','degCell','mvl','meanAn','oo','MeanF','NSpikes','totalUse');
elseif MeanF>0.2 && MeanF>6.9  && WidthD>300 % Pyr
Type='PyrUnc';
    ExperCell_Cell=[experiment,'_',nrtet,'_c',num2str(ll),'_',num2str(group),'_',Type]
 save(['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\OutPutDegress_Pyr\',ExperCell_Cell] ,'meanAn','mu','ul','ll','p','z','degCell','mvl','meanAn','oo','MeanF','NSpikes','totalUse');
elseif MeanF>7 && WidthD<300 && p<0.05 % Inter
    Type='Inter'
       ExperCell_Cell=[experiment,'_',nrtet,'_c',num2str(ci2),'_',num2str(group),'_',Type];
    save(['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\OutPutDegress_Inter\',ExperCell_Cell] ,'meanAn','mu','ul','ll','p','z','degCell','mvl','meanAn','oo','MeanF','NSpikes','totalUse');   
elseif MeanF<7 && WidthD<300% Inter
 Type='InterUnc'
       ExperCell_Cell=[experiment,'_',nrtet,'_c',num2str(ci2),'_',num2str(group),'_',Type];
    save(['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\OutPutDegress_Inter\',ExperCell_Cell] ,'meanAn','mu','ul','ll','p','z','degCell','mvl','meanAn','oo','MeanF','NSpikes','totalUse');   
 
 end
     
     
 OutDeg(ci2,1:end )=[ ci2 SpikeAmpD WidthD  MeanF meanAn InfInter SupInter mvl z p  NSpikes totalUse ];
clear data2 mu ul ll p z  hist NSpikes  totalUse
 end
 
%  filename4=[num2str(Args.date),'-',num2str(Args.mouse), num2str(Args.Mutant),Args.ExP,'_',num2str(Args2.cond)];

%  HeaderDeg ={'ci2','SpikeAmpD','WidthD','MeanF','meanAn','InfInter','SupInter','mvl','z','p','NSpikes','totalUse'  }
%   xlswrite( ['Phase_deg',num2str(Args.Mutant)], HeaderDeg,[filename4,'_out'],['B1']) % Power Speed Correlation
%  xlswrite( ['Phase_deg',num2str(Args.Mutant)], OutDeg,[filename4,'_out'],['B2']) % Power Speed Correlation
 
% experiment='2411CCK-R229_OF3'

% xlswrite('Spatial_CCK1_sta',  HeaderDeg,experiment, ['BC1'])
% xlswrite('Spatial_CCK1_sta', OutDeg,experiment, ['BC2'])
% %------
% if lmk>1 && lmk+1<length(indx3)
%     count4=count4+1; 
% kh=round(abs((length(th2)-320))/2);
% 
% thMean2{count4,1}=th((Negppth(lmk))-kh:(Negppth(lmk+1)+kh)); 
% meany{count4,1}=y((Negppth(lmk))-kh:(Negppth(lmk+1)+kh)); 
% 
% end
% 
% % Mat variables vectors
% Mat1=((S2(:,indx3(lmk):indx3(lmk+1))));     % Original Spectogram Matrix  
% % Mat2=(10*log10(I(indx4(lmk):indx4(lmk+1),:)*10e12)');   % Resized MATRIX
% 
% 
% thres1=mean(Mat1(:));                                 % Power spectra threshold for CM calculation (Original)
% [x1 y1]=centro_masa2(Mat1,thres1);                       % CM calculation of the Spectral Matrix (Original)
% % in=find(Mat1<thres1);
% % Mat1(in)=0;
% % figure;imagesc(Mat1)
% % thres2=mean(mean (Mat2));hold;                          % Power spectra threshold for CM calculation (Resized)
% % [xx2 yy2]=centro_masa(Mat2,thres2)                      % CM calculation of the Spectral Matrix (Resized)
% [fil col]=size(Mat1);
% VarFil=(30:1:fil);
% [bbin5]=histc(VarFil,f2);
% [bbin4] =histc(tt22,EptInt2);                            % [When USING MULTITAPER] bining time resolution from spectogram interval, using eeg resolution times
% DDgradRed=gradthe2(find(bbin4==1));                      % Index finding of the corresponding values to each interval 
% Fre=f2(round(y1));
% 
% %%%%%%%% Continuar aqui
% % repasar el bineado para encontrar CgMC y retomar para encontrar la
% % frequencia asociada al CM
% % Fbin=
% % Degrees calculation of the CM in the power spectra matrix 
%                      % Index finding of the corresponding values to each interval 
% 
% % Multitaper=1;
% % if Multitaper==1
% 
% % [bbin4] =histc(tt22,EptInt2);                            % [When USING MULTITAPER] bining time resolution from spectogram interval, using eeg resolution times
% DDgradRed=gradthe2;%(find(bbin4==1));                      % Index finding of the corresponding values to each interval 
% timeRed=EptInt2(find(bbin4==1));                       
% GcMG1=DDgradRed(round(x1)+1);                              % Ocurrence of CM in degrees
% % else
% timeRed=tt22;
% GcMG1=DDgradRed(round(x1)); % using the degress for the short epoch calculated ad hoc generates a shwift towards the right
% GcMG2=Dgradthe2(round(x1)); % using the corresponding epoch  obtained from the peaks intervals seems more correct
% % end
% % Plotting values
% figure;
% subplot(2,2,1);hold;plot(EptInt2,y2);   ylabel('Amp');                    % RAW data 
% subplot(2,2,1);plot(EptInt2,th2,'r') ;plot(EptInt2,Dgradthe2/1000,'c');     % Theta and corresponding degrees 
% axis([EptInt2(lmk)  EptInt2(end) -0.8 0.8]);              % RAW data 
% text(EptInt2(ppthS),mmthS,num2str(Dgradthe2(ppthS))); 
% text(timeRed(round(x1)),Dgradthe2(round(x1))/1000,['G ',num2str(GcMG2)]);  % Plotting CM degrees !!! 
% subplot(2,2,3);imagesc(Mat1);hold;text(x1,y1,num2str(GcMG2));ylabel('Power (Au');
% 
% %%%%%%%%%% data for calculated on short
% subplot(2,2,2);hold;plot(EptInt2,y2);   ylabel('Amp')   
% subplot(2,2,2);plot(EptInt2,th2,'r') ;plot(EptInt2,DDgradRed/1000,'g')     % Theta and corresponding degrees 
% text(EptInt2(ppthS),mmthS,num2str(Dgradthe2(ppthS))); 
% text(timeRed(round(x1)),DDgradRed(round(x1))/1000,['G ',num2str(GcMG1)]);  % Plotting CM degrees !!! 
% axis([EptInt2(lmk)  EptInt2(end) -0.8 0.8]); 
% % subplot(2,1,2);figure; imagesc(tt22,f2,10*log10(Mat1*10e12));
% figure;imagesc(Mat1);text(x1,y1,num2str(GcMG1));ylabel('Power (Au');
% %   set(gca,'XTickLabel',num2str(DDgradRed ))
%   %
%   figure; imagesc(tt22,f2,10*log10(Mat1*10e12));
% %   figure;imagesc(Mat1);text(x1,y1,'XX');ylabel('Power (Au');
% %   set(gca,'XTickLabel',num2str(DDgradRed ))
% %  pause;
% % pause %(0.2)
% close  all
% Mean Power for the different gamma subbands based on the wavelet Matrix
% Matlow=Mat1(1:31,:);%fr=30-60
% MatInter=Mat1(32:61,:);%fr=60-90
% MatHigh=Mat1(62:121,:);%fr=91-125

% if EmptyPeakTheta==0
%    OutPut1{count1}=GcMG1;% Center of mass based on degress calculated  for short epochs
%    OutPut2{count1}=GcMG2;% Center of mass based on degress calculated  for the 5 secs epoch
%    OutPut3{count1}=Fre;       % Frequency of gamma in the center of mass
% 
% 
% end
% 
    end
% end
%
%      end
% end
% end
%  end
 % CM: Degrees vs Frequency
% scatter(cell2mat(OutPut1),cell2mat(OutPut3));axis([0  360 30 150]);  hold; text(mean(cell2mat(OutPut2)),mean(cell2mat(OutPut3)),'O'), title([num2str(count1), ' Cicles / ',num2str(count2) ' s'])
% edges=(-pi:0.3307:pi);
% data2=deg2rad(cell2mat(OutPut2)-180);
% % R109=[cell2mat(OutPut1);cell2mat(OutPut2);cell2mat(OutPut3)];
% hist=histc(data2,edges); hist(end)=[];
% % Rayleight test
% [p,z]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist);  %% 
% % mean circular 
% [mu ul ll] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist);  meanAn=rad2deg(mu)+180  ; text(meanAn,mean(cell2mat(OutPut3)),'X')
% mvl=circ_r((edges(1:end-1)+edges(2:end))/2,hist);% mean vector length
% 
% mu=circ_rad2ang(mu)+180;
% ul=circ_rad2ang(ul)+180;
% ll=circ_rad2ang(ll)+180;
% figure; hold
% waveM=zeros(1,320);
%  for i=1:length(thMean2);
%     wav=cell2mat(thMean2(i));
%     
%     if   length(wav)>=321
%         wav=wav(1,1:320);
%        
% 
%     end
%     plot(wav)
%     waveM(i,1:end)=wav;% waveform cleaned from longer waveforms
%      
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%% MeanRAW
% RawwaveM=zeros(1,320);
% 
% figure;hold
% theR_mean=zeros(124,320,length(meany));
% for   ii=1:length(meany);
%     Rawmeany=cell2mat( meany(ii));
%     if   length(Rawmeany)>=321
%        Rawmeany=Rawmeany(1,1:320);
%     end
% %     subplot(2,1,1);plot(Rawmeany)  
%     [theR,theR_t,theR_F] = nbt_TF_plot2( Rawmeany,0.1,30,150,1,0,Args);%[S2,t,f2] = nbt_TF_plot2(y,0.1,30,150,1,0,Args);
%     RawwaveM(ii,1:end)=Rawmeany;
%     theR=sqrt(abs(theR));theR=flipud(theR);
% %     subplot(2,1,2);imagesc(theR_t,theR_F,theR)%imagesc(tt22,f2,Mat1)
%     theR_mean(:,:,ii) =theR;
% %     pause
%    
% end
% figure;imagesc(theR_t,theR_F, mean(theR_mean,3))
% Args3.MeanWavelete=mean(theR_mean,3);
% close all
% % figure;imagesc(10*log10(Args2.MeanWavelete*10e12))
%  end
% Figure



    

% Args3.thetagamma=[cell2mat(OutPut1);cell2mat(OutPut2);cell2mat(OutPut3)]
% Args4.ThetaMeanWave=thMean2;
% Args5.CleanThetaWave=waveM;
% % Args2.MeanWavelete=mean(theR_mean,3);
% save(['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\LAB\R',num2str(Args.mouse),'_LFP\',num2str(Args.date),'-',num2str(Args.Mutant),'-',num2str(Args2.mouse),'_',Args2.ExP,num2str(Args2.cond),'_Pos'], 'Args3','Args2');
% ['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\LAB\R',num2str(Args.mouse),'_LFP\',num2str(Args.date),'-',num2str(Args.Mutant),'-',num2str(Args2.mouse),'_',Args2.ExP,num2str(Args2.cond),'_Pos2']
% ['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\LAB\R',num2str(Args.mouse),'_LFP\',num2str(Args.date),'-',num2str(Args.Mutant),'-',num2str(Args.mouse),'_',Args.ExP,num2str(Args.ExP),'_Header']





% end



% pause
% close
% end

%  [Sc,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMatc(CA1',0.5,params)
%  
%  figure;imagesc(Sc');hold
% tstep