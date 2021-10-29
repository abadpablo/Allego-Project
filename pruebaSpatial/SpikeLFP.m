
% load(['C:\Users\Jorge\Documents\MATLAB\gSpike5\LAB\LFPAnalysis192\1511ERB4-192_OF1_Pos.mat'])
% CA1=LFP(1,4800:4800+2400);
% params= struct('Fs',2400,'pad',0,'fpass',[30 90],'err',0,'tapers',[0.1
% 100 (2*(0.1* 100)-1)]) 
% T=0.06%input('time interval ?');
% W=input('Frequency  ?');
% P=(2*(T*W)-1)
%% Filtering and cleaning
disp( ' Warning, check recording duration in script ')
pause
CA1=Args2.LFP(ch,1:Args.Fs*295); 

%  [Hd b a] = NotchIrFir100 ; % In case of need 50Hz deletion
%  [Hd b a] = NotchIrFir50 ; % In case of need 50Hz deletion
 CA1=Notched(CA1,Args,50,100); % NOTCH FILTER FOR 50 an 100 Hz Arguments (Signal, Frequency 1 , Frequency 2)
CA1=locdetrend(CA1',Args.Fs,[2 0.25])';

% CA1Notch=filter(b,a,CA1);
% [CA1]=ArtifactDetection(CAR,2400,0); % Atifact elimination
% figure; pwelch(CA1,[],[],[],2400,'twosided'); % helps to understand the
%%% effect of the filter

% Time vector construction

k=length(CA1);
Fs=Args.Fs;
EEG_ts_step=1/Fs;
tstep=EEG_ts_step*(1:k);
% 
% B = fir1(400,[25/(Fs/2) 125/(Fs/2)],'bandpass'); %gamma filter convolution as Buzsaki 2003
% gamNoNotch= filtfilt(B,[1],CA1); %gamma filtered

% Ttime_ts_step=1/100;
% timeT=Ttime_ts_step*(1:k1);
 count1=0;% for phase calculation
 count2=0; % for epochs of 2 secs meeting theta/Delta ratio
 count3=0; % cicles meeting criteria for analysis
 count4=0; % counter for theta avergage cicle
% figure 
% for kl=14
% thmean=zeros(1,2401);
 for kl=2:2:(length(tstep))/Args.Fs-2
    t1=kl ;                                             % Time 1 of epoch 
    t2=kl+1;                                            % Time 2 of epoch
    y=CA1(1,Fs*t1:Fs*t2);                               % Signal epoch 
    Eptime=tstep(1,Fs*t1:Fs*t2);                        % Epoch timed vector
    EptInt=EEG_ts_step*(1:length(Eptime));              % helps to find interval for degress CM position
%     y=detrend(y);                                       % Data detrended                                      
    [Eband Eband2 RatTD p  f] =MaKPower2(y,Args,0);            % Power spectra calculation for theta/delta ratio
%% Theta Filter
BC = fir1(200,[5/(Fs/2) 18/(Fs/2)],'bandpass');         % theta filter convolution bands as in Monyer paper Neuron
th=filtfilt(BC,[1],y) ;%theta band filtered

%% Gamma
B = fir1(400,[25/(Fs/2) 150/(Fs/2)],'bandpass');        %gamma filter convolution as Buzsaki 2003
gam= filtfilt(B,[1],y); %gamma filtered
%% Processing of data 
 % counter for saving CM location 
 
if RatTD >3 % Only Epochs of strong theta are included for analysis

    count2=count2+1;
%     thmean(count2,1:end)=th;
% params= struct('Fs',2400,'pad',0,'fpass',[4 25],'err',0);    % Params for gamma detection %'tapers',[T W P]) 
% [S2,tt2,f2]=mtspecgramc(y',[0.5 0.01],params);               % good results with default tapers[0.5 0.01] % gamma detection
[S2,t,f2] = nbt_TF_plot2(y,0.1,30,150,1,0,Args);
% params= struct('Fs',2400,'pad',0,'fpass',[30 90],'err',0,'tapers',[T W P]) 
S2=sqrt(abs(S2));S2=flipud(S2);
% Mat2=(10*log10(Mat*10e12));
% params= struct('Fs',2400,'pad',0,'fpass',[25 125],'err',0,'tapers',[2 3]);  % Params for theta detection 
% [S,t,f]=mtspecgramc(y',[0.05 0.005],params);                                % spectogram for theta
tt2=t;

% I = interp2(S,5);                                                           % Power matrix interpolated by a factor of 5 to visualize better

%% theta trhough detection
[Negmmth Negppth] = findpeaks((-th),'minpeakdistance',round(Fs*0.09));                                         % Theta Valley 
%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 

% [mmth ppth] = findpeaks(th,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.09));                                                % Theta Peaks 
%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 
ThAmp=mean(abs(th)); Stdth=std((abs(th)));                                          % Mean values of amplituded 
AmpThr=ThAmp;%-Stdth;                                                         % Theta Amplitude Signal thresholded
[mmGa ppGa] = findpeaks(gam);                                               % Gamma Peaks detection
%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 


%% Theta phase detection and degress epoch detection

 X = hilbert(th) ;                                                   % Hilbert transformed
 phi = angle(X);                                                     % angle of the different signal samples
%     y = abs(phi);
gradthe=rad2deg(phi)+180;                                                     % corrected to 0-360º
[mmGradth ppGradth] = findpeaks(gradthe,'minpeakdistance',round(Fs*0.07));     % finding neg peaks so transitons between cicles


%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 
%% PLoting data
plott=0;
if plott==1
figure;                                                              % plot signal and components
subplot(2,1,1);hold; plot( Eptime,y,'k','linewidth',1.50 );ylabel('Amp (micv)') % plot Raw Signal
axis([Eptime(1)  Eptime(end) -0.8 0.8]);

%axis([tstep(Fs*t1) tstep(Fs*t2) (min(CA1)-std(abs(CA1))) (max(CA1)+std(abs(CA1)))]); 
plot( Eptime,gradthe/1000,'c','linewidth',1.50 );                    % plot Degrees Signal 
text( Eptime(ppGradth),mmGradth/1000,'PeakGr');                      % texting Neg Degrees Peaks  
plot( Eptime,th,'linewidth',1.25 );                                  % Plot theta filtered signal
% text( Eptime(ppth),mmth,'T');                                        % Texting theta Peaks 
text( Eptime(Negppth),-Negmmth,'NegT');                              % Texting theta valleys
plot( Eptime,gam,'r','linewidth',1.25 );                             % Plot gamma signal  
%text( Eptime(ppGa),mmGa,'G'); % hold;plot(Eptime,gamNoNotch(1,Fs*t1:Fs*t2));

% subplot(4,1,2);imagesc(tt2+kl,f2,10*log10(S2*10e12)');
% figure;imagesc(tt2,f2,flipud(10*log10(S2*10e12))); colorbar('SouthOutside' );axis xy % Plots waveltes
subplot(2,1,2); imagesc(t,f2,S2); colorbar('SouthOutside' );axis xy % Plots waveltes
% subplot(2,1,2); imagesc(tt2,f2,S3); colorbar('SouthOutside' );axis xy % Plots waveltes

ylabel('G Power (Au)')      % Power spectrum in Db plotted
% [Szi SziX]=size(S);
% subplot(4,1,3);imagesc(t+kl,f,10*log10(S*10e12)');hold; 
% ylabel('G Power (Au)')%text( Eptime(ppth),round(SziX/2),'T');
% subplot(4,1,4);
% imagesc(t+kl,f,10*log10(I*10e12)');
% hold;ylabel('G Power (Db)');;xlabel('time (sec)')
% %
% [bbin2] =histc(Eptime(ppth),t+kl)
% subplot(4,1,4);
% text(t( bbin2==1), t( bbin2==1),'xx');
end
%% Construsticting time variable for upsampled matrix
% [SzI SzIX]=size(I);
% Fs2=1/(max(t)/SzI);        % theoretical frequency sampling rate.
% t_I_step=1/Fs2;
% timeI=t_I_step*(0:SzI)+kl; %Time vector upsampled for upsample
%% Binning spectogram time bins vs eeg timestamps for upsample
% [bbin4] =histc(Eptime(Negppth),timeI)                   %bining epoch times  of the EEG to time intervals of t (Upsample Matrix))
% % the result but for the UPsampled I matrix
% indx4=find(bbin4==1);                                   % finding positions of the time vectors for the negative peaks (Upsample Matrix))   
%% Binning times spectograms
[bbin3] =histc(Eptime(Negppth),t+kl);                   %bining epoch times to negative degress phase
%                                                         % Also I have proved against  the EEG filtered negative peaks  
indx3=find(bbin3==1);                                   % finding positions of the time vectors for the negative peaks 


%% Cicles analysis
%  for lmk=4  

    for  lmk= 1:length(indx3)-1 
count1=count1+1;
count3=count3+1
        % Variable generation for each cicle
        % Signal vectors
y2=(y(Negppth(lmk):(Negppth(lmk+1))));                % Raw epoched signal
th2=(th(Negppth(lmk):(Negppth(lmk+1))));              % theta filtered epoched signal
Pole=round(200*length(y2)/2400);
% testing for refiltering after cutting 
pole=round(length(y2)/5);
BC = fir1(pole,[5/(Fs/2) 18/(Fs/2)],'bandpass');         % theta filter convolution bands as in Monyer paper Neuron
th2_2=filtfilt(BC,[1],y2) ;%theta band filtered



gam2=(gam(Negppth(lmk):(Negppth(lmk+1))));            % Gamma filtered epoched signal
% Eptime2=Eptime(Negppth(lmk):(Negppth(lmk+1)));        % Epoched time 
Dgradthe2=gradthe(Negppth(lmk):(Negppth(lmk+1)));      % Corresponding degress
EptInt2=EptInt(Negppth(lmk):(Negppth(lmk+1)));        % Times for the cicle (0-1s no corrected by epoch number)
Intt22=find(t >=(min(EptInt2)) & t <=(max(EptInt2))); % Times intervals of the spectogram,
tt22=t(Intt22);
    % Peak detection for short signals
[mmthS ppthS] = findpeaks(th2,'minpeakheight',mean(abs(th2)));%,'minpeakdistance',round(length(th2)/2)-20);                      % Theta Peaks
%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 

[mmGaS ppGaS] = findpeaks(gam2);                      % Gamma Peaks                   
%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascend'); % Peaks Above 20 Hz 
% Faulty epochs are left out !!!!
%% Theta phase detection and degress epoch detection  for short epochs

 X2 = hilbert(th2) ;                                                   % Hilbert transformed
 phi2 = angle(X2);                                                     % angle of the different signal samples
%     y = abs(phi);
gradthe2=rad2deg(phi2)+180;                                                     % corrected to 0-360º

[mmGradth2 ppGradth2] = findpeaks(-gradthe2,'minpeakdistance',round(Fs*0.06));     % finding neg peaks so transitons between cicles

%,'minpeakheight',0.0025,'minpeakdistance',round(Fs*0.095),'sortstr','ascen
%d'); % Peaks Above 20 Hz 
EmptyPeakTheta=isempty(ppthS);                                       % Cover for error of non peak detected
if EmptyPeakTheta==0
    if (length(th2)>80 && length(th2)<350) && max(mmthS)>AmpThr/2 ;     % Amplitude and too long or too short periods                     
%   for kjl=1:length(indx3)-1
vallim=indx3(end);

if lmk>1 && lmk+1<length(indx3)
    count4=count4+1; 
kh=round(abs((length(th2)-320))/2);

thMean2{count4,1}=th((Negppth(lmk))-kh:(Negppth(lmk+1)+kh)); 
meany{count4,1}=y((Negppth(lmk))-kh:(Negppth(lmk+1)+kh)); 

end

% Mat variables vectors
Mat1=((S2(:,indx3(lmk):indx3(lmk+1))));     % Original Spectogram Matrix  
% Mat2=(10*log10(I(indx4(lmk):indx4(lmk+1),:)*10e12)');   % Resized MATRIX


thres1=mean(Mat1(:));                                 % Power spectra threshold for CM calculation (Original)
[x1 y1]=centro_masa2(Mat1,thres1);                       % CM calculation of the Spectral Matrix (Original)
% in=find(Mat1<thres1);
% Mat1(in)=0;
% figure;imagesc(Mat1)
% thres2=mean(mean (Mat2));hold;                          % Power spectra threshold for CM calculation (Resized)
% [xx2 yy2]=centro_masa(Mat2,thres2)                      % CM calculation of the Spectral Matrix (Resized)
[fil col]=size(Mat1);
VarFil=(30:1:fil);
[bbin5]=histc(VarFil,f2);
[bbin4] =histc(tt22,EptInt2);                            % [When USING MULTITAPER] bining time resolution from spectogram interval, using eeg resolution times
DDgradRed=gradthe2(find(bbin4==1));                      % Index finding of the corresponding values to each interval 
Fre=f2(round(y1));

%%%%%%%% Continuar aqui
% repasar el bineado para encontrar CgMC y retomar para encontrar la
% frequencia asociada al CM
% Fbin=
% Degrees calculation of the CM in the power spectra matrix 
                     % Index finding of the corresponding values to each interval 

% Multitaper=1;
% if Multitaper==1

% [bbin4] =histc(tt22,EptInt2);                            % [When USING MULTITAPER] bining time resolution from spectogram interval, using eeg resolution times
DDgradRed=gradthe2;%(find(bbin4==1));                      % Index finding of the corresponding values to each interval 
timeRed=EptInt2(find(bbin4==1));                       
GcMG1=DDgradRed(round(x1)+1);                              % Ocurrence of CM in degrees
% else
timeRed=tt22;
GcMG1=DDgradRed(round(x1)); % using the degress for the short epoch calculated ad hoc generates a shwift towards the right
GcMG2=Dgradthe2(round(x1)); % using the corresponding epoch  obtained from the peaks intervals seems more correct
% end
% Plotting values
figure;
subplot(2,2,1);hold;plot(EptInt2,y2);   ylabel('Amp');                    % RAW data 
subplot(2,2,1);plot(EptInt2,th2,'r') ;plot(EptInt2,Dgradthe2/1000,'c');     % Theta and corresponding degrees 
axis([EptInt2(lmk)  EptInt2(end) -0.8 0.8]);              % RAW data 
text(EptInt2(ppthS),mmthS,num2str(Dgradthe2(ppthS))); 
text(timeRed(round(x1)),Dgradthe2(round(x1))/1000,['G ',num2str(GcMG2)]);  % Plotting CM degrees !!! 
subplot(2,2,3);imagesc(Mat1);hold;text(x1,y1,num2str(GcMG2));ylabel('Power (Au');

%%%%%%%%%% data for calculated on short
subplot(2,2,2);hold;plot(EptInt2,y2);   ylabel('Amp')   
subplot(2,2,2);plot(EptInt2,th2,'r') ;plot(EptInt2,DDgradRed/1000,'g')     % Theta and corresponding degrees 
text(EptInt2(ppthS),mmthS,num2str(Dgradthe2(ppthS))); 
text(timeRed(round(x1)),DDgradRed(round(x1))/1000,['G ',num2str(GcMG1)]);  % Plotting CM degrees !!! 
axis([EptInt2(lmk)  EptInt2(end) -0.8 0.8]); 
% subplot(2,1,2);figure; imagesc(tt22,f2,10*log10(Mat1*10e12));
figure;imagesc(Mat1);text(x1,y1,num2str(GcMG1));ylabel('Power (Au');
%   set(gca,'XTickLabel',num2str(DDgradRed ))
  %
  figure; imagesc(tt22,f2,10*log10(Mat1*10e12));
%   figure;imagesc(Mat1);text(x1,y1,'XX');ylabel('Power (Au');
%   set(gca,'XTickLabel',num2str(DDgradRed ))
%  pause;
% pause %(0.2)
close  all
% Mean Power for the different gamma subbands based on the wavelet Matrix
Matlow=Mat1(1:31,:);%fr=30-60
MatInter=Mat1(32:61,:);%fr=60-90
MatHigh=Mat1(62:121,:);%fr=91-125

if EmptyPeakTheta==0
   OutPut1{count1}=GcMG1;% Center of mass based on degress calculated  for short epochs
   OutPut2{count1}=GcMG2;% Center of mass based on degress calculated  for the 5 secs epoch
   OutPut3{count1}=Fre;       % Frequency of gamma in the center of mass


end

    end
end
%
     end
  end
 end
 % CM: Degrees vs Frequency
scatter(cell2mat(OutPut1),cell2mat(OutPut3));axis([0  360 30 150]);  hold; text(mean(cell2mat(OutPut2)),mean(cell2mat(OutPut3)),'O'), title([num2str(count1), ' Cicles / ',num2str(count2) ' s'])
edges=(-pi:0.3307:pi);
data2=deg2rad(cell2mat(OutPut2)-180);
% R109=[cell2mat(OutPut1);cell2mat(OutPut2);cell2mat(OutPut3)];
hist=histc(data2,edges); hist(end)=[];
% Rayleight test
[p,z]=circ_rtest((edges(1:end-1)+edges(2:end))/2,hist);  %% 
% mean circular 
[mu ul ll] = circ_mean((edges(1:end-1)+edges(2:end))/2,hist);  meanAn=rad2deg(mu)+180  ; text(meanAn,mean(cell2mat(OutPut3)),'X')
mvl=circ_r((edges(1:end-1)+edges(2:end))/2,hist);% mean vector length

mu=circ_rad2ang(mu)+180;
ul=circ_rad2ang(ul)+180;
ll=circ_rad2ang(ll)+180;
figure; hold
waveM=zeros(1,320);
 for i=1:length(thMean2);
    wav=cell2mat(thMean2(i));
    
    if   length(wav)>=321
        wav=wav(1,1:320);
       

    end
    plot(wav)
    waveM(i,1:end)=wav;% waveform cleaned from longer waveforms
     
end
%%%%%%%%%%%%%%%%%%%%%%%%%% MeanRAW
RawwaveM=zeros(1,320);

figure;hold
theR_mean=zeros(124,320,length(meany));
for   ii=1:length(meany);
    Rawmeany=cell2mat( meany(ii));
    if   length(Rawmeany)>=321
       Rawmeany=Rawmeany(1,1:320);
    end
%     subplot(2,1,1);plot(Rawmeany)  
    [theR,theR_t,theR_F] = nbt_TF_plot2( Rawmeany,0.1,30,150,1,0,Args);%[S2,t,f2] = nbt_TF_plot2(y,0.1,30,150,1,0,Args);
    RawwaveM(ii,1:end)=Rawmeany;
    theR=sqrt(abs(theR));theR=flipud(theR);
%     subplot(2,1,2);imagesc(theR_t,theR_F,theR)%imagesc(tt22,f2,Mat1)
    theR_mean(:,:,ii) =theR;
%     pause
   
end
figure;imagesc(theR_t,theR_F, mean(theR_mean,3))
Args3.MeanWavelete=mean(theR_mean,3);
% figure;imagesc(10*log10(Args2.MeanWavelete*10e12))

% Figure


%     

Args2.thetagamma=[cell2mat(OutPut1);cell2mat(OutPut2);cell2mat(OutPut3)]
Args2.ThetaMeanWave=thMean2;
Args2.CleanThetaWave=waveM;
% Args2.MeanWavelete=mean(theR_mean,3);
save(['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\LAB\R',num2str(Args.mouse),'_LFP\',num2str(Args.date),'-',num2str(Args.Mutant),'-',num2str(Args2.mouse),'_',Args2.ExP,num2str(Args2.cond),'_Pos'], 'Args3','Args2');
['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\LAB\R',num2str(Args.mouse),'_LFP\',num2str(Args.date),'-',num2str(Args.Mutant),'-',num2str(Args2.mouse),'_',Args2.ExP,num2str(Args2.cond),'_Pos2']
['C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\LAB\R',num2str(Args.mouse),'_LFP\',num2str(Args.date),'-',num2str(Args.Mutant),'-',num2str(Args.mouse),'_',Args.ExP,num2str(Args.ExP),'_Header']





% end



% pause
% close
% end

%  [Sc,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMatc(CA1',0.5,params)
%  
%  figure;imagesc(Sc');hold
% tstep