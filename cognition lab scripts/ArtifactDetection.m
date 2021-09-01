%(c) Brotons-Mas 2013 and ValElsa Corp
% This function, looks for artificats in the recording trace and correct
% them by detrending data above or below threshold. 
% LFP trace is broken in 0.5 sec epochs and the mean of the signal is calculated.
% If epoch mean is above 0.1 or  below -0.1 artifact is detected and
% detrended.


%%
 function [CleanDa,Appst1]=ArtifactDetection(Epoch,Fs,Plott)
% Artificat detection
Plott=0;
data=Epoch;
data = data';
sam=1:1:length(Epoch);
% ConfValmean=mean(abs(data))+(std(abs(data)));
% Peaks above 6 std are marked for possible artifact identification.
[Ammst1 Appst1] = findpeaks(data,'minpeakheight',mean(abs(data))+(6*std(abs(data))),'minpeakdistance',round(Fs*0.05),'sortstr','ascend'); % Peaks Above 20 Hz 
% [Ammst1 Appst1] = findpeaks(data,'minpeakheight',mean(abs(data))+(1*std(abs(data))),'sortstr','ascend'); % Peaks Above 20 Hz 

% SAppst1=sort(Appst1);
if Plott==1 % Trace is plotted if Plott=1
figure;subplot(3,1,1);plot(data);hold; text( Appst1,data(Appst1),'0')
subplot(3,1,2);plot(data);hold; 
subplot(3,1,3);hold; 
end

% This loop breake the signal into 1/2 second and calculates the need for
% artificat correction
for gbh=1:Fs/2:length(data)-Fs/2
    EpToTest=data(1,gbh:gbh+(Fs/2));    
%     EpToTest=data(gbh:gbh+(Fs/2),1);
    if Plott==1 % Trace is plotted if Plott=1
        figure;
        subplot(3,1,1);plot(data);hold; text( Appst1,data(Appst1),'X'), xlim([gbh gbh+(Fs/2)])
        subplot(3,1,2);plot(EpToTest)
%         subplot(3,1,2);plot(data);hold; 
%         subplot(3,1,3);hold; 
    end

    DetMValEpo=mean(EpToTest);
    DetMValEpo_abs=mean(abs(EpToTest));
   
    if DetMValEpo_abs>200 || DetMValEpo_abs<-200
       
       
        if Plott==1
        subplot(3,1,2);plot(sam(1,gbh:gbh+(Fs/2)),data(1,gbh:gbh+(Fs/2)),'r')
        end
            EpToTest  = detrend(EpToTest);% Best aproximation for detrending
%      RatTD(gbh)=FourierTransform(EpToTestD1,Fs,Plott) ; calculates
%      fourier transform
%          EpToTestD2= detrend(EpToTest,'constant'); % Different detrending
%          aproaches were tested. It seems than the linear one works well.
% %     RatTD(gbh)=FourierTransform(EpToTestD2,Fs,Plott)
%         plot(EpToTestD1,'r')
%         plot(EpToTestD2,'g')

        if Plott==1       
        subplot(3,1,2);plot(sam(1,gbh:gbh+(Fs/2)),data(1,gbh:gbh+(Fs/2)),'r') % Marks dirty epochs
        end
    end
 
    CleanDa(1,gbh:gbh+(Fs/2))=EpToTest; % concatenate the signal after detrending
end
    
  subplot(3,1,3);plot( CleanDa) % plots the clean signal
%% Saving file


  
  %   EpC=Count1;
% ConfValmean2=mean(abs(CleanDa))+(std(abs(CleanDa)));

% [Ammst1 Appst1] = findpeaks(data,'minpeakheight',mean(abs(data))+(6*std(abs(data))),'minpeakdistance',round(Fs*0.05),'sortstr','ascend'); % Peaks Above 20 Hz 
% 
% figure;subplot(3,1,1);plot(CleanDa);hold; text( Appst1,data(Appst1),'0')
% subplot(3,1,2);plot(CleanDa);hold; 
% subplot(3,1,3);hold; 
%  for gbh=1:Fs/2:length(CleanDa)-Fs/2;
%     EpToTest=CleanDa(1,gbh:gbh+(Fs/2));
% %     figure;plot(EpToTest);hold
%     
%     DetMValEpo2=mean(EpToTest);
%     if DetMValEpo2>0.1 || DetMValEpo2<-0.1
%         subplot(3,1,1);plot(sam(1,gbh:gbh+(Fs/2)),data(1,gbh:gbh+(Fs/2)),'r')
% %         y = erf(EpToTest); 
% %          bp= polyfit(EpToTest,y,(length(EpToTest)-1)/2);
%      EpToTest  = detrend(EpToTest);% Best aproximation for detrending
% %      RatTD(gbh)=FourierTransform(EpToTestD1,Fs,Plott)
% %          EpToTestD2= detrend(EpToTest,'constant');
% % %     RatTD(gbh)=FourierTransform(EpToTestD2,Fs,Plott)
% %         plot(EpToTestD1,'r')
% %         plot(EpToTestD2,'g')
%         subplot(3,1,2);plot(sam(1,gbh:gbh+(Fs/2)),CleanDa(1,gbh:gbh+(Fs/2)),'r')
% 
%     end
% % %     if mean(abs(EpToTest))>ConfValmean
% % %         conc1(1,gbh:gbh*Fs)=0;
% %     end
%     
% %     close all
%      CleanDa2(1,gbh:gbh+(Fs/2))=EpToTest;
%  end
%   subplot(3,1,3);plot( CleanDa2)

% figure;plot
% [xA1,yA1] = ginput(1);


% 
% for gbh=1:Fs/2:length(conc1)-Fs/2;
%     mean


% fs=Fs; %sampling rate
% resolHz=1; %spectral resolution in Hz
% nfft=floor(fs/resolHz); %number of points for FFT
% nw=7/2; %multitaper "time-bandwidth product" default:4  
% noverlap=0; %window overlap default:nfft/2
% delta=nfft-noverlap;
% ydata=conc1-mean(conc1); 
% %----calculating sliding power spectrum 
% for i=1:floor((length(conc1)-nfft)/(delta))
%     [Pxx,w] = pmtm(ydata((i-1)*delta+1:(i-1)*delta+nfft),nw,nfft,fs); 
%     areaPxx=1; %sum(Pxx);
%     PxxN(i,:)=Pxx(:)/areaPxx; %normalized power spectrum
% end;
% spec=10*log10(PxxN*10e12); %time-frequency spectrum
% nbin=i;
% time=((1:floor((length(conc1)-nfft)/(delta)))*delta).*(1/fs); %time variable
% 
% % 
% figure;plot(conc1);hold;
% InArtpp1=sort(Artpp1);
% InArt1=sort(Art1);
% text(InArtpp1,InArt1,'x')
% mean(abs(Art1))
%% Function that might be of interest
% mss = dspdata.msspectrum(Data)
% The mean-squared spectrum (MSS) is intended for discrete spectra.
% Unlike the power spectral density (PSD), the peaks in the MSS reflect
% the power in the signal at a given frequency. The MSS of a signal is 
% the Fourier transform of that signal's autocorrelation. 



%% e.g
% 
% t = 0:1e-3:1;
% x = sin(t*(10*2*pi)); % 10Hz
% 
% Fs = 1000;
% t = 0:1/Fs:2.96;
% x = cos(2*pi*t*1.24e3)+ cos(2*pi*t*10e3)+ randn(size(t));
% nfft = 2^nextpow2(length(x));
% Pxx = abs(fft(x,nfft)).^2/length(x)/Fs;
% Hpsd = dspdata.psd(Pxx,'Fs',Fs);   % Create PSD data object
% Hpsd = dspdata.msspectrum (Pxx,'Fs',Fs); % Create PSD data object and express power in DB
% figure; plot(Hpsd);      y = db2pow(Pxx)% decibels to power  The relationship between power and decibels is ydb = 10*log10(y).