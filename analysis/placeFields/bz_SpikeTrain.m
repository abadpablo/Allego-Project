
function [spikeTrain] = bz_SpikeTrain(spikes,varargin)
%spikeTrain - plot histograms of ISI distribution for a cluster
%
% Plots a histogram either representing the inter-spike interval (ISI)
% distribution for the selected spike events or the autocorrelation function
% in Hertz.  Also plotted is a gray area representing the "shadow"
% period during spike detection and a red area representing the user
% defined refractory period.
%
% The width of the histogram bins as well as the maximum time lag displayed
% are set by the following parameters:
%
%    spikes.params.display.isi_bin_size 
%    spikes.params.display.max_isi_to_display
%    spikes.params.display.correlations_bin_size 
%    spikes.params.display.max_autocorr_to_display
%    
% The user can switch between displaying ISIs or autocorrelatoin by right-
% clicking the axes and selecting from a context menu.  The choice can
% also be imposed on all isi plots in the same figure.
%
% On the y-axis is listed the number of refractory period violationss (RPVs)
% along with an estimated contamination percentage and its 95% confidence
% interval under the assumption that contaminating spikes are independent
% events. See poisson_contamination.m for more details.
%
% USAGE:
%       spikeTrain = bz_SpikeTrain( spikes, zoptions> )
%
%
% INPUTS:
%   spikes          a spikes structure
%   options         optional list of property-value pairs (see table below)
%   
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'show'        array describing which events to show in plot
%     'show_isi'    1-show ISI 0- show Autocorrelation
%                   default is set by spikes.params.display.show_isi 
%    =========================================================================
%
% OUTPUT
%
% spikeTrain.expected              expected spikes in the refractory period
% spikeTrain.lb                    superior limit
% spikeTrain.ub                    inferior
% spikeTrain.RPV                   spikes within the refractory period
% spikeTrain.meanF                 mean frequency
% spikeTrain.MeanAutoc             mean(autocorrelogram)
% spikeTrain.ProbBurst             % of spikes within the 10 msec window of the ISI, 
%                                  so % of spikes that could belong to a burst defined as spikes happening within 10 ms windows.
%                                  Staff, N. P., Jung, H. Y., Thiagarajan, T., Yao, M., & Spruston, N. (2000). Resting and
%                                  active properties of pyramidal neurons in subiculum and CA1 of rat hippocampus. J
%                                  Neurophysiol, 84(5), 2398-2408.
% spikeTrain.BurstIndexISI         This is a modification of the  BurstIndexA , it is a ratio obtained from the the ISI between  values in the burst intervals and values in the baseline , 
%                                  actually between 40-50 ms. This is based in a paper Royer et al 2012
% spikeTrain.BurstIndexA           This is the original  of the  BurstIndexA , it is a ratio obtained from the the ACH between  values in the burst intervals and values in the baseline , 
%                                  actually between 40-50 ms. This is based in a paper Royer et al 2012
% spikeTrain.num_burst             Number of doblets or triplets
% spikeTrain.n_spike               mean Number of spikes in each burst
% spikeTrain.prom_isi              mean inter spike interval in each burst 
% spikeTrain.prom_ibi              mean inter burst interval
% spikeTrain.prom_duracion         mean duration of a burst, 
% spikeTrain.numbursrNorm          Number of burst divided by total number of spikes
% spikeTrain.TMI                   This is a theta modulation index based in the autocorrelogram calculating a ratio between the through and the peak of the gamma modulation  mean value of the through  
%                                  ACH 50-70ms and ACH peak 100-124 ms  (TMI2 and TMI3) are modifications using the mean of the ACH or the sum of the intervals. 
%                                  The reason for using different approaches is that I don't have clear if the measurement is correct. The original TMI was used by Caccuci et al ... 
% spikeTrain.TMI2 
% spikeTrain.TMI3
% spikeTrain.S_theta               This is the theta modulation obtained  using the power
% spikeTrain.S_Gamma               The same for the gamma 
% spikeTrain.cross                 this is the cross correlogram for each cell 
% spikeTrain.S                     Power
% spikeTrain.f                     frequency


addpath('C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\OpenGSpike\QualityClustering\helpers')
addpath('C:\Users\Jorge\Documents\MATLAB\GlobAnalysis-1.0\OpenGSpike\QualityClustering\quality_measures')

%% Default Params
p = inputParser;
addParameter(p,'analyzeSubSessions',false,@islogical);
addParameter(p,'show','all',@isstr);
addParameter(p,'show_isi',1,@isnumeric);
addParameter(p,'isi_maxlag',1,@isnumeric);
addParameter(p,'shadow',1,@isnumeric);
addParameter(p,'refractory_period',1,@isnumeric);
addParameter(p,'bin_size',1,@isnmeric);
addParameter(p,'sr',30000,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'showWaveform',true,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'showFigure',false,@islogical);

parse(p,varargin{:});

show = p.Results.show;
show_isi = p.Results.show_isi;
isi_maxlag = p.Results.isi_maxlag;
shadow = p.Results.shadow;
refractory_period = p.Results.refractory_period;
bin_size = p.Results.bin_size;
sr = p.Results.sr;
saveMat = p.Results.saveMat;
basepath = (p.Results.basepath);
showWaveform = p.Results.showWaveform;
forceReload = p.Results.forceReload;
showFigure = p.Results.showFigure;
analyzeSubSessions = p.Results.analyzeSubSessions;

%% In case spikeTrain already exists 
if ~isempty(dir([basepath filesep '*spikeTrain.cellinfo.mat'])) || forceReload
    disp('spikeTrain already detected! Loading file.');
    file = dir([basepath filesep '*spikeTrain.cellinfo.mat']);
    load(file.name);
    return
end

if analyzeSubSessions
    if ~isempty(dir([basepath filesep '*spikeTrain.cellinfo.SubSession.mat'])) || forceReload
        disp('spikeTrain SubSession already detected! Loading file.');
        file = dir([basepath filesep '*spikeTrain.cellinfo.SubSession.mat']);
        load(file.name)
        return
    end
end


session = sessionTemplate(pwd,'showGUI',false);
duration = session.general.duration;
[sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);



if analyzeSubSessions
    try
        disp('Computing Spike Train Analysis for SubSessions...')
        if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
            load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
            count = 1;
            for ii = 1:size(MergePoints.foldernames,2)
                cd([basepath filesep MergePoints.foldernames{ii}]);
                fprintf('Computing Spike Train Analysis in %s folder \n',MergePoints.foldernames{ii});
                timestamps = MergePoints.timestamps(ii,:);
                for i=1:spikes.numcells
                    SpikeWaveform = spikes.filtWaveform{i};
                    spiketime = spikes.times{i};
                    
                    spiketime  = spiketime(timestamps(2) > spiketime & spiketime > timestamps(1)); 
                    
                    %% ELECTRO
                    if isi_maxlag >= 1
                        meanF = length(spiketime)/duration;
                        isis = diff(spiketime);
                        maxlag = isi_maxlag;
                        bins = bin_size;
                        [n,x] = hist(isis*1000,linspace(0,1000*maxlag,bins));
                        if isi_maxlag > 0.5
                            intervalnumber = length(spiketime)-1;
                            [MaxIsh PosMax] = max(n);
                            bins = round(1000* maxlag/bin_size );
                            % make plot
                            isis = isis(isis <= maxlag); 
                            [n,x] =hist(isis*1000,linspace(0,1000*maxlag,bins));
                            ymax = max(n)+1;
                            values1 = n(1:10);
                            values2 = n(40:50);
                            
                            ms10=mean(values1);
                            ms40=mean(values2);
                            if ms10>=ms40
                                BurstIndex=(ms10-ms40)/ms10;

                            else
                                BurstIndex=(ms10-ms40)/ms40;
                            end

                            Values1 = sum(values1);
                            ProbBurst = Values1*100/length(spiketime-1);

                            ll=length(spiketime);   

                            dif_max=0.015;
                            %%% se considera un burst si al menos esta compuesto por 2 spikes.
                            % n_spk=2;

                            %%%1.- calculo el vector de diferencias para determinar en donde ocurre un IBI
                            %%% o hay spikes que no pertenecen a un burst.
                            dif_spk = diff(spiketime);
                            no_burst = find(dif_spk > dif_max);

                            %%%% mi struct variable almacenara todos los datos que quiero por cada
                            %%%% burst encontrado en el registro, llamado reg_burst con campos:
                            %%%% nspk --> tiene el número de spikes por bursts
                            %%%% isi --> el ISI en cada burst
                            %%%% duration --> duracion del burst
                            %%%% burst --> el vector con los spikes que componen el burst.
                            cont=0;

                            ind=1;
                            if  (length(no_burst)+1 ~= length(spiketime))
                                aa = dif_spk;
                                aa(no_burst) = 0;
                                pp = find(aa);
                                if (no_burst(ind)-pp(1) >= 2)
                                    burst_spk=spiketime(pp(1):no_burst(ind));
                                    cont=cont+1;
                                    reg_burst(cont).nspk=length(burst_spk);
                                    reg_burst(cont).isi=mean(diff(burst_spk));
                                    reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                                    reg_burst(cont).burst=burst_spk;
                                    ind=ind+1; 
                                end
                                while (ind < length(no_burst))
                                    if (no_burst(ind+1)-no_burst(ind) > 1) %% este >1 no tiene que ver con el numero de spikes
                                        burst_spk = spiketime(no_burst(ind)+1:no_burst(ind+1));
                                        if (length(burst_spk)>=2)
                                             cont=cont+1;
                                             reg_burst(cont).nspk=length(burst_spk);
                                             reg_burst(cont).isi=mean(diff(burst_spk));
                                             reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                                             reg_burst(cont).burst=burst_spk;
                                        end
                                        ind=ind+1;
                                    else
                                        ind=ind+1;
                                    end
                                end

                                if (ind == length(no_burst))
                                    if (no_burst(ind) < length(spiketime)-1)
                                        burst_spk = spiketime(no_burst(ind)+1:length(spiketime));
                                        cont = cont+1;
                                        reg_burst(cont).nspk=length(burst_spk);
                                        reg_burst(cont).isi=mean(diff(burst_spk));
                                        reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                                        reg_burst(cont).burst=burst_spk;
                                    end
                                end
                                %%%%% Parte que calcula el IBI:
                                if (cont>2)
                                    for jj=1:cont-1
                                        t1=reg_burst(jj).burst(end);
                                        t2=reg_burst(jj+1).burst(1);
                                        dif_ibi(jj)=t2-t1;
                                    end
                                else
                                    dif_ibi=[];
                                end

                                %%%%%%% el promedio del numero de spikes por burst:
                                kk = [reg_burst.nspk];
                                %sprintf('numero promedio de spikes: %2.5g', mean(kk))
                                n_spike = mean(kk);
                                %%%% el promedio de ISI en cada burst es:
                                kk1 = [reg_burst.isi];
                                % sprintf('el promedio de ISI por el registro: %2.5g ms', mean(kk1)*1000)
                                prom_isi = mean(kk1)*1000;
                                %%%% la duracion:
                                kk2 = [reg_burst.duration];
                                % sprintf('la duracion promedio: %2.5g ms', mean(kk2)*1000)
                                prom_duracion = mean(kk2)*1000;
                                %%%% el IBI promedio:
                                % sprintf('el IBI promedio: %2.5g s', mean(dif_ibi))
                                if isempty(dif_ibi)
                                    prom_ibi=0;
                                else
                                    prom_ibi = mean(dif_ibi);
                                end
                                % numero de bursts
                                num_burst = cont;
                            else
                                %   sprintf('%s ','neurona no es burster')
                                num_burst=cont;
                                n_spike= 0;
                                prom_isi=0;
                                prom_duracion=0;
                                prom_ibi=0;
                            end
                            
                            numbursrNorm = (num_burst*100)/length(spiketime);
                            % XCorr

                            [C,C2,lags] = pxcorrV3(spiketime,spiketime, round(1000/bin_size), maxlag);
                            cross=C;
                            cross(1001)=0; % eliminates the peak for 0 bin difference in the ACH
                            autoco=cross;%/length(spiketimes); % Spikes normalized by the number of events, so probability !

                            MeanAutoc = mean(autoco(1001:end));
                            trough1 = mean(autoco(1051:1071)); % Original
                            % trough1=mean(autoco(1045:1081));
                            peak1 = mean(autoco(1101:1121));

                            trough11=sum(autoco(1045:1081));
                            peak11=sum(autoco(1101:1121));
                            TMI=abs(trough1-peak1)/(trough1+peak1);
                            TMI2=abs(MeanAutoc-peak1)/(MeanAutoc+peak1); 
                            TMI3=abs(trough11-peak11)/(trough11+peak11);

                            ms10A=mean(autoco(1,1002:1012));
                            ms40A=mean(autoco(1,1042:1052));
                            if ms10A>=ms40A
                                BurstIndexA = (ms10A-ms40A)/ms10A;
                            else
                                BurstIndexA = (ms10A-ms40A)/ms40A;
                            end

                            CrossSpike=autoco(1,502:1502);% To take the normal part of the autocorrelogram
                            CrossSpike = conv(CrossSpike, ones(1,3)/3);
                            CrossSpike = CrossSpike(500:end-1);% normalized by the total number of spikes and expressed as prob


                            params= struct('Fs',1000,'pad',1,'fpass',[0 120],'err',[1 0.05],'tapers',[3 5]) ;
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Multitaper%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %% Power Spectrum

                            [Sx,f1,Serr]=mtspectrumc((CrossSpike-(mean(CrossSpike)))',params); 
                            % figure;plot(f1,Sx);title([' tapers ', num2str(params.tapers)])%hold;plot(f,Serr)

                            %%%% Power modulation of spikes for theta and gamma 
                            iDx = f1 >= 4 & f1 < 12; % theta index %%% this index depend on the lenght of the ACH and need to be modified if the lenght changes
                            ind_th=find(iDx==1);

                            [The_Peak Th_PeakLoc]=max(Sx( iDx));
                            Th_Peak=f1 ( ind_th(Th_PeakLoc)); 

                            if  Th_PeakLoc==1
                                ThetaVec=[Sx(ind_th(Th_PeakLoc)) Sx(ind_th(Th_PeakLoc+1))]; 
                            elseif Th_PeakLoc>1 && Th_PeakLoc<=7
                               ThetaVec=[Sx(ind_th(Th_PeakLoc-1)) Sx(ind_th(Th_PeakLoc)) Sx(ind_th(Th_PeakLoc+1))]; 
                            elseif  Th_PeakLoc==8
                                ThetaVec=[Sx(ind_th(Th_PeakLoc-1)) Sx(ind_th(Th_PeakLoc)) ]; 
                            end

                            S_theta=mean( ThetaVec);
                            sum_theta= sum(Sx(ind_th))/sum(Sx);
                            in_BackG=find(f1>=100);
                            BackG=mean(Sx(in_BackG));
                            %     thet_Mod=S_theta/  BackG; % old
                            thet_Mod=sum (Sx(  iDx))/sum( Sx); % new 


                            BackG2ind=find(f1>=1 & f1<120);
                            BackG2=mean(Sx(BackG2ind));
                            
                            % Theta modulation was determined from the FFT-based power spectrum of the spike-train
                            % autocorrelation functions of individual cells. The bin size of the autocorrelogram was 2 ms. The
                            % autocorrelogram was truncated at 500 ms and the peak at zero lag was reduced to the maximal value
                            % not including the peak. The autocorrelation function was mean-normalized by subtracting the mean
                            % value from all values. The power spectrum was then generated. Before calculating the FFT from the
                            % autocorrelogram, the signal was tapered with a Hamming window to reduce spectral leakage. The
                            % length of the FFT was set to 216. The FFT was scaled to the length of the signal. The power
                            % spectrum was obtained by taking the square of the magnitude of the FFT. A cell was defined as
                            % theta modulated if the mean power within 1 Hz of each side of the peak in the 4-11 Hz frequency
                            % range was at least 5 times greater than the mean spectral power between 0 Hz and 125 Hz. 

                            % BackG=mean(Sx(1:find(f>100,1,'first')));

                            % S_theta=sum( Sx(ind_th))/sum(Sx(1:find(f>50,1,'first'))); % either the complete periodogram or up to 50 Hz

                            %%% Gamma
                            iDxContrast1 = find(f1 >0 & f1 < 30);
                            iDxContrast2 = find(f1 >100 & f1 < f1(end));
                            iDxContrast3=[ iDxContrast1 iDxContrast2];
                            iDxG1 = f1 >= 30 & f1 < 100; % Gamma index (notching by eliminating 50 hz (49-51))
                            ind_Gamma1=find(iDxG1==1);

                            [Ga_Peak Ga_PeakLoc]=max(Sx( iDxG1));
                            G_Peak = f1( ind_Gamma1(Ga_PeakLoc)); 

                            if   ind_Gamma1(Ga_PeakLoc)>(ind_Gamma1(1)) && ind_Gamma1(Ga_PeakLoc) ==103
                                GaVec = [Sx( ind_Gamma1(Ga_PeakLoc-1)) Sx( ind_Gamma1(Ga_PeakLoc)) ];
                            elseif ind_Gamma1(Ga_PeakLoc)==103
                                GaVec = [Sx( ind_Gamma1(Ga_PeakLoc)) Sx( ind_Gamma1(Ga_PeakLoc-1))];
                            elseif  ind_Gamma1(Ga_PeakLoc) >(ind_Gamma1(1)) && ind_Gamma1(Ga_PeakLoc)<103
                                GaVec = [Sx( ind_Gamma1(Ga_PeakLoc-1)) Sx( ind_Gamma1(Ga_PeakLoc)) Sx( ind_Gamma1(Ga_PeakLoc+1))];
                             elseif  ind_Gamma1(Ga_PeakLoc) ==(ind_Gamma1(1))  
                                GaVec = [ Sx( ind_Gamma1(Ga_PeakLoc)) Sx( ind_Gamma1(Ga_PeakLoc+1))];
                            end

                            S_Gamma=mean(GaVec);
                            S_Gamma2=sum(Sx(ind_Gamma1))/sum(Sx); 

                            S_Gamma2Contrast=sum(Sx(  iDxContrast1 ))/sum(Sx); 

                            gamma_Mod=sum (Sx(  iDx))/sum( Sx((  iDx==0))); % new 

                            GammaRat1=(S_Gamma-BackG)/(S_Gamma+BackG); 
                            GammaRat2=(S_Gamma-BackG2)/(S_Gamma+BackG2);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%FFT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%% Langston et al 2010

                            % Theta modulation was determined from the FFT-based power spectrum of the spike-train
                            % autocorrelation functions of individual cells. The bin size of the autocorrelogram was 2 ms. The
                            % autocorrelogram was truncated at 500 ms and the peak at zero lag was reduced to the maximal value
                            % not including the peak. The autocorrelation function was mean-normalized by subtracting the mean
                            % value from all values. The power spectrum was then generated. Before calculating the FFT from the
                            % autocorrelogram, the signal was tapered with a Hamming window to reduce spectral leakage. The
                            % length of the FFT was set to 216. The FFT was scaled to the length of the signal. The power
                            % spectrum was obtained by taking the square of the magnitude of the FFT. A cell was defined as
                            % theta modulated if the mean power within 1 Hz of each side of the peak in the 4-11 Hz frequency
                            % range was at least 5 times greater than the mean spectral power between 0 Hz and 125 Hz.
                            % 
                            Fs = 1000;                    % Sampling frequency
                            T = 1/Fs;                     % Sample time
                            L=length(CrossSpike);
                            NFFT = 2^nextpow2(L); % Next power of 2 from length of y
                            S = fft(CrossSpike-(mean(CrossSpike)),NFFT)/L; %%%% better after normalization
                            % need to test filtering by speed
                            S2 = abs(S).^2;
                            f2 = Fs/2*linspace(0,1,NFFT/2);
                            % figure
                            % % Plot single-sided amplitude spectrum.
                            % % plot(f2(1,1:125),2*abs(S2(1:125))) ;
                            % title(  ['nrcell=',num2str(Cell1),' Single-Sided Amplitude Spectrum of y(t) NO Hamming'])
                            % xlabel('Frequency (Hz)')
                            % ylabel('|Y(f)|')

                            %%%% Power modulation of spikes for theta and gamma 
                            iDx = f2 >= 4 & f2 < 12; % theta index
                            ind_th = find(iDx==1);


                            [The_Peak Th_PeakLoc]=max(S2( iDx));

                            Th_PeakFF=f2 ( ind_th(Th_PeakLoc)); 
                            if Th_PeakLoc==1
                                ThetaVec=[S2(ind_th(Th_PeakLoc)) S2(ind_th(Th_PeakLoc+1))];

                            elseif Th_PeakLoc>1 && Th_PeakLoc<4
                                ThetaVec=[S2(ind_th(Th_PeakLoc-1)) S2(ind_th(Th_PeakLoc)) S2(ind_th(Th_PeakLoc+1))];

                            elseif Th_PeakLoc==4
                                 ThetaVec=[S2(ind_th(Th_PeakLoc-1)) S2(ind_th(Th_PeakLoc))];

                            else
                                ThetaVec=[ S2(ind_th(Th_PeakLoc)) S2(ind_th(Th_PeakLoc+1))];
                            end

                            S_thetaFFT = mean( ThetaVec);

                            BackGFF = mean(S2(1:find(f2>100,1,'first')));
                            
                            %%___Gamma

                            iDxG1 = f2 >= 30 & f2 < 80; % Gamma index (notching by eliminating 50 hz (49-51))
                            ind_Gamma1=find(iDxG1==1);

                            [Ga_Peak Ga_PeakLoc]=max(S2( iDxG1));
                            G_PeakFF=f2 ( ind_Gamma1(Ga_PeakLoc)); 

                            if Ga_PeakLoc>1 && Ga_PeakLoc<25
                                GaVec=[S2( ind_Gamma1(Ga_PeakLoc-1)) S2( ind_Gamma1(Ga_PeakLoc)) S2( ind_Gamma1(Ga_PeakLoc+1))];
                            elseif Ga_PeakLoc==25
                                GaVec=[S2( ind_Gamma1(Ga_PeakLoc)) S2( ind_Gamma1(Ga_PeakLoc-1))];
                            else
                                GaVec=[S2( ind_Gamma1(Ga_PeakLoc)) S2( ind_Gamma1(Ga_PeakLoc+1))];
                            end

                            S_GammaFFT=mean(GaVec);
                            
                            %% MakPower
                            Makpower=0;
                            
                            if Makpower==1
                                [Eband Eband2 RatTD pm  fm] =MaKPower2(CrossSpike-(mean(CrossSpike)),Args,0,1)
                                % figure, plot(fm,pm)
                                %%%% Power modulation of spikes for theta and gamma 
                                %%%% Power modulation of spikes for theta and gamma 
                                iDxM = fm >= 4 & fm < 12; % theta index
                                ind_thM=find(iDxM==1);


                                [The_PeakM Th_PeakLocM]=max(pm( iDxM));

                                Th_PeakFFM=fm(ind_thM(Th_PeakLocM)); 
                                if Th_PeakLocM==1
                                    ThetaVecM=[pm(ind_thM,(Th_PeakLocM )) pm(ind_thM,(Th_PeakLocM +1))];     %ThetaVec=[S2(ind_th(Th_PeakLoc)) S2(ind_th(Th_PeakLoc+1))];

                                elseif Th_PeakLocM>1 && Th_PeakLocM<4
                                    ThetaVecM=[pm(ind_thM(Th_PeakLocM -1)) pm(ind_thM(Th_PeakLocM )) pm(ind_thM(Th_PeakLocM +1))];

                                elseif Th_PeakLocM ==4
                                    ThetaVecM=[pm(ind_thM(Th_PeakLocM -1)) pm(ind_thM(Th_PeakLocM ))];

                                else
                                   ThetaVecM=[ pm(ind_th(Th_PeakLocM )) pm(ind_th(Th_PeakLocM +1))];
                                end

                                S_thetaFFTM=mean( ThetaVecM);

                                BackGFTh_PeakLocFM=mean(pm(1:find(fm>100 & fm<150 ,1,'first')));

                                theta_BG_Ratio=S_thetaFFTM/BackGFTh_PeakLocFM
                                pause 
 

                                %%%___________________ Gamma

                                iDxG1 = f2 >= 30 & f2 < 80; % Gamma index (notching by eliminating 50 hz (49-51))
                                ind_Gamma1=find(iDxG1==1);

                               [Ga_Peak Ga_PeakLoc]=max(S( iDxG1));
                               G_PeakFF=f2 ( ind_Gamma1(Ga_PeakLoc)); 
                               if Ga_PeakLoc>1 && Ga_PeakLoc<51
                               GaVec=[S( ind_Gamma1(Ga_PeakLoc-1)) S( ind_Gamma1(Ga_PeakLoc)) S( ind_Gamma1(Ga_PeakLoc+1))];
                               elseif Ga_PeakLoc==51

                                   GaVec=[S( ind_Gamma1(Ga_PeakLoc)) S( ind_Gamma1(Ga_PeakLoc-1))];
                               else
                                      GaVec=[S( ind_Gamma1(Ga_PeakLoc)) S( ind_Gamma1(Ga_PeakLoc+1))];
                               end
                               S_GammaFFT=mean(GaVec);

                            end
                        
                    
                        end
                        
                    end
                    
                    %%  estimate poisson contamination for y-axis

                    [expected,lb,ub,RPV] = ss_rpv_contaminationV2(spiketime,spiketime,refractory_period,duration,shadow);
                    if isempty(lb)
                        data.ystr = [num2str(RPV) ' RPVs (' num2str(round(expected*100)) '%)' ];
                    else    
                        data.ystr = [num2str(RPV) ' RPVs (' num2str(round(lb*100)) '-' num2str(round(expected*100)) '-' num2str(round(ub*100)) '%)' ];
                    end
                    data.spiketimes = spiketime;
                    data.show_isi = show_isi;
                    data.isi_maxlag = isi_maxlag;
                    data.autocorr_maxlag = data.isi_maxlag;
                    data.shadow = shadow;
                    data.refractory_period = refractory_period;
                    data.corr_bin_size = bin_size;
                    data.isi_bin_size = bin_size;

                    set(gca,'UserData', data,'Tag','isi' )

                    % write updating method
                    if showFigure
                        update_isi( [], [], show_isi, gca);
                    end

                    spikeTrain.(MergePoints.foldernames{ii}){i}.expected = expected;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.lb = lb;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.ub = ub;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.RPV = RPV;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.meanF = meanF;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.MeanAutoc = MeanAutoc;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.ProbBurst = ProbBurst;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.BurstIndex = BurstIndex;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.BurstIndexA = BurstIndexA;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.num_burst = num_burst;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.n_spike = n_spike;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.prom_isi = prom_isi;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.prom_duration = prom_duracion;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.numBurstNorm = numbursrNorm;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.TMI = TMI;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.TMI2 = TMI2;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.TMI3 = TMI3;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.S_theta = S_theta;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.Th_Peak = Th_Peak;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.S_Gamma = S_Gamma;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.S_Gamma2 = S_Gamma2;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.G_Peak = G_Peak;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.thet_Mod = thet_Mod;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.BackG = BackG;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.S_thetaFFT = S_thetaFFT;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.Th_PeakFF = Th_PeakFF;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.S_GammaFFT = S_GammaFFT;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.G_PeakFF = G_PeakFF;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.BackGFF = BackGFF;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.cross = cross;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.Sx = Sx;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.f1 = f1;
                    spikeTrain.(MergePoints.foldernames{ii}){i}.spiketime = spiketime;
                    
                end
                
                cd(basepath)
                
            end
            
        end
            
    catch
        warning('It has not been possible to run Spike Train Analysis for SubSessions...')      
    end        
else
    try
        disp('Computing Spike Train Analysis for Whole Session..')
        for i=1:spikes.numcells
            
            SpikeWaveform = spikes.filtWaveform{i};
            spiketime = spikes.times{i};

            %% ELECTRO
            if isi_maxlag >=1
                meanF = length(spiketime)/duration;
                isis = diff(spiketime); 
                maxlag = isi_maxlag;
                bins = bin_size;
                [n,x] = hist(isis*1000,linspace(0,1000*maxlag,bins)); 
                if isi_maxlag > 0.5
                    intervalnumber = length(spiketime)-1;
                    [MaxIsh PosMax]=max(n);
                    maxlag = isi_maxlag;
                    bins = round(1000* maxlag/bin_size );

                    % make plot
                    isis = isis(isis <= maxlag); 
                    [n,x] =hist(isis*1000,linspace(0,1000*maxlag,bins));
                    ymax = max(n)+1;
                    values1 = n(1:10);
                    values2 = n(40:50);

                    ms10=mean(values1);
                    ms40=mean(values2);
                    if ms10>=ms40
                        BurstIndex=(ms10-ms40)/ms10;

                    else
                        BurstIndex=(ms10-ms40)/ms40;
                    end

                    Values1=sum(values1);
                    ProbBurst = Values1*100/length(spiketime-1);

                    ll=length(spiketime);   

                    dif_max=0.015; 

                    %%% se considera un burst si al menos esta compuesto por 2 spikes.
                    % n_spk=2;

                    %%%1.- calculo el vector de diferencias para determinar en donde ocurre un IBI
                    %%% o hay spikes que no pertenecen a un burst.
                    dif_spk = diff(spiketime);
                    no_burst = find(dif_spk > dif_max);

                    %%%% mi struct variable almacenara todos los datos que quiero por cada
                    %%%% burst encontrado en el registro, llamado reg_burst con campos:
                    %%%% nspk --> tiene el número de spikes por bursts
                    %%%% isi --> el ISI en cada burst
                    %%%% duration --> duracion del burst
                    %%%% burst --> el vector con los spikes que componen el burst.
                    cont=0;

                    ind=1;
                    if  (length(no_burst)+1 ~= length(spiketime))
                        aa = dif_spk;
                        aa(no_burst) = 0;
                        pp = find(aa);
                        if (no_burst(ind)-pp(1) >= 2)
                            burst_spk=spiketimes(pp(1):no_burst(ind));
                            cont=cont+1;
                            reg_burst(cont).nspk=length(burst_spk);
                            reg_burst(cont).isi=mean(diff(burst_spk));
                            reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                            reg_burst(cont).burst=burst_spk;
                            ind=ind+1; 
                        end
                        while (ind < length(no_burst))
                            if (no_burst(ind+1)-no_burst(ind) > 1) %% este >1 no tiene que ver con el numero de spikes
                                burst_spk = spiketime(no_burst(ind)+1:no_burst(ind+1));
                                if (length(burst_spk)>=2)
                                     cont=cont+1;
                                     reg_burst(cont).nspk=length(burst_spk);
                                     reg_burst(cont).isi=mean(diff(burst_spk));
                                     reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                                     reg_burst(cont).burst=burst_spk;
                                end
                                ind=ind+1;
                            else
                                ind=ind+1;
                            end
                        end
                        
                        if (ind == length(no_burst))
                            if (no_burst(ind) < length(spiketime)-1)
                                burst_spk = spiketime(no_burst(ind)+1:length(spiketime));
                                cont = cont+1;
                                reg_burst(cont).nspk=length(burst_spk);
                                reg_burst(cont).isi=mean(diff(burst_spk));
                                reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                                reg_burst(cont).burst=burst_spk;
                            end
                        end
                        %%%%% Parte que calcula el IBI:
                        if (cont>2)
                            for jj=1:cont-1
                                t1=reg_burst(jj).burst(end);
                                t2=reg_burst(jj+1).burst(1);
                                dif_ibi(jj)=t2-t1;
                            end
                        else
                            dif_ibi=[];
                        end
                        
                        %%%%%%% el promedio del numero de spikes por burst:
                        kk = [reg_burst.nspk];
                        %sprintf('numero promedio de spikes: %2.5g', mean(kk))
                        n_spike = mean(kk);
                        %%%% el promedio de ISI en cada burst es:
                        kk1 = [reg_burst.isi];
                        % sprintf('el promedio de ISI por el registro: %2.5g ms', mean(kk1)*1000)
                        prom_isi = mean(kk1)*1000;
                        %%%% la duracion:
                        kk2 = [reg_burst.duration];
                        % sprintf('la duracion promedio: %2.5g ms', mean(kk2)*1000)
                        prom_duracion = mean(kk2)*1000;
                        %%%% el IBI promedio:
                        % sprintf('el IBI promedio: %2.5g s', mean(dif_ibi))
                        if isempty(dif_ibi)
                            prom_ibi=0;
                        else
                            prom_ibi = mean(dif_ibi);
                        end
                        % numero de bursts
                        num_burst = cont;
                    else
                        %   sprintf('%s ','neurona no es burster')
                        num_burst=cont;
                        n_spike= 0;
                        prom_isi=0;
                        prom_duracion=0;
                        prom_ibi=0;
                    end


                    numbursrNorm = (num_burst*100)/length(spiketime);
                    % XCorr

                    [C,C2,lags] = pxcorrV3(spiketime,spiketime, round(1000/bin_size), maxlag);
                    cross=C;
                    cross(1001)=0; % eliminates the peak for 0 bin difference in the ACH
                    autoco=cross;%/length(spiketimes); % Spikes normalized by the number of events, so probability !

                    MeanAutoc = mean(autoco(1001:end));
                    trough1 = mean(autoco(1051:1071)); % Original
                    % trough1=mean(autoco(1045:1081));
                    peak1 = mean(autoco(1101:1121));

                    trough11=sum(autoco(1045:1081));
                    peak11=sum(autoco(1101:1121));
                    TMI=abs(trough1-peak1)/(trough1+peak1);
                    TMI2=abs(MeanAutoc-peak1)/(MeanAutoc+peak1); 
                    TMI3=abs(trough11-peak11)/(trough11+peak11);

                    ms10A=mean(autoco(1,1002:1012));
                    ms40A=mean(autoco(1,1042:1052));
                    if ms10A>=ms40A
                        BurstIndexA = (ms10A-ms40A)/ms10A;
                    else
                        BurstIndexA = (ms10A-ms40A)/ms40A;
                    end

                    CrossSpike=autoco(1,502:1502);% To take the normal part of the autocorrelogram
                    CrossSpike = conv(CrossSpike, ones(1,3)/3);
                    CrossSpike = CrossSpike(500:end-1);% normalized by the total number of spikes and expressed as prob


                    params= struct('Fs',1000,'pad',1,'fpass',[0 120],'err',[1 0.05],'tapers',[3 5]) ;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Multitaper%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% Power Spectrum

                    [Sx,f1,Serr]=mtspectrumc((CrossSpike-(mean(CrossSpike)))',params); 
                    % figure;plot(f1,Sx);title([' tapers ', num2str(params.tapers)])%hold;plot(f,Serr)

                    %%%% Power modulation of spikes for theta and gamma 
                    iDx = f1 >= 4 & f1 < 12; % theta index %%% this index depend on the lenght of the ACH and need to be modified if the lenght changes
                    ind_th=find(iDx==1);

                    [The_Peak Th_PeakLoc]=max(Sx( iDx));
                    Th_Peak=f1 ( ind_th(Th_PeakLoc)); 

                    if  Th_PeakLoc==1
                        ThetaVec=[Sx(ind_th(Th_PeakLoc)) Sx(ind_th(Th_PeakLoc+1))]; 
                    elseif Th_PeakLoc>1 && Th_PeakLoc<=7
                       ThetaVec=[Sx(ind_th(Th_PeakLoc-1)) Sx(ind_th(Th_PeakLoc)) Sx(ind_th(Th_PeakLoc+1))]; 
                    elseif  Th_PeakLoc==8
                        ThetaVec=[Sx(ind_th(Th_PeakLoc-1)) Sx(ind_th(Th_PeakLoc)) ]; 
                    end

                    S_theta=mean( ThetaVec);
                    sum_theta= sum(Sx(ind_th))/sum(Sx);
                    in_BackG=find(f1>=100);
                    BackG=mean(Sx(in_BackG));
                    %     thet_Mod=S_theta/  BackG; % old
                    thet_Mod=sum (Sx(  iDx))/sum( Sx); % new 


                    BackG2ind=find(f1>=1 & f1<120);
                    BackG2=mean(Sx(BackG2ind));


                    % Theta modulation was determined from the FFT-based power spectrum of the spike-train
                    % autocorrelation functions of individual cells. The bin size of the autocorrelogram was 2 ms. The
                    % autocorrelogram was truncated at 500 ms and the peak at zero lag was reduced to the maximal value
                    % not including the peak. The autocorrelation function was mean-normalized by subtracting the mean
                    % value from all values. The power spectrum was then generated. Before calculating the FFT from the
                    % autocorrelogram, the signal was tapered with a Hamming window to reduce spectral leakage. The
                    % length of the FFT was set to 216. The FFT was scaled to the length of the signal. The power
                    % spectrum was obtained by taking the square of the magnitude of the FFT. A cell was defined as
                    % theta modulated if the mean power within 1 Hz of each side of the peak in the 4-11 Hz frequency
                    % range was at least 5 times greater than the mean spectral power between 0 Hz and 125 Hz. 

                    % BackG=mean(Sx(1:find(f>100,1,'first')));

                    % S_theta=sum( Sx(ind_th))/sum(Sx(1:find(f>50,1,'first'))); % either the complete periodogram or up to 50 Hz

                    %%% Gamma
                    iDxContrast1 = find(f1 >0 & f1 < 30);
                    iDxContrast2 = find(f1 >100 & f1 < f1(end));
                    iDxContrast3=[ iDxContrast1 iDxContrast2];
                    iDxG1 = f1 >= 30 & f1 < 100; % Gamma index (notching by eliminating 50 hz (49-51))
                    ind_Gamma1=find(iDxG1==1);

                    [Ga_Peak Ga_PeakLoc]=max(Sx( iDxG1));
                    G_Peak = f1( ind_Gamma1(Ga_PeakLoc)); 

                    if   ind_Gamma1(Ga_PeakLoc)>(ind_Gamma1(1)) && ind_Gamma1(Ga_PeakLoc) ==103
                        GaVec = [Sx( ind_Gamma1(Ga_PeakLoc-1)) Sx( ind_Gamma1(Ga_PeakLoc)) ];
                    elseif ind_Gamma1(Ga_PeakLoc)==103
                        GaVec = [Sx( ind_Gamma1(Ga_PeakLoc)) Sx( ind_Gamma1(Ga_PeakLoc-1))];
                    elseif  ind_Gamma1(Ga_PeakLoc) >(ind_Gamma1(1)) && ind_Gamma1(Ga_PeakLoc)<103
                        GaVec = [Sx( ind_Gamma1(Ga_PeakLoc-1)) Sx( ind_Gamma1(Ga_PeakLoc)) Sx( ind_Gamma1(Ga_PeakLoc+1))];
                     elseif  ind_Gamma1(Ga_PeakLoc) ==(ind_Gamma1(1))  
                        GaVec = [ Sx( ind_Gamma1(Ga_PeakLoc)) Sx( ind_Gamma1(Ga_PeakLoc+1))];
                    end

                    S_Gamma=mean(GaVec);
                    S_Gamma2=sum(Sx(ind_Gamma1))/sum(Sx); 

                    S_Gamma2Contrast=sum(Sx(  iDxContrast1 ))/sum(Sx); 

                    gamma_Mod=sum (Sx(  iDx))/sum( Sx((  iDx==0))); % new 

                    GammaRat1=(S_Gamma-BackG)/(S_Gamma+BackG); 
                    GammaRat2=(S_Gamma-BackG2)/(S_Gamma+BackG2);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%FFT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%% Langston et al 2010

                    % Theta modulation was determined from the FFT-based power spectrum of the spike-train
                    % autocorrelation functions of individual cells. The bin size of the autocorrelogram was 2 ms. The
                    % autocorrelogram was truncated at 500 ms and the peak at zero lag was reduced to the maximal value
                    % not including the peak. The autocorrelation function was mean-normalized by subtracting the mean
                    % value from all values. The power spectrum was then generated. Before calculating the FFT from the
                    % autocorrelogram, the signal was tapered with a Hamming window to reduce spectral leakage. The
                    % length of the FFT was set to 216. The FFT was scaled to the length of the signal. The power
                    % spectrum was obtained by taking the square of the magnitude of the FFT. A cell was defined as
                    % theta modulated if the mean power within 1 Hz of each side of the peak in the 4-11 Hz frequency
                    % range was at least 5 times greater than the mean spectral power between 0 Hz and 125 Hz.
                    % 
                    Fs = 1000;                    % Sampling frequency
                    T = 1/Fs;                     % Sample time
                    L=length(CrossSpike);
                    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
                    S = fft(CrossSpike-(mean(CrossSpike)),NFFT)/L; %%%% better after normalization
                    % need to test filtering by speed
                    S2 = abs(S).^2;
                    f2 = Fs/2*linspace(0,1,NFFT/2);
                    % figure
                    % % Plot single-sided amplitude spectrum.
                    % % plot(f2(1,1:125),2*abs(S2(1:125))) ;
                    % title(  ['nrcell=',num2str(Cell1),' Single-Sided Amplitude Spectrum of y(t) NO Hamming'])
                    % xlabel('Frequency (Hz)')
                    % ylabel('|Y(f)|')

                    %%%% Power modulation of spikes for theta and gamma 
                    iDx = f2 >= 4 & f2 < 12; % theta index
                    ind_th = find(iDx==1);


                    [The_Peak Th_PeakLoc]=max(S2( iDx));

                    Th_PeakFF=f2 ( ind_th(Th_PeakLoc)); 
                    if Th_PeakLoc==1
                        ThetaVec=[S2(ind_th(Th_PeakLoc)) S2(ind_th(Th_PeakLoc+1))];

                    elseif Th_PeakLoc>1 && Th_PeakLoc<4
                        ThetaVec=[S2(ind_th(Th_PeakLoc-1)) S2(ind_th(Th_PeakLoc)) S2(ind_th(Th_PeakLoc+1))];

                    elseif Th_PeakLoc==4
                         ThetaVec=[S2(ind_th(Th_PeakLoc-1)) S2(ind_th(Th_PeakLoc))];

                    else
                        ThetaVec=[ S2(ind_th(Th_PeakLoc)) S2(ind_th(Th_PeakLoc+1))];
                    end

                    S_thetaFFT = mean( ThetaVec);

                    BackGFF = mean(S2(1:find(f2>100,1,'first')));

      %     S_theta=sum( S2(ind_th))/sum(S2(1:find(f>50,1,'first'))); % either the complete periodogram or up to 50 Hz

      %%%___________________ Gamma

                iDxG1 = f2 >= 30 & f2 < 80; % Gamma index (notching by eliminating 50 hz (49-51))
                ind_Gamma1=find(iDxG1==1);

                [Ga_Peak Ga_PeakLoc]=max(S2( iDxG1));
                G_PeakFF=f2 ( ind_Gamma1(Ga_PeakLoc)); 

                if Ga_PeakLoc>1 && Ga_PeakLoc<25
                    GaVec=[S2( ind_Gamma1(Ga_PeakLoc-1)) S2( ind_Gamma1(Ga_PeakLoc)) S2( ind_Gamma1(Ga_PeakLoc+1))];
                elseif Ga_PeakLoc==25
                    GaVec=[S2( ind_Gamma1(Ga_PeakLoc)) S2( ind_Gamma1(Ga_PeakLoc-1))];
                else
                    GaVec=[S2( ind_Gamma1(Ga_PeakLoc)) S2( ind_Gamma1(Ga_PeakLoc+1))];
                end
                
                S_GammaFFT=mean(GaVec);

                %% MakPower
                Makpower=0;
                if Makpower==1
                    [Eband Eband2 RatTD pm  fm] =MaKPower2(CrossSpike-(mean(CrossSpike)),Args,0,1)
                    % figure, plot(fm,pm)
                    %%%% Power modulation of spikes for theta and gamma 
                    %%%% Power modulation of spikes for theta and gamma 
                    iDxM = fm >= 4 & fm < 12; % theta index
                    ind_thM=find(iDxM==1);


                    [The_PeakM Th_PeakLocM]=max(pm( iDxM));

                    Th_PeakFFM=fm(ind_thM(Th_PeakLocM)); 
                    if Th_PeakLocM==1
                        ThetaVecM=[pm(ind_thM,(Th_PeakLocM )) pm(ind_thM,(Th_PeakLocM +1))];     %ThetaVec=[S2(ind_th(Th_PeakLoc)) S2(ind_th(Th_PeakLoc+1))];

                    elseif Th_PeakLocM>1 && Th_PeakLocM<4
                        ThetaVecM=[pm(ind_thM(Th_PeakLocM -1)) pm(ind_thM(Th_PeakLocM )) pm(ind_thM(Th_PeakLocM +1))];

                    elseif Th_PeakLocM ==4
                        ThetaVecM=[pm(ind_thM(Th_PeakLocM -1)) pm(ind_thM(Th_PeakLocM ))];

                    else
                       ThetaVecM=[ pm(ind_th(Th_PeakLocM )) pm(ind_th(Th_PeakLocM +1))];
                    end

                    S_thetaFFTM=mean( ThetaVecM);

                    BackGFTh_PeakLocFM=mean(pm(1:find(fm>100 & fm<150 ,1,'first')));

                    theta_BG_Ratio=S_thetaFFTM/BackGFTh_PeakLocFM
                    pause 
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%
          %     S_theta=sum( S(ind_th))/sum(S(1:find(f>50,1,'first'))); % either the complete periodogram or up to 50 Hz

          %%%___________________ Gamma

                    iDxG1 = f2 >= 30 & f2 < 80; % Gamma index (notching by eliminating 50 hz (49-51))
                    ind_Gamma1=find(iDxG1==1);

                   [Ga_Peak Ga_PeakLoc]=max(S( iDxG1));
                   G_PeakFF=f2 ( ind_Gamma1(Ga_PeakLoc)); 
                   if Ga_PeakLoc>1 && Ga_PeakLoc<51
                   GaVec=[S( ind_Gamma1(Ga_PeakLoc-1)) S( ind_Gamma1(Ga_PeakLoc)) S( ind_Gamma1(Ga_PeakLoc+1))];
                   elseif Ga_PeakLoc==51

                       GaVec=[S( ind_Gamma1(Ga_PeakLoc)) S( ind_Gamma1(Ga_PeakLoc-1))];
                   else
                          GaVec=[S( ind_Gamma1(Ga_PeakLoc)) S( ind_Gamma1(Ga_PeakLoc+1))];
                   end
                   S_GammaFFT=mean(GaVec);

            end
        end
    end

    %%  estimate poisson contamination for y-axis

        [expected,lb,ub,RPV] = ss_rpv_contaminationV2(spiketime,spiketime,refractory_period,duration,shadow);
        if isempty(lb)
            data.ystr = [num2str(RPV) ' RPVs (' num2str(round(expected*100)) '%)' ];
        else    
            data.ystr = [num2str(RPV) ' RPVs (' num2str(round(lb*100)) '-' num2str(round(expected*100)) '-' num2str(round(ub*100)) '%)' ];
        end
        data.spiketimes = spiketime;
        data.show_isi = show_isi;
        data.isi_maxlag = isi_maxlag;
        data.autocorr_maxlag = data.isi_maxlag;
        data.shadow = shadow;
        data.refractory_period = refractory_period;
        data.corr_bin_size = bin_size;
        data.isi_bin_size = bin_size;

        set(gca,'UserData', data,'Tag','isi' )

        % write updating method
        if showFigure
            update_isi( [], [], show_isi, gca);
        end

        spikeTrain{i}.expected = expected;
        spikeTrain{i}.lb = lb;
        spikeTrain{i}.ub = ub;
        spikeTrain{i}.RPV = RPV;
        spikeTrain{i}.meanF = meanF;
        spikeTrain{i}.MeanAutoc = MeanAutoc;
        spikeTrain{i}.ProbBurst = ProbBurst;
        spikeTrain{i}.BurstIndex = BurstIndex;
        spikeTrain{i}.BurstIndexA = BurstIndexA;
        spikeTrain{i}.num_burst = num_burst;
        spikeTrain{i}.n_spike = n_spike;
        spikeTrain{i}.prom_isi = prom_isi;
        spikeTrain{i}.prom_duration = prom_duracion;
        spikeTrain{i}.numBurstNorm = numbursrNorm;
        spikeTrain{i}.TMI = TMI;
        spikeTrain{i}.TMI2 = TMI2;
        spikeTrain{i}.TMI3 = TMI3;
        spikeTrain{i}.S_theta = S_theta;
        spikeTrain{i}.Th_Peak = Th_Peak;
        spikeTrain{i}.S_Gamma = S_Gamma;
        spikeTrain{i}.S_Gamma2 = S_Gamma2;
        spikeTrain{i}.G_Peak = G_Peak;
        spikeTrain{i}.thet_Mod = thet_Mod;
        spikeTrain{i}.BackG = BackG;
        spikeTrain{i}.S_thetaFFT = S_thetaFFT;
        spikeTrain{i}.Th_PeakFF = Th_PeakFF;
        spikeTrain{i}.S_GammaFFT = S_GammaFFT;
        spikeTrain{i}.G_PeakFF = G_PeakFF;
        spikeTrain{i}.BackGFF = BackGFF;
        spikeTrain{i}.cross = cross;
        spikeTrain{i}.Sx = Sx;
        spikeTrain{i}.f1 = f1;
        spikeTrain{i}.spiketime = spiketime;   
        end
    catch
        warning('It has not been possible to run Spike Train Analysis for Whole Session ...')
    end
end


    % Plot Waveform

    % save mat file
    if saveMat
        basename = session.general.name;
        if analyzeSubSessions
            save(fullfile(basepath,[basename,'.spikeTrain.cellinfo.SubSession.mat']),'spikeTrain');
        else    
            save(fullfile(basepath,[basename,'.spikeTrain.cellinfo.mat']),'spikeTrain');
        end
    end

end








% callback that updates display on axes
function update_isi( hObject, event, displaymode, ax)

        % get display mode
        data = get( ax,'UserData');
        data.show_isi = displaymode;
        set(ax,'UserData',data)
        set(gcf,'CurrentAxes',ax )
        
        % update display
        make_isi;
    
        % set context menu - allow switch between ISI and autocorrelation
        % as well as a global switch for all plot_isi instances on the current figure
        cmenu = uicontextmenu;
        item(1) = uimenu(cmenu, 'Label', 'Show ISI', 'Callback', {@update_isi, 1,ax} );
        item(2) = uimenu(cmenu, 'Label', 'Show autocorr', 'Callback',{@update_isi, 0,ax}  );
        item(3) = uimenu(cmenu, 'Label', 'Use this style on all ISIs in figure', 'Callback',{@impose_all,displaymode,ax},'Separator','on'  );
        set(item(2-displaymode), 'Checked', 'on');    
        set(ax,'UIContextMenu', cmenu )

end

% callback to impose display mode on all plot_isi axes in this figure
function impose_all(hObject,event,displaymode,ax)
        [o,h] = gcbo;
        my_axes = findobj(h,'Tag','isi');
        my_axes = setdiff(my_axes,ax);
        for j = 1:length(my_axes), update_isi([],[],displaymode,my_axes(j)); end
end    
    
% plots the ISI or autocorrelation
function make_isi     

    data = get( gca,'UserData');
    spiketimes = data.spiketimes;
    shadow = data.shadow;
    rp     = data.refractory_period;
    cla reset
    
    % ISI case
    if data.show_isi
        maxlag = data.isi_maxlag;
        bins = round(1000* maxlag/data.isi_bin_size );

        % make plot
        isis = diff(spiketimes);   
        isis = isis(isis <= maxlag); 
        [n,x] =hist(isis*1000,linspace(0,1000*maxlag,bins));
        ymax   = max(n)+1;
%       values1=n(1:10);
% Values1=sum(values1);
% ProbBurst=Values1*100/length(spiketimes-1);
%%        
        
        % make patches to represent shadow and refractory period

        patch([0 shadow shadow 0 ], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none');
        patch([shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        hold on,    b2 = bar(x,n,1.0); hold off
        set(b2,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])

        % update axes
        set(gca,'YLim',[0 ymax],'XLim',[0 1000*maxlag])
        xlabel('Interspike interval (msec)')
       ylabel({'No. of spikes',data.ystr})
 
    else
        maxlag = data.autocorr_maxlag;
        
        % calculate autocorrelation
        if length(spiketimes) > 1
          [cross,lags] = pxcorr(spiketimes,spiketimes, round(1000/data.corr_bin_size), maxlag);
        else
            cross = 0;  lags = 0;
        end
        cross(find(lags==0)) = 0;
        
        % place patches to represent shadow and refractory period
        ymax = max(cross) + 1;
        patch(shadow*[-1 1 1 -1], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none');
        patch([shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        patch(-[shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        
        % plot autocorrelation histogram
        hold on, bb = bar(lags*1000,cross,1.0); hold off;  
        set(bb,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])

        % set axes
        set(gca, 'XLim', maxlag*1000*[-1 1]);
        set(gca,'YLim',[0 ymax])
        xlabel( 'Time lag (msec)')
        ylabel({'Autocorrelation (Hz)',data.ystr})
        
    end
    set(gca,'Tag','isi','UserData',data)
     
end

% if   maxlag=>1
% %         lb1=0 
% %         ub1=0 
% %         RPV1=0 
%         meanF1=0
%         MeanAutoc1=0
%         ProbBurst1=0
%         num_burst1=0
%         n_spike1=0
%         prom_isi1=0
%         prom_ibi1=0
%         prom_duracion1 =0
%         numbursrNorm1 =0
%     TMI1=0
%     end

