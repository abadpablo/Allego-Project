function [frequencyBands] = bz_extractFrequencyBands(a,varargin)
% This function computes the power Spectrum and Coherence (accordingly to
% the input parameters) over different frequency bands.
% Power spectrum will be computed as the area under the curve (integral),
% whilst coherence will be computed ah the mean value for all the
% frequencies included in the band.
%
%
%   USAGE
%   
%   frequencyBands = bz_extractFrequencyBands(varargin)
%
%   % INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     FreqBands      
%     signal        PowerSpectrum or Coherence (default);
%     variable       Name of the variable to compute the frequency Bands
%                   (coherence_Shanks; coherogram; ...)
%
%
%    =========================================================================
%
%
% Created by Pablo Abad, 2021
%
%
%
%% Defaults and Params
p = inputParser;
addParameter(p,'FreqBands',[1 3; 4 12; 13 16; 17 29; 30 65; 66 130; 150 185],@isnumeric);
addParameter(p,'signal','coherence',@isstr);
addParameter(p,'variable',[],@isstr);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'analyzeSubSessions',false,@islogical);

parse(p,varargin{:});
FreqBands = p.Results.FreqBands;
signal = p.Results.signal;
variable = p.Results.variable;
analyzeSubSessions = p.Results.analyzeSubSessions;
basepath = p.Results.basepath;

nBands = size(FreqBands,1);

CV = 1;
if strcmpi(variable,'coherence_Shanks')
    if analyzeSubSessions
        if ~isempty(dir([basepath filesep '*MergePoints.events.mat']))
            disp('Loading MergePoints...')
            file = dir([basepath filesep '*MergePoints.events.mat']);
            load(file.name)
        end
        for ii=1:length(MergePoints.foldernames)
            coherogram = a.(MergePoints.foldernames{ii}).coherogram;
            phase = a.(MergePoints.foldernames{ii}).phase;
            t = a.(MergePoints.foldernames{ii}).t;
            f = a.(MergePoints.foldernames{ii}).f;
            S1 = a.(MergePoints.foldernames{ii}).S1;
            S2 = a.(MergePoints.foldernames{ii}).S2;
            ch1 = a.(MergePoints.foldernames{ii}).ch1;
            ch2 = a.(MergePoints.foldernames{ii}).ch2;
            coherence = [];
            % Power Spectrum
            for i=1:length(S1) 
                for j = 1:length(S1{i})
                    for it = 1:length(t)
                        for k=1:nBands
                            iDx = f >= FreqBands(k,1) & f < FreqBands(k,2);
                            t_Pow_S1 = S1{i}{j}(it,:);
                            t_Pow_S2 = S2{i}{j}(it,:);
                            Eband_S1{i}{j}(k,it) = CV*trapz(f(iDx),t_Pow_S1(iDx));
                            Eband_S2{i}{j}(k,it) = CV*trapz(f(iDx),t_Pow_S2(iDx));

                            c = coherogram{i}{j}(it,:);
                            ph = phase{i}{j}(it,:);
                            coherence{i}{j}(k,it) = mean(c(iDx));
                            Phase{i}{j}(k,it) = mean(ph(iDx));
                        end
                    end
                    Eband_S1_mean{i}{j} = mean(Eband_S1{i}{j},2);
                    Eband_S2_mean{i}{j} = mean(Eband_S2{i}{j},2);
                    phase_mean{i}{j} = mean(Phase{i}{j},2);
                    coherence_mean{i}{j} = mean(coherence{i}{j},2);
                end
            end
            frequencyBands.(MergePoints.foldernames{ii}).coherence_mean = coherence_mean;
            frequencyBands.(MergePoints.foldernames{ii}).phase_mean = phase_mean;
            frequencyBands.(MergePoints.foldernames{ii}).S1_mean = Eband_S1_mean;
            frequencyBands.(MergePoints.foldernames{ii}).S2_mean = Eband_S2_mean;

            frequencyBands.(MergePoints.foldernames{ii}).coherence = coherence;
            frequencyBands.(MergePoints.foldernames{ii}).phase = Phase;
            frequencyBands.(MergePoints.foldernames{ii}).S1 = Eband_S1;
            frequencyBands.(MergePoints.foldernames{ii}).S2 = Eband_S2;  
        end
        
    else
        coherogram = a.coherogram;
        phase = a.phase;
        S1 = a.S1;
        S2 = a.S2;
        t = a.t;
        f = a.f;

        coherence = [];    
        % Power Spectrum
        for i=1:length(S1) 
            for j = 1:length(S1{i})
                for it = 1:length(t)
                    for k=1:nBands
                        iDx = f >= FreqBands(k,1) & f < FreqBands(k,2);
                        t_Pow_S1 = S1{i}{j}(it,:);
                        t_Pow_S2 = S2{i}{j}(it,:);
                        Eband_S1{i}{j}(k,it) = CV*trapz(f(iDx),t_Pow_S1(iDx));
                        Eband_S2{i}{j}(k,it) = CV*trapz(f(iDx),t_Pow_S2(iDx));

                        c = coherogram{i}{j}(it,:);
                        ph = phase{i}{j}(it,:);
                        coherence{i}{j}(k,it) = mean(c(iDx));
                        Phase{i}{j}(k,it) = mean(ph(iDx));
                    end
                end
                Eband_S1_mean{i}{j} = mean(Eband_S1{i}{j},2);
                Eband_S2_mean{i}{j} = mean(Eband_S2{i}{j},2);
                phase_mean{i}{j} = mean(Phase{i}{j},2);
                coherence_mean{i}{j} = mean(coherence{i}{j},2);
            end
        end

        frequencyBands.coherence_mean = coherence_mean;
        frequencyBands.phase_mean = phase_mean;
        frequencyBands.S1_mean = Eband_S1_mean;
        frequencyBands.S2_mean = Eband_S2_mean;

        frequencyBands.coherence = coherence;
        frequencyBands.phase = Phase;
        frequencyBands.S1 = Eband_S1;
        frequencyBands.S2 = Eband_S2;  
    end

elseif strcmpi(variable,'coherogram')
    if analyzeSubSessions
        if ~isempty(dir([basepath filesep '*MergePoints.events.mat']))
            disp('Loading MergePoints...')
            file = dir([basepath filesep '*MergePoints.events.mat']);
            load(file.name)
        end
        
        for ii=1:length(MergePoints.foldernames)
            coherogram = a.(MergePoints.foldernames{ii}).coherogram;
            phase = a.(MergePoints.foldernames{ii}).phase;
            t = a.(MergePoints.foldernames{ii}).t;
            f = a.(MergePoints.foldernames{ii}).f;
            S1 = a.(MergePoints.foldernames{ii}).S1';
            S2 = a.(MergePoints.foldernames{ii}).S2';
            ch1 = a.(MergePoints.foldernames{ii}).ch1;
            ch2 = a.(MergePoints.foldernames{ii}).ch2;
            
            for it = 1:length(t)
                for k = 1:nBands
                    iDx = f >= FreqBands(k,1) & f < FreqBands(k,2);
                    t_Pow_S1 = S1(:,it);
                    t_Pow_S2 = S2(:,it);
                    Eband_S1(k,it) = CV*trapz(f(iDx),t_Pow_S1(iDx));
                    Eband_S2(k,it) = CV*trapz(f(iDx),t_Pow_S2(iDx));

                    c = coherogram(:,it);
                    ph = phase(:,it);
                    coherence(k,it) = mean(c(iDx));
                    Phase(k,it) = mean(ph(iDx));
                end
            end
            
            coherence_mean = mean(coherence,2);
            phase_mean = mean(Phase,2);
            Eband_S1_mean = mean(Eband_S1,2);
            Eband_S2_mean = mean(Eband_S2,2);

            frequencyBands.(MergePoints.foldernames{ii}).coherence_mean = coherence_mean;
            frequencyBands.(MergePoints.foldernames{ii}).phase_mean = phase_mean;
            frequencyBands.(MergePoints.foldernames{ii}).S1_mean = Eband_S1_mean;
            frequencyBands.(MergePoints.foldernames{ii}).S2_mean = Eband_S2_mean;

            frequencyBands.(MergePoints.foldernames{ii}).coherence = coherence;
            frequencyBands.(MergePoints.foldernames{ii}).phase = phase;
            frequencyBands.(MergePoints.foldernames{ii}).S1 = Eband_S1;
            frequencyBands.(MergePoints.foldernames{ii}).S2 = Eband_S2;

            frequencyBands.(MergePoints.foldernames{ii}).ch1 = ch1;
            frequencyBands.(MergePoints.foldernames{ii}).ch2 = ch2;

        end
    else
        coherogram = a.coherogram;
        phase = a.phase;
        t = a.t;
        f = a.f;
        S1 = a.S1';
        S2 = a.S2';
        ch1 = a.ch1;
        ch2 = a.ch2;

        for it = 1:length(t)
            for k = 1:nBands
                iDx = f >= FreqBands(k,1) & f < FreqBands(k,2);
                t_Pow_S1 = S1(:,it);
                t_Pow_S2 = S2(:,it);
                Eband_S1(k,it) = CV*trapz(f(iDx),t_Pow_S1(iDx));
                Eband_S2(k,it) = CV*trapz(f(iDx),t_Pow_S2(iDx));

                c = coherogram(:,it);
                ph = phase(:,it);
                coherence(k,it) = mean(c(iDx));
                Phase(k,it) = mean(ph(iDx));
            end
        end

        coherence_mean = mean(coherence,2);
        phase_mean = mean(Phase,2);
        Eband_S1_mean = mean(Eband_S1,2);
        Eband_S2_mean = mean(Eband_S2,2);

        frequencyBands.coherence_mean = coherence_mean;
        frequencyBands.phase_mean = phase_mean;
        frequencyBands.S1_mean = Eband_S1_mean;
        frequencyBands.S2_mean = Eband_S2_mean;

        frequencyBands.coherence = coherence;
        frequencyBands.phase = phase;
        frequencyBands.S1 = Eband_S1;
        frequencyBands.S2 = Eband_S2;

        frequencyBands.ch1 = ch1;
        frequencyBands.ch2 = ch2;        
    end
    
    
    
    
    
end


end

