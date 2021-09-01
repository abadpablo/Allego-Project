function [data] = removeArtifacts(varargin)

% [data] = sessionRemoveArtifacts(varargin)

% Find electrical artifacts in the session produced when animal drinks from needle in
% LinearMaze and TMaze

% INPUTS
% 
% <OPTIONALS>
% filename
% basepath
% correctDC - Logical variable to indicate if DC is corrected, default
%               false
% ch        - List of channels to clean pulses, default all
% winArt    - Window for artefact removal, in seconds, default 0.0005s
% winDC     - Window for DC removal, in seconds, default 0.005s
% 
% OUTPUTS
% data
% 
% Pablo Abad 2021


% Default parameters
filename = split(pwd,'\'); filename = filename{end};

% Parse options
p = inputParser;
addParameter(p,'filename',filename,@isstr);
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'correctDC',false, @islogical);
addParameter(p,'ch','all');
addParameter(p,'winArt',0.0005,@isnumeric);
addParameter(p,'winDC',0.005,@isnumeric);

parse(p, varargin{:});

filename = p.Results.filename;
basepath = p.Results.basepath;
correctDC = p.Results.correctDC;
ch = p.Results.ch;
winArt = p.Results.winArt;
winDC = p.Results.winDC;

try [sessionInfo] = bz_getSessionInfo(basepath, 'noPrompts', true);
    fs = sessionInfo.rates.wideband;
    nChannels = sessionInfo.nChannels;
catch
    warning('SessionInfo file not found.');
end

if ischar('ch') && strcmpi(ch, 'all')
    ch = 1:nChannels;
else
    ch = ch + 1;
end

% ts = int32(ts * fs);
% winArt = winArt * fs;
% winDC = winDC * fs;


if ~exist([basepath filesep strcat('copy_bin.dat')],'file')
    disp('Creating .dat back up...');
    copyfile(strcat(filename,'.dat'),'copy_bin.dat'); % create .dat back up
end


%% Make LFP
[~,basename] = fileparts(basepath);

if isempty(dir('*.lfp'))
    try 
        bz_LFPfromDat(pwd,'outFs',1250); % generating lfp, NOTE: THIS FUNCTION WILL GENERATE A SESSIONINFO FILE!! WE NEED TO FIX THIS
        disp('Making LFP file ...');
    catch
        disp('Problems with bz_LFPfromDat, resampling...');
%         sessionFile = dir('*session.mat*');
%         if ~isempty(sessionFile)
%             load(sessionFile.name);
%         end
%         ResampleBinary(strcat(basename,'.dat'),strcat(basename,'.lfp'),...
%             session.extracellular.nChannels,1, session.extracellular.sr/session.extracellular.srLfp);
ResampleBinary(strcat(basename,'.dat'),strcat(basename,'.lfp'),...
            sessionInfo.nChannels,1, sessionInfo.rates.wideband/sessionInfo.rates.lfp);
        
    end
end




m = memmapfile(fullfile(basepath,strcat(filename,'.dat')),'Format','int16','Writable', true);
data=reshape(m.Data,nChannels,[]);

tamano = get(0,'ScreenSize');
% figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
% plot(data(1,:));
% for i=1:nChannels
%     subplot(nChannels,1,i)
%     plot(data(i,:));
% end
lfp = bz_GetLFP([1]);
figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
plot(lfp.data);

% [specSlope,spec] = bz_PowerSpectrumSlope(lfp,2,1,'showfig',true);

% Fourier analysis
% T = 1/fs;
% L = length(lfp.data);
% t = (0:L-1)*T;
% Y = fft(lfp.data);
% figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
% plot(abs(Y));

params= struct('Fs',1250,'pad',1,'tapers',[5 3],'fpass',[1 1000],'err',[2 0.05]);
% [S,f,Serr] = mtspectrumc(double(lfp.data),params);
%      
% figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
% plot(f,10*log10(S),'Color', 'b','LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)')



% [S,t,f,Serr] = mtspecgramc(double(lfp.data),[1 0.01], params);
% figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
% plot(f,10*log10(mean(S)),'Color', 'k','LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)')
% figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
% figure,imagesc(t,f,10*log10(S*10e12)'),view(0,-90),colorbar, colormap(jet)


%% Vamos a intentar eliminar los artefactos

Plott = 1;
[CleanDa,timestamps] = ArtifactDetection(double(lfp.data),1250,Plott);

figure('position',[tamano(1) tamano(2) tamano(3) tamano(4)]);
subplot(2,1,1)
plot(lfp.data), hold; text( timestamps,double(lfp.data(timestamps)),'0')

% t = 0:1/1250:size(data,2);
subplot(2,1,2)
plot(CleanDa)
% Parece que con la primera aproximación (detrending) no surge efecto en la
% señal



%% Vamos a probar con la función de Manu cleanPulses ya que ahora tenemos los picos que queremos eliminar
timestamps = sort(timestamps);
ts = timestamps/1250;
[data_clean] = cleanPulses(ts,'winArt',0.5,'ch',0);

%%
filename = [basename,'.dat'];
nChansTotal = 64;

medianTrace = applyCARtoDat(filename,nChansTotal);




