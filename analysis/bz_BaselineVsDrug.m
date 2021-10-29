function bz_BaselineVsDrug(varargin)

% Function to plot in the same figure the baseline vs drug condition

% USAGE 
%   bz_BaselineVsDrug
%   
% INPUTS
%
%
%
%
% OUTPUTS
%
%
%
%
%
%
% Created by Pablo Abad 2021

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'drug',[],@isstr);
addParameter(p,'openField',true,@islogical);
addParameter(p,'yMaze',true,@islogical);
addParameter(p,'saveFig',true,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
drug = p.Results.drug;
openField = p.Results.openField;
yMaze = p.Results.yMaze;
saveFig = p.Results.saveFig;

mkdir('BaselineVSDrug')

%% Load sessionInfo, session, MergePoints and Tracking

sessionInfo = bz_getSessionInfo(basepath);
session = loadSession();

if ~isempty('*.MergePoints.events.mat')
    file = dir('*MergePoints.events.mat');
    load(file.name)
end

if ~isempty('*Coherogram.SubSession.lfp.mat')
    file = dir('*Coherogram.SubSession.lfp.mat');
    load(file.name)
end


count_openField = 0;
count_YMaze = 0;

for i=1:length(MergePoints.foldernames)
    
    if ~isempty('*Tracking.Behavior.mat')
            file = dir('*Tracking.Behavior.mat');
            load(file.name)
    end
        
    if ismember(MergePoints.foldernames{i},tracking.folders)
        cd(MergePoints.foldernames{i})
        if ~isempty('*Tracking.Behavior.mat')
            file = dir('*Tracking.Behavior.mat');
            load(file.name)
        end
        
        if openField
            if strcmpi(tracking.apparatus.name, 'Open Field')
                count_openField = count_openField + 1;
                data_openField{count_openField} = coherogram.(MergePoints.foldernames{i});
            end
        end
        cd(basepath)
        
        if yMaze
            if strcmpi(tracking.apparatus.name,'YMaze Apparatus')
                count_YMaze = count_YMaze +1;
                data_YMaze{count_YMaze} = coherogram.(MergePoints.foldernames{i});
            end
        end
        cd(basepath)
        
    end
end

%% Plot Figures Baseline vs Drug Condition
col_baseline = [0 0 0];
switch drug
    case 'vehicle'
        col_def = [0.5 0.5 0.5];
        legend_name = 'vehicle';
    case 'ketamine'
        col_def = [0 0.7 1];
        legend_name = 'ketamine';
    case 'mk801'
        col_def = [1 0 0];
        legend_name = 'MK801';
end

if openField
    figure,
    set(gcf,'Position',get(0,'ScreenSize'))
    
    for i = 1:length(data_openField)
        if i == 1
            subplot(1,4,1), plot(data_openField{i}.f,mean(data_openField{i}.coherogram,2),'Color',col_baseline,'LineWidth',1),ylabel('coherence (r)'),xlabel ('Frequency (f)'), title(['Coherence Ch: ', num2str(data_openField{i}.ch1), 'and Ch: ', num2str(data_openField{i}.ch2),' filtered']),axis([0 200 0 1]); hold on;
            subplot(1,4,2),plot(data_openField{i}.f,mean(data_openField{i}.phase,2),'Color',col_baseline,'LineWidth',1), ylabel(' Coherence in frequency(r)'), xlabel('Frequency(f)'), title(['Coherence Ch: ', num2str(data_openField{i}.ch1), 'and Ch: ', num2str(data_openField{i}.ch2), 'filtered']), axis([0 200 -pi pi]), hold on;
            subplot(1,4,3),plot(data_openField{i}.f,10*log10(mean(data_openField{i}.S1)),'Color',col_baseline,'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(data_openField{i}.ch1), 'filtered']), hold on;
            subplot(1,4,4),plot(data_openField{i}.f,10*log10(mean(data_openField{i}.S2)),'Color',col_baseline,'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(data_openField{i}.ch1), 'filtered']), hold on;
        elseif i == 2
            subplot(1,4,1), plot(data_openField{i}.f,mean(data_openField{i}.coherogram,2),'Color',col_def,'LineWidth',1),ylabel('coherence (r)'),xlabel ('Frequency (f)'), title(['Coherence Ch: ', num2str(data_openField{i}.ch1), 'and Ch: ', num2str(data_openField{i}.ch2),' filtered']),axis([0 200 0 1]); hold on;
            subplot(1,4,2),plot(data_openField{i}.f,mean(data_openField{i}.phase,2),'Color',col_def,'LineWidth',1), ylabel(' Coherence in frequency(r)'), xlabel('Frequency(f)'), title(['Coherence Ch: ', num2str(data_openField{i}.ch1), 'and Ch: ', num2str(data_openField{i}.ch2), 'filtered']), axis([0 200 -pi pi]), hold on;
            subplot(1,4,3),plot(data_openField{i}.f,10*log10(mean(data_openField{i}.S1)),'Color',col_def,'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(data_openField{i}.ch1), 'filtered']), hold on;
            subplot(1,4,4),plot(data_openField{i}.f,10*log10(mean(data_openField{i}.S2)),'Color',col_def,'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(data_openField{i}.ch1), 'filtered']), hold on;
        end
    end
    legend('baseline',legend_name)
    
    if saveFig
        saveas(gcf,['BaselineVSDrug\OpenField.png']);
    end
end


if yMaze
    figure,
    set(gcf,'Position',get(0,'ScreenSize'))
    
    for i = 1:length(data_openField)
        if i == 1
            subplot(1,4,1), plot(data_YMaze{i}.f,mean(data_YMaze{i}.coherogram,2),'Color',col_baseline,'LineWidth',1),ylabel('coherence (r)'),xlabel ('Frequency (f)'), title(['Coherence Ch: ', num2str(data_YMaze{i}.ch1), 'and Ch: ', num2str(data_YMaze{i}.ch2),' filtered']),axis([0 200 0 1]); hold on;
            subplot(1,4,2),plot(data_YMaze{i}.f,mean(data_YMaze{i}.phase,2),'Color',col_baseline,'LineWidth',1), ylabel(' Coherence in frequency(r)'), xlabel('Frequency(f)'), title(['Coherence Ch: ', num2str(data_YMaze{i}.ch1), 'and Ch: ', num2str(data_YMaze{i}.ch2), 'filtered']), axis([0 200 -pi pi]), hold on;
            subplot(1,4,3),plot(data_YMaze{i}.f,10*log10(mean(data_YMaze{i}.S1)),'Color',col_baseline,'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(data_YMaze{i}.ch1), 'filtered']), hold on;
            subplot(1,4,4),plot(data_YMaze{i}.f,10*log10(mean(data_YMaze{i}.S2)),'Color',col_baseline,'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(data_YMaze{i}.ch1), 'filtered']), hold on;
        elseif i == 2
            subplot(1,4,1), plot(data_YMaze{i}.f,mean(data_YMaze{i}.coherogram,2),'Color',col_def,'LineWidth',1),ylabel('coherence (r)'),xlabel ('Frequency (f)'), title(['Coherence Ch: ', num2str(data_YMaze{i}.ch1), 'and Ch: ', num2str(data_YMaze{i}.ch2),' filtered']),axis([0 200 0 1]); hold on;
            subplot(1,4,2),plot(data_YMaze{i}.f,mean(data_YMaze{i}.phase,2),'Color',col_def,'LineWidth',1), ylabel(' Coherence in frequency(r)'), xlabel('Frequency(f)'), title(['Coherence Ch: ', num2str(data_YMaze{i}.ch1), 'and Ch: ', num2str(data_YMaze{i}.ch2), 'filtered']), axis([0 200 -pi pi]), hold on;
            subplot(1,4,3),plot(data_YMaze{i}.f,10*log10(mean(data_YMaze{i}.S1)),'Color',col_def,'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(data_YMaze{i}.ch1), 'filtered']), hold on;
            subplot(1,4,4),plot(data_YMaze{i}.f,10*log10(mean(data_YMaze{i}.S2)),'Color',col_def,'LineWidth',1),ylabel('10*10log(X)'), xlabel('Frequency(f)'), title(['Power Ch: ', num2str(data_YMaze{i}.ch1), 'filtered']), hold on;
        end
    end
    legend('baseline',legend_name)
    
    if saveFig
        saveas(gcf,['BaselineVSDrug\YMaze.png']);
    end
end




end

