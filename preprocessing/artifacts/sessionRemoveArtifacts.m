function [outputArg1,outputArg2] = sessionRemoveArtifacts(varargin)

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
filename = split(pwd,'\');
filename = filename(end);

% Pase options
p = inputParser;
addParameter(p,'filename',filename,@isstr);
addParameter(p,'basepath',pwd,@isstr);

parse(p,varargin{:});

filename = p.Results.filename;
basepath = p.Results.basepath;


try [sessionInfo] = bz_getSessionInfo(basepath,'noPrompts',true);
    fs = sessionInfo.rates.wideband;
    nChannels = sessionInfo.nChannels;
catch
    warning('SessionInfo file not found')
end

% Move to the folder where are located TMaze and linearMaze files

if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
        count = 1;
        for ii = 1:size(MergePoints.foldernames,2)
            %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
             if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*.csv*']))   
                cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
                fprintf('Removing electrical artifacts in %s folder \n',MergePoints.foldernames{ii});
                 tempRemovingArtifacts{count} = removeArtifacts();
                 artifactsFolder(count) = ii;
                count = count + 1;
            end
        end
        cd(basepath);
    else
        error('missing MergePoints, quiting...');
     end





end

