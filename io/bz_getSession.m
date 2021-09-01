function [session] = bz_getSession(basePath,varargin)
%[session] = bz_getSession(basePath) loads the session metadata
%for the recording in basePath. basePath should be in the format:
%       /whateverPath/baseName/
%           a file  basePath/baseName.session.mat
%           or      basePath/baseName.xml
%           should exist.
%If no baseName.session.mat exists, loads from the xml.
%
%INPUT
%   basePath            directory: '/whatevetPath/baseName/'
%   (options)
%       'saveMat'       (default: prompt)
%       'noPrompts'     (default: false) prevents prompts about
%                       saving/adding metadata
%       'editGUI'       (default: false) opens a GUI to edit select
%                       sessionInfo fields (beta, please add/improve!)
%
%OUTPUT
%   session         metadata structure
%
%2017 DLevenstein and DTingley
%% inputs and defaults
p = inputParser;
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'editGUI',false,@islogical);
addParameter(p,'saveMat',false,@islogical);
parse(p,varargin{:})
noPrompts = p.Results.noPrompts;
editGUI = p.Results.editGUI;
saveMat = p.Results.saveMat;

if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);
filename = fullfile(basePath,[baseName,'.session.mat']);

baseName_prev = strsplit(fullfile(pwd),['\' baseName]);
baseName_prev = baseName_prev{1};
sessionName_prev = strsplit(baseName_prev,filesep);
sessionName_prev = sessionName_prev{end-1};
filename_prev = fullfile(baseName_prev,[sessionName_prev,'.session.mat']);
%% Load the stuff
%d = dir('*sessionInfo*'); %all files with sessioninfo in the name
if exist(filename,'file')
    sessionInfostruct = load(filename);
    %Checks that there is a single structure in the sessionInfo file
    varsInFile = fieldnames(sessionInfostruct); 
    if numel(varsInFile)==1
        sessionInfo = sessionInfostruct.(varsInFile{1});
    else
        warning('Your .session.mat has multiple variables/structures in it... wtf.')
        sessionInfo = sessionInfostruct;
    end
    SIexist = true;  %Marks that session info exists as expected
else
    cd(baseName_prev)
    if exist(filename_prev,'file')
        sessionInfostruct = load(filename_prev);
        %Checks that there is a single structure in the sessionInfo file
        varsInFile = fieldnames(sessionInfostruct); 
        if numel(varsInFile)==1
            sessionInfo = sessionInfostruct.(varsInFile{1});
        else
            warning('Your .sessionInfo.mat has multiple variables/structures in it... wtf.')
            sessionInfo = sessionInfostruct;
        end
        SIexist = true;  %Marks that session info exists as expected
        cd(basePath)
    else
        cd(basePath)
       warning(['could not find file ',baseName,'.sessionInfo.mat ',...
           'running LoadParameters instead..']) 
       sessionInfo = LoadParameters(basePath);
       SIexist = false; 
    end
end

%% Check sessionInfo using bz_isSessionInfo and update if necessary
bz_isSessionInfo(sessionInfo);

%Here: prompt user to add any missing sessionInfo fields and save
if editGUI
    [ sessionInfo ] = bz_sessionInfoGUI(sessionInfo);
    SIexist = false;
% elseif ~isfield(sessionInfo,'region') && ~noPrompts
%     regionadd = questdlg(['Your sessionInfo is missing regions, ',...
%         'would you like to add them?'],'Add Regions?','Yes');
%     switch regionadd
%         case 'Yes'
%             [sessionInfo] = bz_sessionInfoGUI(sessionInfo,'Regions');
%             SIexist = false; 
%         case 'Cancel'
%             return
%     end
end

%Should check that sessionInfo.session.name and sesioninfo.session.path
%match basePath....  if not, prompt the user to save with the correct
%information.
    
%% Save sessionInfo file   
%Prompt user to save basePath/baseName.sessionInfo.mat 
%if loaded from xml or changed
if ~noPrompts && ~SIexist %Inform the user that they should save a file for later
    savebutton = questdlg(['Would you like to save your sessionInfo in ',...
        filename '?'],'Save sessionInfo?','Yes');
    switch savebutton
        case 'Yes'
            saveMat = true;
        case 'No'
            saveMat = false;
        case 'Cancel'
            return
    end
end

if saveMat
    disp(['saving ',baseName,'.sessionInfo.mat'])
    save(filename,'sessionInfo'); 
end
end


