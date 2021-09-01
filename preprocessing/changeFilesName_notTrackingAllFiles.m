function [outputArg1,outputArg2] = changeFilesName_notTrackingAllFiles(basepath)

% USAGE
%   changeFilesName(varargin)
% 
% INPUT (OPTIONAL)
% basepath      By default pwd

%% Defaults and Parms
if nargin < 1
    basepath = pwd;
end
rename = 1;
% prevPath = pwd;

% Reading name of the files and being sure that not files are mixed and
% there is the same number of the different type

cd(basepath);
animalName = strsplit(basepath,filesep);
animalName = animalName{end};

% Find all the _data.xdat in the folder
xdat_files = dir('*_data.xdat*');
timestamp_xdat_file = dir('*_timestamp.xdat*');
json_file = dir('*.xdat.json*');
tracking_xml_file = dir('*.xml*');
tracking_csv_file = dir('*.csv*');
txt_file = dir('*.txt*'); % File only present in TMaze recordings

% xml_file = dir('*global.xml*');

% We need to check that global.xml file is not included in
% tracking_xml_file
for i=1:length(tracking_xml_file)    
    globalXMLPosition(i) = strcmp(tracking_xml_file(i).name,'global.xml');
    if exist('globalXMLPosition','var')
        if globalXMLPosition(i) == 1
            disp('Excluding global.xml file from the .xml tracking files');
            globalXMLPositionToDelete = i;
        end
    end
end


% tracking_xml_file(globalXMLPositionToDelete) = [];


if isequal(length(json_file),length(timestamp_xdat_file),length(xdat_files))
    disp('Correct number of ephys files. Changing the name ...')
else
    disp('There are ephys files missing. Quiting...');
end

if isequal(length(tracking_csv_file),length(tracking_xml_file))
    disp('Correct numver of tracking files')
else
    disp('There are tracking files missing. Quiting...')
    return
end

if ~isequal(length(xdat_files),length(tracking_csv_file))
    disp('Not all recordings have tracking')
    trackingANDephys = 0;
end
%% Changing name of the ephys files
if ~isempty(xdat_files)
% Change name of all the files
        for i=1:size(xdat_files,1)

            file{i} = xdat_files(i).name;
            file_date{i} = xdat_files(i).date;
            
            file_timestamp{i} = timestamp_xdat_file(i).name;
            file_timestamp_date{i} = timestamp_xdat_file(i).date;
            
            file_json{i} = json_file(i).name;
            file_json_date{i} = json_file(i).date;
            
            if rename == 1
                a = find(file{i} == '-');
                newFileName = file{i};
                spaceToDelete = [a(1)-10:a,a(2),a(3)]
                newFileName(spaceToDelete) = [];

                [filepath,filename_,ext] = fileparts(newFileName);
                filename_(end-4:end) = [];
                filename{i} = filename_;

            end

        end
        
        for i=1:size(tracking_xml_file,1)
            file_tracking_xml{i} = tracking_xml_file(i).name;
            file_tracking_xml_date{i} = tracking_xml_file(i).date;
            
            file_tracking_csv{i} = tracking_csv_file(i).name;
            file_tracking_csv_date{i} = tracking_csv_file(i).date;
        end
        
        for i=1:size(txt_file,1)
            if size(txt_file,1) == 1
                file_txt{i} = txt_file.name;
            else
                file_txt{i} = txt_file{i}.name;
            end
        end
     
        % Now we find which ephyhs recording corresponds to the tracking
        % recording
        
        
        
            % Now we need to change the name of the files and the folder
        % accordingly
        if rename == 1
            for i=1:size(file,2)
                movefile(file{i},strcat(filename{i},'_data.xdat'));
                movefile(file_timestamp{i},strcat(filename{i},'_timestamp.xdat'));
                movefile(file_json{i},strcat(filename{i},'.xdat.json'));
            end
        end
        
        for i=1:size(filename,2)
            if i <= 5
                if rem(i,2) == 0
                    numberOddRecordings(i) = i;
                end
            elseif i > 5
                if rem(i,2) ~= 0
                    numberOddRecordings(i) = i;
                end
            end
        end
        
        if exist('numberOddRecordings','var')
            numberOddRecordings(numberOddRecordings == 0) = [];
        end
        
        
        if rename == 1
            for i=1:size(file_tracking_csv,2)
                movefile(file_tracking_csv{i},strcat(filename{numberOddRecordings(i)},'_tracking.csv'));
                movefile(file_tracking_xml{i},strcat(filename{numberOddRecordings(i)},'_tracking.xml'));
            end
        end
        
        if exist('file_txt','var')
            if rename==1
                for i=1:size(file_txt)
                    movefile(file_txt{i},strcat(filename{numberOddRecordings},'_TMaze.txt'));
                end
            end
        end

        % Folder
        % Now we need to copy each file pertaining to the corresponding
        % folder
        %%        
       
            for i=1:size(file,2)
               basepath_prev = strsplit(basepath,animalName);
               cd(basepath_prev{1})
               mkdir(strcat(basepath,'_',filename{i}(end-5:end)))
               cd(basepath)

                movefile(strcat(filename{i},'_data.xdat'),strcat(basepath,'_',filename{i}(end-5:end)));
%                 cd(basepath)
                movefile(strcat(filename{i},'.xdat.json'),strcat(basepath,'_',filename{i}(end-5:end)));
%                 cd(basepath)
                movefile(strcat(filename{i},'_timestamp.xdat'),strcat(basepath,'_',filename{i}(end-5:end)));
%                 cd(basepath)
                if i <= 5
                    if rem(i,2) == 0
                        movefile(strcat(filename{i},'_tracking.csv'),strcat(basepath,'_',filename{i}(end-5:end)));
        %                 cd(basepath)
                        movefile(strcat(filename{i},'_tracking.xml'),strcat(basepath,'_',filename{i}(end-5:end)));
        %                 cd(basepath)
                        if exist('file_txt','var')

                            movefile(strcat(filename{i},'_TMaze.txt'),strcat(basepath,'_',filename{i}(end-5:end)));
                        end
                    end
                elseif i > 5
                    if rem(i,2) ~= 0
                        movefile(strcat(filename{i},'_tracking.csv'),strcat(basepath,'_',filename{i}(end-5:end)));
        %                 cd(basepath)
                        movefile(strcat(filename{i},'_tracking.xml'),strcat(basepath,'_',filename{i}(end-5:end)));
        %                 cd(basepath)
                        if exist('file_txt','var')

                            movefile(strcat(filename{i},'_TMaze.txt'),strcat(basepath,'_',filename{i}(end-5:end)));
                        end
                        
                    end        
                end
                
            end
       %%
        
        
        basepath_prev = strsplit(basepath,animalName);
        cd(basepath_prev{1})
        
%         rmdir HPS22_200421 s
        basedatapath = strsplit(basepath_prev{1},animalName(1:5))
        movefile(basepath,basedatapath{1})
        
end



end