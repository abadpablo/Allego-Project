function tracking = anyMazeTracking(csvfile,xmlfile,varargin)

% Need to include the zones for the tracking

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'fs',30,@isnumeric);
addParameter(p,'artifactThreshold',10,@isnumeric);
addParameter(p,'convFact',[],@isnumeric); % 0.1149
addParameter(p,'roiTracking',[],@ismatrix);
addParameter(p,'roiLED',[],@ismatrix);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'saveFrames',true,@islogical)
addParameter(p,'verbose',false,@islogical);
addParameter(p,'thresh',.98,@isnumeric)
addParameter(p,'bazlerTTL',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'orderKalmanVel',2,@isnumeric);

% addParameter(p,'RGBChannel',[],@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
fs = p.Results.fs;
artifactThreshold = p.Results.artifactThreshold;
convFact = p.Results.convFact;
roiTracking = p.Results.roiTracking;
roiLED = p.Results.roiLED;
forceReload = p.Results.forceReload;
saveFrames = p.Results.saveFrames;
verbose = p.Results.verbose;
thresh = p.Results.thresh;
bazlerTtl = p.Results.bazlerTTL;
saveMat = p.Results.saveMat;
csvfile = [];
xmlfile = [];
order = p.Results.orderKalmanVel;
%% Deal with inputs
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) || forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end


if ~exist('csvfile') || isempty(csvfile)
    if ~isempty(dir([basepath filesep '*csv']))
        csvfile = dir([basepath filesep '*csv*']); 
        csvfile = erase(csvfile.name,'.csv');
    else
        warning('No csv file ( needed for position variables)!!');
        tracking = [];
        return
    end
end

% Load the csv file
if ~isempty(csvfile)
    CSV = readtable(strcat(csvfile,'.csv'));
else
    disp('Tracking csv file not found')
    return
end

    variableNames = CSV.Properties.VariableNames;
    
    time_column = strcmp(variableNames,'Time');
    centreX_column = strcmp(variableNames,'CentrePositionX');
    centreY_column = strcmp(variableNames,'CentrePositionY');
    speed_column = strcmp(variableNames,'Speed');
    headX_column = strcmp(variableNames,'HeadPositionX');
    headY_column = strcmp(variableNames,'HeadPositionY');
    tailX_column = strcmp(variableNames,'TailPositionX');
    tailY_column = strcmp(variableNames,'TailPositionY');
    % YMaze
    inArm1_column = strcmp(variableNames,'InArm1');
    inArm2_column = strcmp(variableNames,'InArm2');
    inArm3_column = strcmp(variableNames,'InArm3');
    inCenter_column = strcmp(variableNames,'InCenter');
    inNoZone_column = strcmp(variableNames,'InNoZone');
    YMazeSequenceStarts_column = strcmp(variableNames,'YMazeSequenceStarts');
    YMazeSequenceEnds_column = strcmp(variableNames,'YMazeSequenceEnds');
    % Open Field
    inMain_column = strcmp(variableNames,'InMain');
    % TMaze
    inLeftReward_column = strcmp(variableNames,'InLeftReward');
    inRightReward_column = strcmp(variableNames,'InRightReward');
    inDecision_column = strcmp(variableNames,'InDecision');
    inStarting_column = strcmp(variableNames,'InStarting');
    
    
    prueba_column = strcmp(variableNames,'prueba');
    
    
    % We need to be sure that the column of time is not a duration variable
      
    % If the time variable is already in seconds there are no problems in
    % the table2array conversion, but if it is in HH:MM:SS we need to
    % convert the data
    
    CSVArray = table2array(CSV);
    
    if isduration(CSVArray)
        num_variables = length(variableNames);
        for i=1:num_variables
            varTypes{i} = 'double';
        end       
        CSV_table = table('Size',size(CSV),'VariableTypes',varTypes,'VariableNames',variableNames);
        for i=1:num_variables
            if i==1
            CSV_table.(variableNames{i}) = seconds(CSV.(variableNames{i}));
            else
                CSV_table.(variableNames{i}) = CSV.(variableNames{i});
            end
        end
        clear CSVArray
        CSVArray = table2array(CSV_table);
    end
    
% headX_column = strcmp(variableNames,'HeadPositionX');


timesamples = CSVArray(:,time_column);
xPos = CSVArray(:,centreX_column);
yPos = CSVArray(:,centreY_column);
speed = CSVArray(:,speed_column);

head_xPos = CSVArray(:,headX_column);
head_yPos = CSVArray(:,headY_column);
tail_xPos = CSVArray(:,tailX_column);
tail_yPos = CSVArray(:,tailY_column);
% YMaze
inArm1 = CSVArray(:,inArm1_column);
inArm2 = CSVArray(:,inArm2_column);
inArm3 = CSVArray(:,inArm3_column);
inCenter = CSVArray(:,inCenter_column);
inNoZone = CSVArray(:,inNoZone_column);
% Open Field
inMain = CSVArray(:,inMain_column);
% TMaze
inLeftReward = CSVArray(:,inLeftReward_column);
inRightReward = CSVArray(:,inRightReward_column);
inDecision = CSVArray(:,inDecision_column);
inStarting = CSVArray(:,inStarting_column);

%% Load the tracking_xmlFile

if ~exist('xmlfile') || isempty(xmlfile)
    if ~isempty(dir([basepath filesep '*xml']))
        xmlfile = dir([basepath filesep '*xml*']); 
        xmlfile = erase(xmlfile.name,'.xml');
    else
        warning('No xml file (needed for apparatus bounding)!!');
        tracking = [];
        return
    end
end

if ~isempty(xmlfile)
    XML = xml2struct(strcat(xmlfile,'.xml'));
else
    disp('Tracking xml file not found')
    return
end

% Need to extract the values of the xml file
apparatus_name = [];
[m,n] = size(XML.experiment.animal.test.apparatus);
if n == 1
    pixels_metre = str2num(XML.experiment.animal.test.pixelspermetre.Text);
    apparatus = XML.experiment.animal.test.apparatus;
    if isfield(XML.experiment.animal.test,'zone')
        zone = XML.experiment.animal.test.zone;
    end
else
    pixels_metre = str2num(XML.experiment.animal.test.pixelspermetre.Text);
    apparatus = XML.experiment.animal.test.apparatus;
    if isfield(XML.experiment.animal.test,'zone')
        zone = XML.experiment.animal.test.zone;
    end
end

[m_apparatus,n_apparatus] = size(apparatus);

if n_apparatus == 2
    apparatus_name = apparatus{1}.Text;
    centre = apparatus{2}.centre;
    boundingbox = apparatus{2}.boundingbox;
else
    if isfield(apparatus,'Text')
        apparatus_name = apparatus.Text;
    end
    
    if isfield(apparatus,'centre')
        centre = apparatus.centre;
    end
    
    if isfield(apparatus,'boundingbox')
        boundingbox = apparatus.boundingbox;
    end
    disp('No apparatus name saved')
end

%% APPARATUS

% centre = apparatus{2}.centre;
% centre = apparatus.centre;
centre_x = str2num(centre.x.Text);
centre_y = str2num(centre.y.Text);

% boundingbox = apparatus{2}.boundingbox;
% boundingbox = apparatus.boundingbox;
boundingbox_x = str2num(boundingbox.x.Text);
boundingbox_y = str2num(boundingbox.y.Text);
boundingbox_w = str2num(boundingbox.w.Text);
boundingbox_h = str2num(boundingbox.h.Text);

%% CONVERT BOUNDING BOX AND CENTRE TO CM
centre_x = centre_x*100/pixels_metre;
centre_y = centre_y*100/pixels_metre;
boundingbox_x = boundingbox_x*100/pixels_metre;
boundingbox_y = boundingbox_y*100/pixels_metre;
boundingbox_w = boundingbox_w*100/pixels_metre;
boundingbox_h = boundingbox_h*100/pixels_metre;

%%
boundingbox_X = boundingbox_x + boundingbox_w;
boundingbox_Y = boundingbox_y + boundingbox_h;

boundingbox_xmin = boundingbox_x;
boundingbox_xmax = boundingbox_X;
boundingbox_ymin = boundingbox_y;
boundingbox_ymax = boundingbox_Y;
boundingbox_xy = [boundingbox_xmax-boundingbox_xmin;boundingbox_ymax-boundingbox_ymin];

xMaze = [boundingbox_xmin boundingbox_xmax];
yMaze = [boundingbox_ymin boundingbox_ymax];


%% ZONES

if exist('zone','var')
    num_zones = length(zone);
    for i=1:num_zones
        if num_zones == 1
            name_zones{i} = zone.name.Text;
            center_zones_x{i} = str2num(zone.centre.x.Text)*100/pixels_metre;
            center_zones_y{i} = str2num(zone.centre.y.Text)*100/pixels_metre;
            boundingbox_zones_x{i} = str2num(zone.boundingbox.x.Text)*100/pixels_metre;
            boundingbox_zones_y{i} = str2num(zone.boundingbox.y.Text)*100/pixels_metre;
            boundingbox_zones_w{i} = str2num(zone.boundingbox.w.Text)*100/pixels_metre;
            boundingbox_zones_h{i} = str2num(zone.boundingbox.h.Text)*100/pixels_metre;
            
            

            boundingbox_zones_X{i} = boundingbox_zones_x{i} + boundingbox_zones_w{i};
            boundingbox_zones_Y{i} = boundingbox_zones_y{i} + boundingbox_zones_h{i};

            boundingbox_zones_xmin{i} = boundingbox_zones_x{i};
            boundingbox_zones_xmax{i} = boundingbox_zones_X{i};
            boundingbox_zones_ymin{i} = boundingbox_zones_y{i};
            boundingbox_zones_ymax{i} = boundingbox_zones_Y{i};

            xMaze_zones{i} = [boundingbox_zones_xmin{i} boundingbox_zones_xmax{i}];
            yMaze_zones{i} = [boundingbox_zones_ymin{i} boundingbox_zones_ymax{i}];
        else
            name_zones{i} = zone{i}.name.Text;
            center_zones_x{i} = str2num(zone{i}.centre.x.Text)*100/pixels_metre;
            center_zones_y{i} = str2num(zone{i}.centre.y.Text)*100/pixels_metre;
            boundingbox_zones_x{i} = str2num(zone{i}.boundingbox.x.Text)*100/pixels_metre;
            boundingbox_zones_y{i} = str2num(zone{i}.boundingbox.y.Text)*100/pixels_metre;
            boundingbox_zones_w{i} = str2num(zone{i}.boundingbox.w.Text)*100/pixels_metre;
            boundingbox_zones_h{i} = str2num(zone{i}.boundingbox.h.Text)*100/pixels_metre;

            boundingbox_zones_X{i} = boundingbox_zones_x{i} + boundingbox_zones_w{i};
            boundingbox_zones_Y{i} = boundingbox_zones_y{i} + boundingbox_zones_h{i};

            boundingbox_zones_xmin{i} = boundingbox_zones_x{i};
            boundingbox_zones_xmax{i} = boundingbox_zones_X{i};
            boundingbox_zones_ymin{i} = boundingbox_zones_y{i};
            boundingbox_zones_ymax{i} = boundingbox_zones_Y{i};

            xMaze_zones{i} = [boundingbox_zones_xmin{i} boundingbox_zones_xmax{i}];
            yMaze_zones{i} = [boundingbox_zones_ymin{i} boundingbox_zones_ymax{i}];
        end

    end    
else
    name_zones = [];
    center_zones_x{i} = [];
    center_zones_y{i} = [];
    boundingbox_zones_x{i} = [];
    boundingbox_zones_y{i} = [];
    boundingbox_zones_w{i} = [];
    boundingbox_zones.h{i} = [];
    boundingbox_zones_X{i} = [];
    boundingbox_zones_Y{i} = [];
        
    boundingbox_zones_xmin{i} = [];
    boundingbox_zones_xmax{i} = [];
    boundingbox_zones_ymin{i} = [];
    boundingbox_zones_ymax{i} = [];
        
    xMaze_zones{i} = [];
    yMaze_zones{i} = [];
        
end


%% CONVERTING POSITION VARIABLES INTO CM
xPos = xPos*100/pixels_metre;
yPos = yPos*100/pixels_metre;

head_xPos = head_xPos*100/pixels_metre;
head_yPos = head_yPos*100/pixels_metre;

tail_xPos = tail_xPos*100/pixels_metre;
tail_yPos = tail_yPos*100/pixels_metre;

%% Figure no filtered tracking data
% Body tracking
figure,
plot(xPos,yPos,'k'), title('No filtered body tracking data')
hold on;
plot(centre_x,centre_y,'b+'),
xlim([min(xPos) max(xPos)])
ylim([min(yPos) max(yPos)])


% Head tracking
if ~isempty('head_xPos')
figure,
plot(head_xPos,head_yPos,'k'), title('No filtered head tracking data')
hold on;
plot(centre_x,centre_y,'b+'),
xlim([min(head_xPos) max(head_xPos)])
ylim([min(head_yPos) max(head_yPos)])
end

% Tail tracking
if ~isempty('tail_xPos')
figure,
plot(tail_xPos,tail_yPos,'k'), title('No filtered tail tracking data')
hold on;
plot(centre_x,centre_y,'b+'),
xlim([min(tail_xPos) max(tail_xPos)])
ylim([min(tail_yPos) max(tail_yPos)])
end  





%% Filtering tracking data
pos = [xPos'; yPos']';
if ~isempty(head_xPos)
    head_pos = [head_xPos'; head_yPos']';
else
    head_pos = [];
end
if ~isempty(tail_xPos)
    tail_pos = [tail_xPos';head_yPos']';
else
    tail_pos = [];
end

% xt = linspace(0,size(pos,1)/fs,size(pos,1));
% Post Processing of tacking

[t,x,y,~,~,~,~] = trajectory_kalman_filter(pos(:,1)',pos(:,2)',timesamples,0);

[~,~,~, vx, vy, ax, ay] = KalmanVel(pos(:,1)',pos(:,2)',timesamples,order);
v = sqrt(vx.^2+vy.^2);
ace = sqrt(ax.^2+ay.^2);
if length(v) < length(pos(:,1))
    dif = length(pos(:,1)) - length(v);
    s = v;
    sx = vx;
    sy = vy;
    a = ace;
    acex = ax;
    acey = ay;
    clear v
    clear vx
    clear vy
    clear ace
    clear ax
    clear ay
    
    v = zeros(1,length(pos(:,1)));
    vx = zeros(1,length(pos(:,1)));
    vy = zeros(1,length(pos(:,1)));
    ace = zeros(1,length(pos(:,1)));
    ax = zeros(1,length(pos(:,1)));
    ay = zeros(1,length(pos(:,1)));

       
    v = [nan(dif,1);s];
    vx = [nan(dif,1); sx];
    vy = [nan(dif,1); sy];
    ace = [nan(dif,1); a];
    ax = [nan(dif,1); acex];
    ay = [nan(dif,1); acey];
    mspeed = nanmean(v);
    mace = nanmean(ace);

end

if length(x) < length(pos(:,1))
    dif = length(pos(:,1)) - length(x);
    pos_x = x;
    pos_y = y;
    clear x
    clear y
    x = zeros(1,length(pos(:,1)));
    y = zeros(1,length(pos(:,2)));    
    x = [nan(dif,1);pos_x];
    y = [nan(dif,1);pos_y];
end

if ~isempty(head_pos)
    [t,head_x,head_y,~,~,~,~] = trajectory_kalman_filter(head_pos(:,1)',head_pos(:,2)',timesamples,0);
    
    [~,~,~, head_vx, head_vy, head_ax, head_ay] = KalmanVel(head_pos(:,1)',head_pos(:,2)',timesamples,order);
    head_v = sqrt(head_vx.^2+head_vy.^2);
    head_ace = sqrt(head_ax.^2+head_ay.^2);
    if length(head_v) < length(head_pos(:,1))
        dif = length(head_pos(:,1)) - length(head_v);
        head_s = head_v;
        head_sx = head_vx;
        head_sy = head_vy;
        head_a = head_ace;
        head_acex = head_ax;
        head_acey = head_ay;
        clear head_v
        clear head_vx
        clear head_vy
        clear head_ace
        clear head_ax
        clear head_ay


        head_v = zeros(1,length(head_pos(:,1)));
        head_vx = zeros(1,length(head_pos(:,1)));
        head_vy = zeros(1,length(head_pos(:,1)));
        head_ace = zeros(1,length(head_pos(:,1)));
        head_ax = zeros(1,length(head_pos(:,1)));
        head_ay = zeros(1,length(head_pos(:,1)));


        head_v = [nan(dif,1);head_s];
        head_vx = [nan(dif,1); head_sx];
        head_vy = [nan(dif,1); head_sy];
        head_ace = [nan(dif,1); head_a];
        head_ax = [nan(dif,1); head_acex];
        head_ay = [nan(dif,1); head_acey];
        head_mspeed = nanmean(head_v);
        head_mace = nanmean(head_ace);

        
    end
    
    if length(head_x) < length(head_pos(:,1))
    
        head_pos_x = head_x;
        head_pos_y = head_y;
        clear head_x
        clear head_y
        head_x = zeros(1,length(head_pos(:,1)));
        head_y = zeros(1,length(head_pos(:,1)));
        head_x = [nan(dif,1);head_pos_x];
        head_y = [nan(dif,1);head_pos_y];
    end

else
    head_x = [];
    head_y = [];
    head_vx = [];
    head_vy = [];
    head_ax = [];
    head_ay = [];
    head_mace = [];
    head_mspeed = [];
    head_v = [];
    head_ace = [];
end

if ~isempty(tail_pos)
    [t,tail_x,tail_y,~,~,~,~] = trajectory_kalman_filter(tail_pos(:,1)',tail_pos(:,2)',timesamples,0);
    [~,~,~, tail_vx, tail_vy, tail_ax, tail_ay] = KalmanVel(tail_pos(:,1)',tail_pos(:,2)',timesamples,order);
    tail_v = sqrt(tail_vx.^2+tail_vy.^2);
    tail_ace = sqrt(tail_ax.^2+tail_ay.^2);
    if length(tail_v) < length(tail_pos(:,1))
        dif = length(tail_pos(:,1)) - length(tail_v);
        tail_s = tail_v;
        tail_sx = tail_vx;
        tail_sy = tail_vy;
        tail_a = tail_ace;
        tail_acex = tail_ax;
        tail_acey = tail_ay;
        clear tail_v
        clear tail_vx
        clear tail_vy
        clear tail_ace
        clear tail_ax
        clear tail_ay

        tail_v = zeros(1,length(tail_pos(:,1)));
        tail_vx = zeros(1,length(tail_pos(:,1)));
        tail_vy = zeros(1,length(tail_pos(:,1)));
        tail_ace = zeros(1,length(tail_pos(:,1)));

        tail_v = [nan(dif,1);tail_s];
        tail_vx = [nan(dif,1); tail_sx];
        tail_vy = [nan(dif,1); tail_sy];
        tail_ace = [nan(dif,1); tail_a];
        tail_ax = [nan(dif,1); tail_acex];
        tail_ay = [nan(dif,1); tail_acey];
        tail_mspeed = nanmean(tail_v);
        tail_mace = nanmean(tail_ace);

    end
    
    if length(tail_x) < length(pos(:,1))
        tail_pos_x = tail_x;
        tail_pos_y = tail_y;
        clear tail_x
        clear tail_y
        tail_ax = zeros(1,length(tail_pos(:,1)));
        tail_ay = zeros(1,length(tail_pos(:,1)));
        tail_x = [nan(dif,1);tail_pos_x];
        tail_y = [nan(dif,1);tail_pos_y];
    end
else
    tail_x = [];
    tail_y = [];
    tail_vx = [];
    tail_vy = [];
    tail_ax = [];
    tail_ay = [];
    tail_mace = [];
    tail_mspeed = [];
    tail_ace = tail_ace;
    tail_v = tail_v;
end




%% Plotting filtered tracking data

folder = pwd;
folder = strsplit(folder,filesep);
folder = folder{end};

figure,
plot(x,y,'k'), title('Filtered body tracking data')
hold on
plot([boundingbox_xmin boundingbox_xmin boundingbox_xmax boundingbox_xmax boundingbox_xmin], [boundingbox_ymin boundingbox_ymax, boundingbox_ymax, boundingbox_ymin boundingbox_ymin],'r')
if exist('zone','var')
    hold on;
    for i=1:num_zones
    plot([boundingbox_zones_xmin{i} boundingbox_zones_xmin{i} boundingbox_zones_xmax{i} boundingbox_zones_xmax{i} boundingbox_zones_xmin{i}],[boundingbox_zones_ymin{i} boundingbox_zones_ymax{i}, boundingbox_zones_ymax{i}, boundingbox_zones_ymin{i} boundingbox_zones_ymin{i}],'b')
    end  
end
view(0,-90)

saveas(gcf,[apparatus_name,'_',folder,'.png']);


if ~isempty(head_pos)
    figure,
    plot(head_x,head_y,'k'), title('Filtered head tracking data')
    hold on
    plot([boundingbox_xmin boundingbox_xmin boundingbox_xmax boundingbox_xmax boundingbox_xmin], [boundingbox_ymin boundingbox_ymax, boundingbox_ymax, boundingbox_ymin boundingbox_ymin],'r')
    if exist('zone','var')
        hold on;
        for i=1:num_zones
        plot([boundingbox_zones_xmin{i} boundingbox_zones_xmin{i} boundingbox_zones_xmax{i} boundingbox_zones_xmax{i} boundingbox_zones_xmin{i}],[boundingbox_zones_ymin{i} boundingbox_zones_ymax{i}, boundingbox_zones_ymax{i}, boundingbox_zones_ymin{i} boundingbox_zones_ymin{i}],'b')
        end  
    end
end
view(0,-90)

if ~isempty(tail_pos)
    figure,
    plot(head_x,head_y,'k'), title('Filtered tail tracking data')
    hold on
    plot([boundingbox_xmin boundingbox_xmin boundingbox_xmax boundingbox_xmax boundingbox_xmin], [boundingbox_ymin boundingbox_ymax, boundingbox_ymax, boundingbox_ymin boundingbox_ymin],'r')
    if exist('zone','var')
        hold on;
        for i=1:num_zones
        plot([boundingbox_zones_xmin{i} boundingbox_zones_xmin{i} boundingbox_zones_xmax{i} boundingbox_zones_xmax{i} boundingbox_zones_xmin{i}],[boundingbox_zones_ymin{i} boundingbox_zones_ymax{i}, boundingbox_zones_ymax{i}, boundingbox_zones_ymin{i} boundingbox_zones_ymin{i}],'b')
        end  
    end
end
view(0,-90)
% Maybe we need to do something to sync with the digitalIn data

%% TRYING FILTER REFINEMENT
% pos = [xPos,yPos];
% pos_filtered = [x,y];

% position = trackingRefinement(pos,pos_filtered,'pixels_metre',pixels_metre);
% pos = [x,y];
% position = trackingRefinement_v2(pos,'pixels_metre',pixels_metre);

%%
[~,fbasename,~] = fileparts(pwd);

tracking.position.x = x; % x filtered
tracking.position.y = y; % y filtered
tracking.position.z = [];
tracking.position.vx = vx/100;
tracking.position.vy = vy/100;
tracking.position.ax = ax/100;
tracking.position.ay = ay/100;
tracking.position.ace = ace/100;
tracking.position.speed = v/100;
tracking.position.mspeed = mspeed/100;
tracking.position.mace = mace/100;

tracking.description = 'AnyMaze Tracking';
tracking.timestamps = timesamples;
tracking.originalTimestamps = [];
tracking.folder = fbasename;
tracking.sync.sync = [];
tracking.sync.timestamps = [];
tracking.samplingRate = fs;
% Average frame in our pipeline would be the apparatus map
average_frame = [];
tracking.avFrame.r = average_frame;
tracking.avFrame.xSize = xMaze;
tracking.avFrame.ySize = yMaze;

% tracking.roi.roiTracking = roiTracking;
% tracking.roi.roiLED = roiLED;
if exist('apparatus_name','var')
    tracking.apparatus.name = apparatus_name;
end
tracking.apparatus.centre.x = centre_x;
tracking.apparatus.centre.y = centre_y;
tracking.apparatus.boundingbox.xmin = boundingbox_xmin;
tracking.apparatus.boundingbox.xmax = boundingbox_xmax;
tracking.apparatus.boundingbox.ymin = boundingbox_ymin;
tracking.apparatus.boundingbox.ymax = boundingbox_ymax;
% tracking.apparatus.name = apparatus_name;


% head and tail tracking and zones
if ~isempty(head_xPos)
    tracking.headposition.x = head_x;
    tracking.headposition.y = head_y;
    tracking.headposition.z = [];
    tracking.headposition.vx = head_vx/100;
    tracking.headposition.vy = head_vy/100;
    tracking.headposition.ax = head_ax/100;
    tracking.headposition.ay = ay/100;
    tracking.headposition.ace = head_ace/100;
    tracking.headposition.speed = head_v/100;
    tracking.headposition.mspeed = head_mspeed/100;
    tracking.headposition.mace = head_mace/100;
end

if ~isempty(tail_xPos)
    tracking.tailposition.x = tail_x;
    tracking.tailposition.y = tail_y;
    tracking.tailposition.z = [];
    tracking.tailposition.vx = tail_vx/100;
    tracking.tailposition.vy = tail_vy/100;
    tracking.tailposition.ax = tail_ax/100;
    tracking.tailposition.ay = tail_ay/100;
    tracking.tailposition.ace = tail_ace/100;
    tracking.tailposition.speed = tail_v/100;
    tracking.tailposition.mspeed = tail_mspeed/100;
    tracking.tailposition.mace = tail_mace/100;
end

if ~isempty(inArm1)
    tracking.zone.inArm1 = inArm1;
end

if ~isempty(inArm2)
    tracking.zone.inArm2 = inArm2;
end

if ~isempty(inArm2)
    tracking.zone.inArm3 = inArm3;
end

if ~isempty(inCenter)
    tracking.zone.inCenter = inCenter;
end

if ~isempty(inNoZone)
    tracking.zone.inNoZone = inNoZone;
end

if ~isempty(inMain)
    tracking.zone.inMain = inMain;
end

if exist('zone','var')
    for i=1:num_zones
        tracking.zone.name{i} = name_zones{i};
        tracking.zone.centre.x{i} = center_zones_x{i};
        tracking.zone.centre.y{i} = center_zones_y{i};
        tracking.zone.boundingbox.x{i} = boundingbox_zones_x{i};
        tracking.zone.boundingbox.y{i} = boundingbox_zones_y{i};
        tracking.zone.boundingbox.w{i} = boundingbox_zones_w{i};
        tracking.zone.boundingbox.h{i} = boundingbox_zones_h{i};
               
    end
end

tracking.pixelsmetre = pixels_metre;

% Compute distance travelled
distance = computeDistance(tracking.position.x, tracking.position.y);
distance = cumsum(distance(~isnan(distance)));
distance = max(distance);

tracking.distance = distance/100;

if saveMat
    save([basepath filesep fbasename '.Tracking.Behavior.mat'],'tracking');
end

end

