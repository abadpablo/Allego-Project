function [pos] = trackingRefinement_v2(position,varargin)
% Refinement of the tracking based on Juan Pablo Quintanilla
% TrackingRefinementPMAZE
%
%
%
%
%
%
%
%


%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'reboundThresholdp',25,@isnumeric);
addParameter(p,'reboundThresholdr',25,@isnumeric);
addParameter(p,'jumpThreshold',60,@isnumeric);
addParameter(p,'pixels_metre',[],@isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
reboundThresholdp = p.Results.reboundThresholdp;
reboundThresholdr = p.Results.reboundThresholdr;
jumpThreshold = p.Results.jumpThreshold;
pixels_metre = p.Results.pixels_metre;

video = 'Test 1.avi';
obj = VideoReader(video);
frame = read(obj,1);
figure,
imagesc(frame)
hold on
plot(position(:,1),position(:,2));

%% RESOLVER SALTOS
position_out = position;
count = 1;
figure,
imagesc(frame)
hold on;
% view(0,-90)

for i=2:length(position(:,1))-1
    difXp = abs(position(i,1)-position(i-1,1));
    difYp = abs(position(i,2)-position(i-1,2));
    Vdistancep(count) = sqrt((difXp^2)+difYp^2);
    if i== 2540
        disp('Hola')
    end
    if Vdistancep(count) > jumpThreshold
        plot(position(1:i,1),position(1:i,2))
        xlim([min(position(:,1)) max(position(:,1))])
        ylim([min(position(:,2)) max(position(:,2))])
        hold on;
        scatter(position(i,1),position(i,2),'r')
%         view(0,-90)
        position_out(i,1) = 0;
        position_out(i,2) = 0;
        position(i,1) = position(i-1,1);
        position(i,2) = position(i-1,2);
    else
        plot(position(1:i,1),position(1:i,2))
        xlim([min(position(:,1)) max(position(:,1))])
        ylim([min(position(:,2)) max(position(:,2))])
%         view(0,-90)
        
        
    end
    count = count+1;
end











end
