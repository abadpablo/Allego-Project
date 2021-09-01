function [pos] = trackingRefinement(position,position_2,varargin)
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
addParameter(p,'jumpThreshold',20,@isnumeric);
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

figure,
imagesc(frame)
hold on
plot(position_2(:,1),position_2(:,2))


%% ELIMINAR COORDENADAS DE PUNTOS CONFUNDIDOS
% figure,
% ptsToCut = readPoints(frame,2);
% 
% for i=1:length(position_2(:,1))
%     if position_2(i,1)>ptsToCut(1,1) && position_2(i,1)<ptsToCut(1,2) &&position_2(i,2)>ptsToCut(2,1) && position_2(i,2)<ptsToCut(2,2)
%         position_2(i,:)=0;
%     end
% end
% close

%% EVALUAR CANTIDAD DE FRAMES SIN POSICIÓN DEL ANIMAL E INTERPOLAR
for i=1:length(position(:,1))
    if isnan(position(i,1))
        lostPosition(i) = 1;
    else
        lostPosition(i) = 0;
    end
end

percentageLost = round(sum(lostPosition)/length(position_2(:,1))*100);
Lost=0;
for i=1:length(lostPosition)
    if lostPosition (i)==1
        Lost=Lost+lostPosition (i);
    elseif lostPosition (i)==0
        Lost=0;
    end
    accumulatedLost(i) =Lost;
end

figure
subplot(2,2,1:2)
plot(lostPosition);xlabel('FRAMES'); ylim([0 2]); xlim([0 length(position(:,1))])
title(['POSITION LOST IN ',num2str(percentageLost),'% OF FRAMES']);

subplot(2,2,3:4)
plot(accumulatedLost);xlabel('FRAMES'); xlim([0, length(position(:,1))])

% INTERPOLAR POSICIÓN EN FRAMES PERDIDOS HASTA X UMBRAL
frameThreshold = 85;
Positionv2 = position;
if lostPosition(1) == 1
    Positionv2(1,:) = position(2,:);
end

for i=2:length(Positionv2(:,1))- frameThreshold %NO se trabaja sobre las últimas coordenadas (tamaño del umbral)
    thresholdVector=i:(i+frameThreshold-1);
       if lostPosition(i)==1 && lostPosition(i-1)==0
         if sum(lostPosition(thresholdVector))<frameThreshold%si hay coordenada en los siguientes 10 frames
              for ii=1:frameThreshold%buscar coordenada siguiente más cercana
                if lostPosition(thresholdVector(ii))==1
                    closerCoordinateTemplate(ii)=0;
                else
                    closerCoordinateTemplate(ii)=1;
                end
              end
               closerCoordinate=min(find(closerCoordinateTemplate));%posiciones con coordenada dentro de los 10 sgtes
               lastCoordinate=position(i-1,:); 
               nextCoordinate=position(i+closerCoordinate-1,:) ;
               
               
              Positionv2(i-1:i+closerCoordinate-1,1)=linspace(position(i-1,1),position(i+closerCoordinate-1,1),closerCoordinate+1); 
              Positionv2(i-1:i+closerCoordinate-1,2)=linspace(position(i-1,2),position(i+closerCoordinate-1,2),closerCoordinate+1);   
         end
       end
end

%% PLOTEAR COORDENADAS DE POSICIÓN SOBRE PRIMER FRAME DE VIDEO DE CONDUCTA
close all;
figure,
imagesc(frame);hold on;plot(Positionv2(:,1),Positionv2(:,2))
%% EVALUAR CANTIDAD DE FRAMES SIN POSICIÓN DEL ANIMAL E INTERPOLAR

for i=1:length(position_2(:,1))
    if position_2(i,1) == 0
        lostPosition(i) = 1;
    else
        lostPosition(i) = 0;
    end
end

percentageLost = round(sum(lostPosition)/length(position_2(:,1))*100);
Lost=0;
for i=1:length(lostPosition)
    if lostPosition (i)==1
        Lost=Lost+lostPosition (i);
    elseif lostPosition (i)==0
        Lost=0;
    end
    accumulatedLost(i) =Lost;
end

figure
subplot(2,2,1:2)
plot(lostPosition);xlabel('FRAMES'); ylim([0 2]); xlim([0 length(position_2(:,1))])
title(['POSITION LOST IN ',num2str(percentageLost),'% OF FRAMES']);

subplot(2,2,3:4)
plot(accumulatedLost);xlabel('FRAMES'); xlim([0, length(position_2(:,1))])

% INTERPOLAR POSICIÓN EN FRAMES PERDIDOS HASTA X UMBRAL
frameThreshold = 85;
Positionv2 = position_2;
if lostPosition(1) == 1
    Positionv2(1,:) = position_2(2,:);
end

for i=2:length(Positionv2(:,1))- frameThreshold %NO se trabaja sobre las últimas coordenadas (tamaño del umbral)
    thresholdVector=i:(i+frameThreshold-1);
       if lostPosition(i)==1 && lostPosition(i-1)==0
         if sum(lostPosition(thresholdVector))<frameThreshold%si hay coordenada en los siguientes 10 frames
              for ii=1:frameThreshold%buscar coordenada siguiente más cercana
                if lostPosition(thresholdVector(ii))==1
                    closerCoordinateTemplate(ii)=0;
                else
                    closerCoordinateTemplate(ii)=1;
                end
              end
               closerCoordinate=min(find(closerCoordinateTemplate));%posiciones con coordenada dentro de los 10 sgtes
               lastCoordinate=Position(i-1,:); 
               nextCoordinate=Position(i+closerCoordinate-1,:) ;
               
               
              Positionv2(i-1:i+closerCoordinate-1,1)=linspace(Position(i-1,1),Position(i+closerCoordinate-1,1),closerCoordinate+1); 
              Positionv2(i-1:i+closerCoordinate-1,2)=linspace(Position(i-1,2),Position(i+closerCoordinate-1,2),closerCoordinate+1);   
         end
       end
end

%% PLOTEAR COORDENADAS DE POSICIÓN SOBRE PRIMER FRAME DE VIDEO DE CONDUCTA
close all;
figure,
imagesc(frame);hold on;plot(Positionv2(:,1),Positionv2(:,2))


%% RESOLVER REBOTES

Positionv3=Positionv2;
for i=2:length(Positionv3(:,1))-1
    difXp=abs(Positionv3(i,1)-Positionv3(i-1,1));
    difYp=abs(Positionv3(i,2)-Positionv3(i-1,2));
    Vdistancep=sqrt((difXp^2)+(difYp^2));
    
    difXr=abs(Positionv3(i+1,1)-Positionv3(i-1,1));
    difYr=abs(Positionv3(i+1,2)-Positionv3(i-1,2));
    Vdistancer=sqrt((difXr^2)+(difYr^2));
    if Vdistancep>reboundThresholdp && Vdistancer<reboundThresholdr
        Positionv3(i,1)=mean([Positionv3(i-1,1), Positionv3(i+1,1)]);
        Positionv3(i,2)=mean([Positionv3(i-1,2), Positionv3(i+1,2)]);
    end    
end

%% PLOTEAR COORDENADAS DE POSICIÓN SOBRE PRIMER FRAME DE VIDEO DE CONDUCTA
close all;
figure,
imagesc(frame);hold on;plot(Positionv3(:,1),Positionv3(:,2))


%% RESOLVER SALTOS
position_out = Positionv3;
for i=2:length(Positionv3(:,1))-1
    difXp=abs(Positionv3(i,1)-Positionv3(i-1,1));
    difYp=abs(Positionv3(i,2)-Positionv3(i-1,2));
    Vdistancep=sqrt((difXp^2)+(difYp^2));
    
    difXr=abs(Positionv3(i+1,1)-Positionv3(i-1,1));
    difYr=abs(Positionv3(i+1,2)-Positionv3(i-1,2));
    Vdistancer=sqrt((difXr^2)+(difYr^2));
    
    if Vdistancep>jumpThreshold
        Positionv3(i,1)= Positionv3(i-1,1);
        Positionv3(i,2)= Positionv3(i-1,2);
        frames_remove(i) = i;
        position_out(i,1) = 0;
        position_out(i,2) = 0;
%         figure,
%         plot(Positionv3(:,1),Positionv3(:,2))
%         view(0,-90)
%         hold on
%         scatter(Positionv3(i,1),Positionv3(i,2),'r')
    end    
end

frames_remove2 = frames_remove;
frames_remove2(frames_remove2 == 0) = [];
%Corregir ùltimo frame
Positionv3(length(Positionv3),:)=Positionv3(length(Positionv3)-1,:);




%% PLOTEAR COORDENADAS DE POSICIÓN SOBRE PRIMER FRAME DE VIDEO DE CONDUCTA
close all;
figure,
imagesc(frame);hold on;plot(Positionv3(:,1),Positionv3(:,2));
figure,
imagesc(frame); hold on;
plot(position_2(:,1),position_2(:,2))
%% PLOTEAR NUEVAS POSICIONES SOBRE CADA FRAME







end
