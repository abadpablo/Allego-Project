close all; clear all; clc
data=importdata('GC2LINEARTRACKEPIPRE3 10072021 REGISTRODLC_resnet50_PMAZEtrackingJul21shuffle1_300000.csv');%ARCHIVO OUTPUT DLC
video='GC2LINEARTRACKEPIPRE3 10072021 REGISTRO.avi' %video de conducta de la sesión
likehoodDLC=0.9;% punto de corte likelihood dlc
reboundThresholdp=25;
reboundThresholdr=25;
jumpThreshold=80;%80
fileName='tracking'


data=data.data;
obj = VideoReader(video);
frame = read(obj,1);
back=data(:,(5:7));back(back(:, 3)<likehoodDLC, :)= 0;back(:,3)= [];
 neck=data(:,(2:4));neck(neck(:, 3)<likehoodDLC, :)= 0;neck(:,3) = [];
% for i=1:length(back(:,1))
%     if back(i,1)==0
%         back(i,:)=neck(i,:);
%     end
% end
Position=back;
%% ELIMINAR EQUIVALENTES A DROPPED FRAMES DE MINISCOPE SI LOS HAY
toErase=[]
Position(toErase,:) = [];
%% PLOTEAR COORDENADAS DE POSICIÓN SOBRE PRIMER FRAME DE VIDEO DE CONDUCTA
close all;
imagesc(frame);hold on;plot(Position(:,1),Position(:,2))


%% ELIMINAR CORDENADAS DE PUNTOS CONFUNDIDOS
ptsToCut = readPoints(frame,2)
for i=1:length(Position(:,1))
    if Position(i,1)>ptsToCut(1,1) && Position(i,1)<ptsToCut(1,2) &&Position(i,2)>ptsToCut(2,1) && Position(i,2)<ptsToCut(2,2)
        Position(i,:)=0;
    end
end
close
%% EVALUAR CANTIDAD DE FRAMES SIN POSICIÓN DEL ANIMAL e interpolar
for i=1:length(Position(:,1))% crear vector con frames sin posición
    if Position(i,1)==0
        lostPosition(i)=1;
    else
        lostPosition(i)=0;
    end
end
percentageLost=round(sum(lostPosition)/length(Position(:,1))*100);%PORCENTAJE PERDIDO
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
plot(lostPosition);xlabel('FRAMES'); ylim([0 2]); xlim([0 length(Position(:,1))])
title(['POSITION LOST IN ',num2str(percentageLost),'% OF FRAMES']);

subplot(2,2,3:4)
plot(accumulatedLost);xlabel('FRAMES'); xlim([0, length(Position(:,1))])

% INTERPOLAR POSICIÓN EN FRAMES PERDIDOS HASTA X UMBRAL
frameThreshold=85;% 10 equivale a medio segundo
Positionv2=Position;
if lostPosition(1)==1
   Positionv2(1,:)=Position(2,:);
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
imagesc(frame);hold on;plot(Positionv3(:,1),Positionv3(:,2))

%% RESOLVER SALTOS
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
end    
end

%Corregir ùltimo frame
Positionv3(length(Positionv3),:)=Positionv3(length(Positionv3)-1,:);

%% PLOTEAR COORDENADAS DE POSICIÓN SOBRE PRIMER FRAME DE VIDEO DE CONDUCTA
close all;
imagesc(frame);hold on;plot(Positionv3(:,1),Positionv3(:,2));
%% PLOTEAR NUEVAS POSICIONES SOBRE CADA FRAME
for k=1:floor(length(data(:,1))/1000)

        for i=(1000*k)+[-999:10:1]
            frame = read(obj,i);imagesc(frame);hold on;
            scatter(Positionv3(i,1),Positionv3(i,2),5, '*','r');set(gca, 'YDir','reverse');hold on
            title(i)
            if i>50
                axish=i-[(1:50)];
                plot(Positionv3(axish,1),Positionv3(axish,2), 'g')
            end
             pause(0.00015)
        end
       close
end
%% GENERAR VIDEO PARA REVISAR TRACKING
writerObj = VideoWriter(fileName);  writerObj.FrameRate = Fs;open(writerObj);

        for i=1:2:5000
            frame = read(obj,i);imagesc(frame);hold on;
            scatter(Positionv3(i,1),Positionv3(i,2),5, '*','r');set(gca, 'YDir','reverse');hold on
            if i>50
                axish=i-[(1:50)];
                plot(Positionv3(axish,1),Positionv3(axish,2), 'g');set(gca,'xtick',[]);set(gca,'ytick',[]);
            end
            % pause(0.0015)
            frame=getframe(gcf);    writeVideo(writerObj, frame);
        end

close(writerObj);
%% OBTENER COORDENADA DE ESQUINA SUPERIOR IZQUIERDA DE CAMPO DE EXPLORACIÓN Y 
%AJUSTAR COORDENADAS A ESTE PUNTO COMO 0,0
[ceroXY]=readPoints(frame,1)
save('ceroXYneck','ceroXY')
Positionv4(:,1)=Positionv3(:,1)-ceroXY(1);
Positionv4(:,2)=Positionv3(:,2)-ceroXY(2);

%% PLOTEAR COORDENADAS DE POSICIÓN 
close all;
plot(Positionv4(:,1),Positionv4(:,2));set(gca, 'YDir','reverse');set(gcf, 'color', 'white')

 %% DOWNSAMPLEAR
Positionv5=Positionv4(1:2:length(Positionv4),1:2);
%% PLOTEAR COORDENADAS DE POSICIÓN 
close all;
plot(Positionv5(:,1),Positionv5(:,2));set(gca, 'YDir','reverse');set(gcf, 'color', 'white')

%% GUARDAR
Position=Positionv5;
Position=Position(1:end-2,:);
save(fileName,'Position', 'sigfn')
tracking=Position;
save('trackingBACK','tracking')

