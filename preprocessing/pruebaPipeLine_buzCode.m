%% Prueba script para correr el pipeline de BuzCode

% El primer paso es convertir los ficheros .xdat en ficheros .dat y añadir
% el fichero xml
% La función coge todas las subcarpetas dentro de la carpeta de la sesión y
% genera los .dat y los xml

% This function renames the files in order to obtain a buzcode logic

% changeFilesName('\\DESKTOP-IORIG9S\data\HPS22\HPS22_130521');

changeFilesName_notTrackingAllFiles('C:\data\HPR21409');

%%
updateExpFolder('\\DISCOVERY_ONE\Sub\HPS23','F:\data\HPS23');
arrangeSessionFolder('C:\data\HPR21409');
bpath= 'C:\data\HPR21409';
createFiles('basepath',bpath);

bpath = 'G:\HPS23\HPS23_040621_sess1';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[],'medianSubstr',false,'runSummary',true);

bpath = 'H:\data\Project_GLUN3\HPS24\HPS24_290621_sess4';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[],'medianSubstr',false,'runSummary',true);

%%
updateExpFolder('\\DESKTOP-IORIG9S\data\HPS23','F:\data\HPS23');
arrangeSessionFolder('F:\data\HPS23');

bpath = 'F:\data\HPS22';
createFiles('basepath',bpath);
bpath = 'F:\data\HPS23';
createFiles('basepath','F:\data\HPS22');

bpath = 'F:\data\HPS22\HPS22_180521_sess15';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS22\HPS22_180521_sess15';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);


bpath = 'F:\data\HPS23\HPS23_170521_sess8';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS22\HPS22_140521_sess12';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);



%%
bpath = 'F:\data\HPS23\HPS23_040621_sess8';
computeSessionSummary_abad('basepath',bpath,'digitalChannelsList',[3,4],'exclude',{'analogPulses','digitalPulses','tMazeBehaviour','YMazeBehaviour'});
% computeSessionSummary('basepath',bpath,'exclude',{'spikes','analogPulses','digitalPulses','downStates','ripples','tMazeBehaviour','linearMazeBehaviour'});

bpath = 'F:\data\HPS22\HPS22_170521_sess14';
computeSessionSummary_abad('basepath',bpath,'exclude',{'spikes','analogPulses','digitalPulses','downStates','ripples','tMazeBehaviour','YMazeBehaviour'});

%% 

bpath = 'F:\data\HPS22\HPS22_210421_sess4';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS22\HPS22_220421_sess5';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS22\HPS22_230421_sess6';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS22\HPS22_270421_sess7';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS22\HPS22_280421_sess8';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS22\HPS22_290421_sess9';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);
%%
bpath = 'F:\data\HPS22\HPS22_300421_sess10';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS23\HPS23_260421_sess1';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS23\HPS23_270421_sess2';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS23\HPS23_290421_sess3';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);


%%
bpath = 'F:\data\HPS23\HPS23_040621_sess8';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);
bpath = 'F:\data\HPS23\HPS23_090621_sess9';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);
bpath = 'F:\data\HPS23\HPS23_110621_sess10';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'G:\HPS23\HPS23_040621_sess1';
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'diffLFPs',true,'analyzeSubSessions',true,'exclude',{'spikes','ripples','thetaModulation','behaviour','spikeTrain','performance'});

bpath = 'F:\data\HPS23\HPS23_090621_sess9';
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'diffLFPs',true,'analyzeSubSessions',true,'exclude',{'spikes','ripples','thetaModulation','behaviour','spikeTrain'});



%%
bpath = 'G:\HPS23\HPS23_040621_sess1';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'diffLFPs',true,'analyzeSubSessions',false,'exclude',{'spikes','ripples','thetaModulation','behaviour','performance','spikeTrain','excel','CellExplorer'});


bpath = 'F:\data\HPS23\HPS23_300421_sess4';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);


%%
bpath = 'F:\data\HPS23\HPS23_020621_sess13';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS23\HPS23_030621_sess14';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);


% HEAD
%% HPS23 090621 Prueba para SubSessions
bpath = 'G:\HPS23\HPS23_090621_sess9';
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'diffLFPS',true,'analyzeSubSessions',true,'exclude',{'spikes','ripples','thetaModulation','behaviour','performance','spikeTrain','CellExplorer'});


%%

changeFilesName_notTrackingAllFiles('H:\data\HPS25\HPS25_050721');
changeFilesName_notTrackingAllFiles('H:\data\HPS25\HPS25_280621');
arrangeSessionFolder('H:\data\HPS25');
updateExpFolder('\\DISCOVERY_ONE\Sub\HPS23','F:\data\HPS23');

bpath= 'H:\data\HPS25';
createFiles('basepath',bpath);
%%
bpath = 'H:\data\HPS25\HPS25_050721_sess1';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);


%% HPS25 - Drugs
bpath = 'F:\data\HPS25\HPS25_240621_sess2';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS25\HPS25_280621_sess3';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'F:\data\HPS25\HPS25_050721_sess4';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);
