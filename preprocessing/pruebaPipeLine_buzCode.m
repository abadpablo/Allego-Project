%% Prueba script para correr el pipeline de BuzCode

% El primer paso es convertir los ficheros .xdat en ficheros .dat y añadir
% el fichero xml
% La función coge todas las subcarpetas dentro de la carpeta de la sesión y
% genera los .dat y los xml

% This function renames the files in order to obtain a buzcode logic

% changeFilesName('\\DESKTOP-IORIG9S\data\HPS22\HPS22_130521');

changeFilesName_notTrackingAllFiles('\\EXPLORER\data2\Project_GLUN3\HPS25\HPS25_240621');

%%
updateExpFolder('\\DISCOVERY_ONE\Sub\HPS23','F:\data\HPS23');
arrangeSessionFolder('F:\data\HPS23');
bpath= 'F:\data\HPS23';
createFiles('basepath',bpath);
bpath = 'G:\HPS23\HPS23_040621_sess1';
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
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'diffLFPs',true,'analyzeSubSessions',true,'exclude',{'spikes','ripples'});



%%
bpath = 'F:\data\HPS23\HPS23_010621_sess12';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'diffLFPs',true,'analyzeSubSessions',false,'exclude',{'spikes','digitalPulses','ripples','tMazeBehaviour','linearMazeBehaviour','YMazeBehaviour','OpenFieldBehaviour','thetaModulation','behaviour','spikeTrain','performance','excel'});


bpath = 'F:\data\HPS23\HPS23_300421_sess4';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);


%%
bpath = 'F:\data\HPS23\HPS23_020621_sess13';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);

bpath = 'G:\HPS23\HPS23_030621_sess14';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);




%% 

% changeFilesName_notTrackingAllFiles('\\EXPLORER\data2\Project_GLUN3\HPS24\HPS24_230621');
% changeFilesName_notTrackingAllFiles('\\EXPLORER\data2\Project_GLUN3\HPS24\HPS24_240621');


updateExpFolder('\\EXPLORER\data2\Project_GLUN3\HPS24','F:\data\HPS24');
arrangeSessionFolder('F:\data\HPS24');
bpath= 'F:\data\HPS24';
createFiles('basepath',bpath);


updateExpFolder('\\EXPLORER\data2\Project_GLUN3\HPS25','F:\data\HPS25');
arrangeSessionFolder('F:\data\HPS25');
bpath= 'F:\data\HPS25';
createFiles('basepath',bpath);

bpath = 'F:\data\HPS24\HPS24_230621_sess1';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);


bpath = 'F:\data\HPS25\HPS25_040621_sess2';
bz_PreprocessSession('basepath',bpath,'getPos',false,'analogCh',[]);
%%
bpath = 'G:\HPS23\HPS23_040621_sess1';
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'diffLFPs',true,'analyzeSubSessions',false,'exclude',{'spikes','ripples','thetaModulation','behaviour','performance','spikeTrain','CellExplorer','lfp_analysis'});


bpath = 'G:\HPS23\HPS23_090621_sess9';
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'diffLFPs',false,'analyzeSubSessions',true,'exclude',{'spikes','ripples','thetaModulation','behaviour'});


%% OPTO
% No cleanArtifacts
bpath = 'D:\opto\HPR21409\HPR214049_230921_sess3';
bz_PreprocessSession('basepath',bpath);
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'digitalChannelsList',7,'exclude',{'spikes','digitalPulses','ripples','thetaModulation','behaviour','performance','spikeTrain','CellExplorer'});

% Clean Artifacts
bpath = 'D:\opto\HPR21409\HPR214049_230921_sess4';
bz_PreprocessSession('basepath',bpath,'cleanArtifacts',{7});
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'digitalChannelsList',7,'exclude',{'behaviour','performance','lfp_analysis'});

bpath = 'D:\opto\HPR21409\HPR21409_240921_sess4';
bz_PreprocessSession('basepath',bpath,'cleanArtifacts',{7});
bz_analyzeSession('basepath',bpath,'getWaveformsFromDat',true,'digitalChannelsList',7,'exclude',{'behaviour','performance','lfp_analysis','excel'});


