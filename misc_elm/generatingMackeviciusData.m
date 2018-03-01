%% How to load the data
clear all;
DataFolder = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\IsolatePreTut\6991'; 
cnmfeFilePath = fullfile('\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\IsolatePreTut\6991', 'CNMFE_BatchVer.mat');
load(fullfile(DataFolder, 'singinginfo.mat')) 
load(cnmfeFilePath)
%% compile 
L = 1; % units of seconds
SONGfs = 1/diff(singinginfo(1).SpecTime(1:2));
VIDEOfs = singinginfo(1).VIDEOfs;
GCDfs = gcd(SONGfs,VIDEOfs); % for alignment
Lsong = L*SONGfs;
Lneural = L*VIDEOfs;
% for each file
% drop long gaps from spectrogram and neural
% concatenate song and neural, and keep track of new seg times
% do extraction

Fnums = 1:length(singinginfo); 
compSONG = [];
compNEURO = [];

 
for fi = 1:length(Fnums)
    filei = Fnums(fi); 
    RawSong = singinginfo(filei).SongSpec;
    RawNeural =  neuron_batch(filei).DeconvSpiketrain(:,:);
    RawNeural(isnan(RawNeural)) = 0; % deconv sometimes makes nans

    [N,T] = size(RawNeural);
    
    if length(singinginfo(filei).segs)>0
        % take out long gaps, only truncating at GCDfs (so alignment
        % consistent)
        bouts = SegsToBouts(singinginfo(filei).segs,.05*singinginfo(filei).SOUNDfs);
        onsetTimes = ceil(bouts(:,1)*GCDfs/singinginfo(filei).SOUNDfs); 
        offsetTimes = floor(bouts(:,2)*GCDfs/singinginfo(filei).SOUNDfs); 
        
        
        % take out long gaps in neural
        NEURAL1 = [];
        for segi = 1:length(onsetTimes)
            onsetTime =onsetTimes(segi)*VIDEOfs/GCDfs+1; 
            offsetTime =offsetTimes(segi)*VIDEOfs/GCDfs;
            % expand spectrogram at syll onset
            NEURAL1 = [NEURAL1 RawNeural(:,onsetTime:offsetTime)];
        end
        NEURAL = NEURAL1; 
%         FirstOfFileNeuro(fi) = size(compNEURO,2)+1; 
        compNEURO = [compNEURO NEURAL]; 
    
        % take out long gaps in song
        SONG1 = [];
        for segi = 1:length(onsetTimes)
            onsetTime = onsetTimes(segi)*SONGfs/GCDfs+1; 
            offsetTime = offsetTimes(segi)*SONGfs/GCDfs; 
            % expand spectrogram at syll onset
            SONG1 = [SONG1  RawSong(:,onsetTime:offsetTime)];
        end
        SONG = SONG1; 
%         FirstOfFileSong(fi) = size(compSONG,2)+1; 
        compSONG = [compSONG SONG]; 
        figure(5); imagesc(compSONG); drawnow; shg
    end
end
NEURAL = compNEURO; 
SONG = compSONG; 
%
NEURAL = NEURAL./(prctile(NEURAL,100,2)+prctile(NEURAL(:),95));
SONG = compSONG; 
SONG = SONG-prctile(SONG(:),70); SONG(SONG<0)=0; 
SONG = SONG/max(SONG(:));

%%
save('C:\Users\emackev\Documents\MATLAB\seqnmf\MackeviciusData.mat', ...
    'NEURAL', 'SONG', 'SONGfs', 'VIDEOfs')