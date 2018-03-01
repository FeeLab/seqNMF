%% load data

%% load data
clear all;
DataFolder = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\6992\ForStabilityTest\';
% DataFolder = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\6992SingingAuto\';
% cnmfeFilePath = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\6992SingingAuto\CNMFE_BatchVerMO_4days'
cnmfeFilePath = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\6992\ForStabilityTest\CNMFE_full_moBatchVer_updated_partial'


% DataFolder = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\6938JustSinging';
% cnmfeFilePath = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\6938JustSinging\CNMFE_BatchVer_updated';

% DatenumOfTutoring = datenum('04/06/17 12:00pm')

load(fullfile(DataFolder, 'singinginfo.mat')) 
load(cnmfeFilePath)
neuron_batch = neuron_batch_partial;
%% Keep track of which day each neuron was from
MyDatenum = [];
for fi = 1:length(singinginfo)
    MyDatenum(fi) = (datenum(singinginfo(fi).name((end-18):(end-4)), 'yyyymmdd_HHMMSS'));
%     nSegs(fi) = size(singinginfo(fi).segs,1); 
end
UniqueDays = unique(floor(MyDatenum));
%% compile 
L = 1/3; % units of seconds
SONGfs = 1/diff(singinginfo(1).SpecTime(1:2));
VIDEOfs = singinginfo(1).VIDEOfs;
GCDfs = gcd(SONGfs,VIDEOfs); % for alignment
Lsong = L*SONGfs;
Lneural = L*VIDEOfs;
% for each file
% drop long gaps from spectrogram and neural
% concatenate song and neural, and keep track of new seg times
% do extraction

% Fnums =  find(MyDatenum<DatenumOfTutoring); %[6     7     8     9    10    11    12    20    21    22    23    24    25    26    27    32    33    34    35    36    37]; % just 2 days
% Fnums = 20:37;
Fnums = 1:length(singinginfo); 
compSONG = [];
compNEURO = [];
compFnumNeuro = [];
compFnumSong = [];

 
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
        % BUFFER = nan(size(NEURAL,1),Lneural);%.*mean(SPEC,2); 
        NEURAL1 = [];
        for segi = 1:length(onsetTimes)
            onsetTime =onsetTimes(segi)*VIDEOfs/GCDfs+1; 
            offsetTime =offsetTimes(segi)*VIDEOfs/GCDfs;
            % expand spectrogram at syll onset
            NEURAL1 = [NEURAL1 RawNeural(:,onsetTime:offsetTime)];
        end
        NEURAL = NEURAL1; 
%         NEURAL(NEURAL>=1) = 1; % binarize, and scale for plotting
        NEURAL(:,end-min(size(NEURAL,2)-1, Lneural):end) = 0; % can't fit last L timepoints, so pad with zero
%         FirstOfFileNeuro(fi) = size(compNEURO,2)+1; 
        compNEURO = [compNEURO NEURAL]; 
        compFnumNeuro = [compFnumNeuro filei*ones(1,size(NEURAL,2))]; 
    
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
        compFnumSong = [compFnumSong filei*ones(1,size(SONG,2))]; 
        figure(5); imagesc(compSONG); drawnow; shg
    end
end
UniqueDays = unique(floor(MyDatenum(compFnumNeuro)));
NEURAL = compNEURO; 
SONG = compSONG; 

DatenumOfTutoring = datenum('May 20, 2017, 10:00am');


%% choose which day to look at
% NEURAL = NEURAL(:,MyDatenum(compFnumNeuro)>DatenumOfTutoring+2); 
% SONG = SONG(:,(MyDatenum(compFnumSong)>DatenumOfTutoring+2)); 
NEURAL = NEURAL(:,floor(MyDatenum(compFnumNeuro))==floor(DatenumOfTutoring)+2); 
SONG = SONG(:,floor(MyDatenum(compFnumNeuro))==(floor(DatenumOfTutoring)+2)); 
%% renormalize
NEURAL = NEURAL./(prctile(NEURAL,100,2)+prctile(NEURAL(:),95));
% NEURAL(max(NEURAL,[],2)<.5,:) = [];
% NEURAL = double(NEURAL>0.1);

SONG = SONG-prctile(SONG(:),70); SONG(SONG<0)=0; 
SONG = SONG/max(SONG(:));
%%
figure; subplot(4,1,1); imagesc((SONG)); axis tight
subplot(4,1,2:4); helper.grayscalePatchPlot(NEURAL); axis tight
% dayind = find(MyDatenum(compFnumNeuro)<DatenumOfTutoring); 
% dayind = find(floor(MyDatenum(compFnumNeuro))<DatenumOfTutoring); 
% dayind = find(floor(MyDatenum(compFnumNeuro))>=DatenumOfTutoring); 
%% break data into training set and test set
splitN = floor(size(NEURAL,2)*.75); 
splitS = floor(size(SONG,2)*.75); 
trainNEURAL = NEURAL(:,1:splitN); 
trainSONG = SONG(:,1:splitS); 
testNEURAL = NEURAL(:,(splitN+1):end); 
testSONG = SONG(:,(splitS+1):end); 
%
%% plot one example factorization
% rng(122); % fixed rng seed for reproduceability
X = trainNEURAL;
K = 1;
L =1; % units of seconds
Lneural = L*VIDEOfs;
Lsong = L*SONGfs;
display('Fitting seqNMF on real neural data (from songbird HVC, recorded by Emily Mackevicius, Fee Lab)')
[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .0, 'lambda', .004, 'maxiter', 100, 'showPlot', 1); 
p = .05; % desired p value for factors

display('Testing significance of factors on held-out data')
[pvals,is_significant] = test_significance(testNEURAL,W,p);

W = W(:,is_significant,:); 
H = H(is_significant,:); 

% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 1; % plot data starting at this timebin
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 1,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
title('Significant seqNMF factors, with raw data')
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), ...
    helper.reconstruct(W(indSort,:,:),H(:,tstart:end)),...
    1,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
title('SeqNMF reconstruction')









%% run seqNMF
rng(1000);
% dayind = find(MyDatenum(compFnumNeuro)<DatenumOfTutoring); 
% dayind = find(floor(MyDatenum(compFnumNeuro))<DatenumOfTutoring); 
% dayind = find(floor(MyDatenum(compFnumNeuro))>=DatenumOfTutoring); 

NEURAL = double(compNEURO(:,:)>0);

% dayind = find(MyDatenum(compFnumNeuro)<DatenumOfTutoring); 
% dayind = find(floor(MyDatenum(compFnumNeuro))<DatenumOfTutoring); 
% dayind = find(floor(MyDatenum(compFnumNeuro))<=UniqueDays(end)); 

[N,T] = size(NEURAL);
K = 2;
config = struct();
config.ploton = 1; 
config.maxiter = 100;  
config.lambdaW = .002; 
config.competeW = 0; 
config.lambdaL1W = 1; 
config.lambdaL1H = 0; 
[W,H,errneural] = algorithms.seq_NMF(NEURAL,K,Lneural,config); 
W = W(:,[2 1 ],:);
H = H([2 1],:);
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W,1);
indSort = hybrid(:,3);
compNEURO = compNEURO(indSort,:);
NEURAL = NEURAL(indSort,:);
W = W(indSort,:,:);
grayscaleWHPlot(W,H,NEURAL,SONG,[]); 
%% save data
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\EmilySfNposter';
saveTitle = 'SortAndSaveData'; %READ THIS LINE AND NEXT CAREFULLYYY DONT COPY BEFORE RENAMING
% thisFile =  fullfile(savedir, 'SortAndSaveData.m'); 
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '.pdf'])); 
    save(fullfile(savedir, [saveTitle '.mat']))
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end

readytosave = 0; 
%% check that it works
clear all
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\EmilySfNposter';
saveTitle = 'SortAndSaveData'; 
load(fullfile(savedir, [saveTitle '.mat']))
clf; grayscaleWHPlot(W,H,NEURAL,SONG,[]); 



%% propagate W to each day separately
clf;set(gcf, 'color', 'w')
fac = 1; 
rng(1000);

plotSONG = SONG-prctile(SONG(:),60); 
plotSONG(plotSONG<0) = 0;  
D = length(UniqueDays); 
compileWs = []; 
compileWsShuff = []; 
for dayi = 1:D
    axStart = .1 + (dayi-1)*.8/D; 
    axW = .8/D*.85; 
    axD(dayi) = subplot('Position', [axStart .1 axW .8]);
dayind = find(floor(MyDatenum(compFnumNeuro))==UniqueDays(dayi)); 
A = [];
for iteri = 1:100; %INCREASE IF ACTUALLY SHUFFLING
shuffle = iteri>1;
if shuffle
    V = NEURAL(:,dayind);
    for ni = 1:N
        indEmpty = sum(NEURAL(:,dayind),1)==0;
        V(ni,~indEmpty) = circshift(NEURAL(ni,dayind(~indEmpty)),randi(length(dayind)),2); 
    end
else
    V = NEURAL(:,dayind); 
end
H = zeros(K,length(dayind)); 
for l = 1 : Lneural
    V_shifted = [V(:, l:length(dayind)) zeros(N, l-1)];
    H = H + W(:, :, l)' * V_shifted;
end   
if iteri == 1; 
    [m,indH] =  sort(H(fac,:), 'descend');
    plotS = [];
    if m(1)>0
        for ei = 1:15
            tSample = dayind(1)-1+indH(ei); 
            if (round(tSample*SONGfs/VIDEOfs) +  Lsong)<=size(plotSONG,2)
                plotS = [plotSONG(:,round(tSample*SONGfs/VIDEOfs):(round(tSample*SONGfs/VIDEOfs) +  Lsong)); plotS]; 
            end
        end
    end
    imagesc(plotS, ...
        'xdata',[1 Lneural], 'ydata', [-N -1])
    cmap = jet; 
    cmap(1,:) = 0 
        colormap(cmap)
     set(gca, 'ydir', 'normal')
end
    runningW = [];
    for li = 1:Lneural
        H_shifted = [zeros(1, li-1) H(fac, 1:end-li+1)];
        runningW(:,li) = V * H_shifted';
    end
    runningW = runningW-prctile(runningW,50,2); 
    runningW(runningW<0)=0;


    if iteri==1
        compileWs(dayi,:,:) = runningW;
        grayscalePatchPlot(flipud(runningW)+eps,[])
        xlim([0 Lneural+2]); ylim([-N N])

        axis off;
        shg; drawnow
    else
        compileWsShuff(dayi,iteri,:,:) = runningW;
    end

end
% num2str(A(1)/mean(A(2:end)))
% title(['p\leq ' num2str((sum(A(2:end)>A(1))+1)/iteri(end))], 'interpreter', 'tex')
% drawnow
end

papersize = [10 12];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\EmilySfNposter';
saveTitle = 'TrackSeqPlusSigTesting'; 
thisFile = 'C:\Users\emackev\Documents\MATLAB\Calmaing\mcn_code\misc_elm\TestbedSfN.m'; 
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '.pdf'])); 
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end

for dayi = 1:D
    subplot(axD(dayi));
    xl = xlim; yl = ylim; cla
    xlim([0 Lneural+2]); ylim([-N N])
%     normfac = squeeze(sum(compileWsShuff(dayi,:,:,:),[],2),4)+1);
    wFlat = squeeze(sum(compileWs(dayi,:,:),3));
    wFlatShuff = squeeze(sum(compileWs(dayi,:,:),4));
    normfac = prctile(wFlatShuff(:),100-(5/N/D));
    scatter(fliplr(wFlat)./normfac + Lneural, (1:N)-.5,5,'o','cdata', [0 0 0] + (fliplr(wFlat)'./normfac < 1)*.75, ...
        'markerfacecolor', 'flat'); 
    hold on
    for ni = 1:N
        plot([Lneural ((wFlat(:,N-ni+1))./normfac + Lneural)],(ni-.5)*[1 1], ...
            'color', [0 0 0] + ((wFlat(:,N-ni+1))'./normfac < 1)*.75)
    end
    %     bar(wFlatShuff./normfac, '.','color', (daycolors(dayi,:)+1)/2)
    plot( [1 1]+ Lneural, [1 N ]-.5,'r')
    set(gca, 'color', 'none')% 'ydir', 'reverse')
    axis off
    drawnow; 
%     grayscalePatchPlot(squeeze(sum(compileWs,3))./squeeze(sum(mean(compileWsShuff,2),4)+eps),10)
end

papersize = [10 12];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\EmilySfNposter';
saveTitle = 'TrackSeqPlusSigTesting_auxiliary'; 
thisFile = 'C:\Users\emackev\Documents\MATLAB\Calmaing\mcn_code\misc_elm\TestbedSfN.m'; 
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '.pdf'])); 
%     copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end
readytosave = 0; 
% clf
% clf

% NOTE, COMMENTED OUT HORIZONTAL LINES IN GRAYSCALEWHPLOT SO THAT IT WASN'T
% TOO BIG TO SAVE AS VECTOR GRAPHICS

