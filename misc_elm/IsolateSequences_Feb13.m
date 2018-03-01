%% How to load the data, deconvolve it
clear all;
DataFolder = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\IsolatePreTut\6991'; 
cnmfeFilePath = fullfile('\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\IsolatePreTut\6991', 'CNMFE_BatchVer.mat');
load(fullfile(DataFolder, 'singinginfo.mat')) 
load(cnmfeFilePath)
%% compile 
L =1; % units of seconds
SONGfs = 1/diff(singinginfo(1).SpecTime(1:2));
VIDEOfs = singinginfo(1).VIDEOfs;
Lsong = L*SONGfs;
Lneural = L*VIDEOfs;

% for each file
% drop long gaps from spectrogram and neural
% concatenate song and neural, and keep track of new seg times
% do extraction

Fnums = 1:length(neuron_batch); 
compSONG = [];
compNEURO = [];
FirstOfFileNeuro = [];
FirstOfFileSong = [];
for fi = 1:length(Fnums)
    filei = Fnums(fi); 
    RawSong = singinginfo(filei).SongSpec;
    RawNeural =  neuron_batch(filei).DeconvSpiketrain(:,:);
    RawNeural(isnan(RawNeural)) = 0; % deconv sometimes makes nans
    [N,T] = size(RawNeural);
    
    if length(singinginfo(filei).segs)>0
        % take out long gaps in neural
        bouts = SegsToBouts(singinginfo(filei).segs,.05*singinginfo(filei).SOUNDfs);
        onsetTimes = floor(bouts(:,1)*singinginfo(filei).VIDEOfs/singinginfo(filei).SOUNDfs); 
        offsetTimes = ceil(bouts(:,2)*singinginfo(filei).VIDEOfs/singinginfo(filei).SOUNDfs); 
        % BUFFER = nan(size(NEURAL,1),Lneural);%.*mean(SPEC,2); 
        NEURAL1 = [];
        for segi = 1:size(bouts,1)
            onsetTime =max(onsetTimes(segi),1); 
            offsetTime =min(offsetTimes(segi),size(RawNeural,2));
            % expand spectrogram at syll onset
            NEURAL1 = [NEURAL1 RawNeural(:,onsetTime:offsetTime)];
        end
        NEURAL = NEURAL1; 
%         NEURAL(NEURAL>=1) = 1; % binarize, and scale for plotting
%         NEURAL(:,end-L:end) = 0; % can't fit last L timepoints, so pad with zero
        FirstOfFileNeuro(fi) = size(compNEURO,2)+1; 
        compNEURO = [compNEURO NEURAL]; 
    
        % take out long gaps in song
        onsetTimes = ceil((onsetTimes)/singinginfo(filei).VIDEOfs/diff(singinginfo(filei).SpecTime(1:2))); 
        offsetTimes = floor((offsetTimes+1)/singinginfo(filei).VIDEOfs/diff(singinginfo(filei).SpecTime(1:2))); 
        SONG1 = [];
        for segi = 1:size(bouts,1)
            onsetTime =max(onsetTimes(segi),1); 
            offsetTime =min(offsetTimes(segi), size(RawSong,2)); 
            % expand spectrogram at syll onset
            SONG1 = [SONG1  RawSong(:,onsetTime:offsetTime)];
        end
        SONG = SONG1; 
        FirstOfFileSong(fi) = size(compSONG,2)+1; 
        compSONG = [compSONG SONG]; 
        figure(5); imagesc(compSONG); drawnow; shg
    end
end
NEURAL = compNEURO; 
%
NEURAL = NEURAL./(prctile(NEURAL,100,2)+prctile(NEURAL(:),95));
NEURAL(max(NEURAL,[],2)<.5,:) = [];
SONG = compSONG; 
SONG = SONG-prctile(SONG(:),70); SONG(SONG<0)=0; 
SONG = SONG/max(SONG(:));
%% break data into training set and test set
splitN = floor(size(NEURAL,2)*.75); 
splitS = floor(size(SONG,2)*.75); 
trainNEURAL = NEURAL(:,1:splitN); 
trainSONG = SONG(:,1:splitS); 
testNEURAL = NEURAL(:,(splitN+1):end); 
testSONG = SONG(:,(splitS+1):end); 


%% choosing lambda
K = 10; 
X = trainNEURAL;
lambdas = sort([logspace(-1,-5,100)], 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    X = trainNEURAL; 
    [N,T] = size(X);
    [W, H, ~,loadings(li,:),power]= seqNMF(X,'K',K,'L',Lneural,...
        'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    li
end
%% plot costs as a function of lambda
R = (regularization-min(regularization)); R = R./mean(R); 
C = (cost -min(cost)); C = C./mean(C);

clf; hold on
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Regularization cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')

savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'IsolateSequences_Feb13'; 
thisFile = 'C:\Users\emackev\Documents\MATLAB\seqnmf\TestbedFigs_elm.m'; 
papersize = [4 3];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_raw_costvslambda.pdf'])); 
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end
readytosave = 0
%% chose lambda=.004; run multiple times, see loadings
loadings = [];
pvals = []; 
is_significant = []; 
X = trainNEURAL;
for iteri = 1:100
    [W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .1, 'lambda', .004, 'maxiter', 100, 'showPlot', 0); 
    p = .01;
    [pvals(iteri,:),is_significant(iteri,:)] = test_significance(testNEURAL,W,p, 5000)
    W = W(:,is_significant(iteri,:)==1,:); 
    H = H(is_significant(iteri,:)==1,:); 
    [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
    indSort = hybrid(:,3);
    tstart = 300; 
    clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
    iteri
end
%% plot significant loadings
clf; hold on
h = histogram(sum(is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]); 
h.BinCounts = h.BinCounts/sum(h.BinCounts)*100; 
xlim([0 10]); 
xlabel('# significant factors')
ylabel('% seqNMF runs')
% %%

% XX = repmat(1:size(loadings,1),size(loadings,1),1)+.25*rand(size(loadings,1))-.125;
% XX = XX(:); 
% LL = loadings(:); 
% 
% scatter(XX(is_significant(:)==1),LL(is_significant(:)==1), 'k*')
% scatter(XX(is_significant(:)==0),LL(is_significant(:)==0), 'ko')
% xlabel('Sorted factor #'); ylabel('Loading')

savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'IsolateSequences_Feb13'; 
thisFile = 'C:\Users\emackev\Documents\MATLAB\seqnmf\TestbedFigs_elm.m'; 
papersize = [3 2];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_NumSigFacs.pdf'])); 
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end
readytosave = 0

%% plot one example factorization
rng(122); %
X = trainNEURAL;
[W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .1, 'lambda', .004, 'maxiter', 100, 'showPlot', 1); 
p = .01;
[pvals,is_significant] = test_significance(testNEURAL,W,p, 5000)
W = W(:,is_significant,:); 
H = H(is_significant,:); 
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 300; 
clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
%
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'IsolateSequences_Feb13'; 
thisFile = 'C:\Users\emackev\Documents\MATLAB\seqnmf\TestbedFigs_elm.m'; 
papersize = [8 6];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_RawSorted.pdf'])); 
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end

clf; WHPlot(W(indSort,:,:),H(:,tstart:end), [], 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'IsolateSequences_Feb13'; 
thisFile = 'C:\Users\emackev\Documents\MATLAB\seqnmf\TestbedFigs_elm.m'; 
papersize = [8 6];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_Reconstruction.pdf'])); 
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end

indperm = indSort(randperm(N));
clf; WHPlot(W(indperm,:,:),H(:,tstart:end), X(indperm,tstart:end), 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'IsolateSequences_Feb13'; 
thisFile = 'C:\Users\emackev\Documents\MATLAB\seqnmf\TestbedFigs_elm.m'; 
papersize = [8 6];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_RawUnsorted.pdf'])); 
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end
readytosave = 0
% %% run seqNMF on neural data training set
% rng(111)
% X = trainNEURAL; 
% [N,T] = size(X);
% K = 3;
% clf; shg
% [W,H,errneural] = seqNMF(X,'K',K,'L',Lneural,...
%     'lambdaL1W', .1, 'lambda', .01, 'maxiter', 150); 
% [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
% indSort = hybrid(:,3);
% 
% tstart = 300; 
% clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
% % clf; WHPlot(W(indSort,:,:),H, X(indSort,:), 0,trainSONG)
% savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
% saveTitle = 'IsolateSequences_Feb13'; 
% thisFile = 'C:\Users\emackev\Documents\MATLAB\seqnmf\TestbedFigs_elm.m'; 
% papersize = [6 4];
% set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
% if readytosave
%     set(gcf, 'color', 'none')
%     saveas(gcf, fullfile(savedir, [saveTitle '_raw_sorted.pdf'])); 
%     copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
% end
% readytosave = 0


% grayscaleWHPlot(WTWO,HTWO, double(NEURAL>0), SONG, VIDEOfs/SONGfs)

% %% testing significance
% p = .05;
% [pvals,is_significant] = test_significance(testNEURAL,W,p, 5000)

%% Factor-triggered song examples and rastors
HTriggeredSpec(H,trainSONG,VIDEOfs,SONGfs,Lsong); 

thisFile = 'C:\Users\emackev\Documents\MATLAB\seqnmf\TestbedFigs_elm.m'; 
papersize = [3 2];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_HtriggeredSpec.pdf'])); 
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end


HTriggeredRaster(H,trainNEURAL(indSort,:),Lneural)
thisFile = 'C:\Users\emackev\Documents\MATLAB\seqnmf\TestbedFigs_elm.m'; 
papersize = [3 3];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
if readytosave
    set(gcf, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_HtriggeredRaster.pdf'])); 
    copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
end

readytosave = 0

% 
% 
% 
% 
% 
% HTriggeredRaster(HTWO,NEURAL(indSort,:),Lneural); 
% papersize = [3 7]; %[185.5 116.17*3/4]/30;
% set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
% if readytosave
%     set(gcf, 'color', 'none')
%     saveas(gcf, fullfile(savedir, [saveTitle '_rasters.pdf'])); 
%     copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
% end
% clf
% HTriggeredSpec(HTWO,SONG,VIDEOfs,SONGfs,Lsong)
% papersize = [3 4]; %[185.5 116.17*3/4]/30;
% set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
% if readytosave
%     set(gcf, 'color', 'none')
%     saveas(gcf, fullfile(savedir, [saveTitle '_specs.pdf'])); 
%     copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
%     set(gcf, 'color', 'w')
% end
% readytosave = 0;
% %%
% helper.AutoSelectLambda(NEURAL, 6, Lneural)
% papersize = [4 3];
% set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
% if readytosave
%     set(gcf, 'color', 'none')
%     saveas(gcf, fullfile(savedir, [saveTitle '_choosinglambda.pdf'])); 
%     copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
%     set(gcf, 'color', 'w')
% end
% readytosave = 0;
% %% compute loadings for multiple runs
% rng(3333)
% 
% [N,T] = size(NEURAL);
% K = 15;
% config = struct();
% config.ploton = 0; 
% config.maxiter = 100;  
% config.lambdaW = .002; 
% config.competeW = 0; 
% config.lambdaL1W = .25; 
% config.lambdaL1H = 0; 
% config.shift = 1;
% clf; shg
% for iteri = 1:10
%     config.lambdaW = .002; 
%     config.lambdaL1W = .25; 
% %     [Wneural,Hneural,errneural] = algorithms.seq_NMF(NEURAL,K,Lneural,config); 
% %     LOADING(iteri,:) = helper.computeLoadingPercentPower(NEURAL,Wneural,Hneural); 
%     config.lambdaW = 0; 
%     config.lambdaL1W = 0; 
%     [Wneural,Hneural,errneural] = algorithms.seq_NMF(NEURAL,K,Lneural,config); 
%     LOADING_cnmf(iteri,:) = helper.computeLoadingPercentPower(NEURAL,Wneural,Hneural); 
%     iteri
% end
% %%
% clf
% X = repmat(1:K,iteri,1) + .5*rand(iteri,K)-.25; 
% scatter(X(:), LOADING_cnmf(:), 'k', 'markerfacecolor', 'flat')
% hold on
% scatter(X(:), LOADING(:), 'r', 'markerfacecolor', 'flat')
% legend('CNMF', 'seqNMF')
% xlabel('Factor #')
% ylabel('Loading (% power explained)'); 
% %
% papersize = [3 2.5];
% set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
% if readytosave
%     set(gcf, 'color', 'none')
%     saveas(gcf, fullfile(savedir, [saveTitle '_Loadings.pdf'])); 
%     copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
%     set(gcf, 'color', 'w')
% end
% readytosave = 0; 
% %% compute rmse & pow explained for multiple K's
% rng(3333)
% 
% [N,T] = size(NEURAL);
% K = 10;
% config = struct();
% config.ploton = 0; 
% config.maxiter = 100;  
% config.lambdaW = .002; 
% config.competeW = 0; 
% config.lambdaL1W = .25; 
% config.lambdaL1H = 0; 
% config.shift = 1;
% clf; shg
% RMSE = []; 
% RMSE_cnmf = [];
% POWEXP = []; 
% POWEXP_cnmf = [];
% W = {};
% H = {};
% W_cnmf = {};
% H_cnmf = {};
% Score = [];
% Score_cnmf = [];
% 
%     
% for Ki = 1:K
%     for iteri = 1:10
%         config.lambdaW = .002; 
%         config.lambdaL1W = .25; 
%         [W{iteri,Ki},H{iteri,Ki},RMSE(iteri,Ki),~,POWEXP(iteri,Ki)] = algorithms.seq_NMF(NEURAL,Ki,Lneural,config); 
%         config.lambdaW = 0; 
%         config.lambdaL1W = 0; 
%         [W_cnmf{iteri,Ki},H_cnmf{iteri,Ki},RMSE_cnmf(iteri,Ki),~,POWEXP_cnmf(iteri,Ki)] = algorithms.seq_NMF(NEURAL,Ki,Lneural,config); 
%         iteri
%         RMSE
%         RMSE_cnmf
%     end
%     Score(:,Ki) = helper.consistency(W(:,Ki), H(:,Ki)); 
%     Score_cnmf(:,Ki) = helper.consistency(W_cnmf(:,Ki), H_cnmf(:,Ki)); 
% end
% %% plot it
% clf
% subplot(1,3,1); cla
% X = repmat(1:K,iteri,1) + .5*rand(iteri,K)-.25; 
% scatter(X(:), RMSE_cnmf(:), 5,'k', 'markerfacecolor', 'flat')
% hold on;
% scatter(X(:), RMSE(:), 5,'r', 'markerfacecolor', 'flat')
% % legend('CNMF', 'seqNMF')
% axis tight; 
% ylims = ylim;
% ylim([0 ylims(2)])
% xlabel('K')
% ylabel('RMSE'); 
% 
% subplot(1,3,2); cla
% X = repmat(1:K,iteri,1) + .5*rand(iteri,K)-.25; 
% scatter(X(:), POWEXP_cnmf(:), 5,'k', 'markerfacecolor', 'flat')
% hold on
% scatter(X(:), POWEXP(:), 5,'r', 'markerfacecolor', 'flat')
% % legend('CNMF', 'seqNMF')
% axis tight; 
% ylims = ylim;
% ylim([0 ylims(2)])
% xlabel('K')
% ylabel('Power explained (%)'); 
% %
% subplot(1,3,3); cla
% patches = plot.violin(Score_cnmf+rand(size(Score))*.00001, 'facecolor','k');
% for ki = 1:K
%     patches(ki).Vertices(:,1) = patches(ki).Vertices(:,1) -.25;
% end
% hold on
% plot.violin(Score+rand(size(Score))*.00001, 'facecolor','r');
% ylims = ylim;
% ylim([0 ylims(2)])
% % X = repmat(1:K,size(Score,1),1) + .5*rand(size(Score,1),K)-.25; 
% % scatter(X(:), Score_cnmf(:), 2,'k', 'markerfacecolor', 'flat')
% % hold on
% % scatter(X(:), Score(:), 2,'r', 'markerfacecolor', 'flat')
% % legend('CNMF', 'seqNMF')
% legend off; box off
% xlabel('K')
% ylabel('Consistency'); 
% 
% papersize = [8 2];
% set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
% if readytosave
%     set(gcf, 'color', 'none')
%     saveas(gcf, fullfile(savedir, [saveTitle '_chooseK.pdf'])); 
%     copyfile(thisFile, fullfile(savedir, [saveTitle '.m'])); 
%     set(gcf, 'color', 'w')
% end
% readytosave = 0; 