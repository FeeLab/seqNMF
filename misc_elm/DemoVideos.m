%% Generate the synthetic data
% number_of_seqences = 4;
% T = 5000; % length of data to generate
% Nneurons = [10,10,10,10]'; % number of neurons in each sequence
% Dt = 3.*ones(number_of_seqences,1); % gap between each member of the sequence
% NeuronNoise = 0.000; 
% SeqNoiseTime = zeros(number_of_seqences,1); % Jitter parameter = 0%
% SeqNoiseNeuron = 1.*ones(number_of_seqences,1); % Participation parameter = 100%
% NEURAL DATA (ISOLATE W STUTTERED SEQUENCES
clear all; 
DataFolder = '\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\IsolatePreTut\7187'; 
cnmfeFilePath = fullfile('\\feevault\data0-shared\shared\EmilyShijieShared\CaExtraction\IsolatePreTut\7187', 'CNMFE_BatchVer.mat');
load(fullfile(DataFolder, 'singinginfo.mat')) 
load(cnmfeFilePath)

%%
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
NEURAL = NEURAL./(prctile(NEURAL,100,2)+prctile(NEURAL(:),99)+eps);
% NEURAL(max(NEURAL,[],2)<.5,:) = [];
SONG = compSONG; 
SONG = SONG-prctile(SONG(:),70); SONG(SONG<0)=0; 
SONG = SONG/max(SONG(:));
%% break into training set

if size(NEURAL,2) == 2098 % 6961 first files are super faint signal
    NEURAL = NEURAL(:,(floor(size(NEURAL,2)*.75)+1):end); 
    SONG = SONG(:,(floor(size(SONG,2)*.75)+1):end); 
end

splitN = floor(size(NEURAL,2)*.75); 
splitS = floor(size(SONG,2)*.75); 
trainNEURAL = NEURAL(:,1:splitN); 
trainSONG = SONG(:,1:splitS); 
testNEURAL = NEURAL(:,(splitN+1):end); 
testSONG = SONG(:,(splitS+1):end); 
%%
% rng(123); %
X = trainNEURAL;
K = 3; 
shg

[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .25, 'lambda', .01, 'maxiter', 100, 'showPlot', 1); 
p = .01;
[pvals,is_significant] = test_significance(testNEURAL,W,p);
W = W(:,is_significant,:); 
H = H(is_significant,:); 
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
%
tstart = 321; 
clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end))
shg
%
%%
136;
rng(seed)
savepath = fullfile('C:\Users\emackev\Downloads\tmpPlots\', ['seqNMFDemo' num2str(seed) '_' datestr(now, 'mmddyyyyMMSS'),'.avi']);
obj = vision.VideoFileWriter(savepath, 'AudioInputPort', 0);%,  'fps', 20);
obj.FrameRate = 4; 

K = 3; L = Lneural; 
X = trainNEURAL(randperm(N),:);
[N,T] = size(X);
tstart = 321; 
clf; WHPlot(zeros(N,K,L),zeros(K,T-tstart), X(:,tstart:end), 0,[],0); drawnow
saveas(gcf, [savepath(1:end-4) '_pre.pdf']);
cdata = print(gcf, '-RGBImage', '-r300');
step(obj, cdata)
[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .2, 'lambda', .01, 'maxiter', 1, 'showPlot', 0); 

for iteri = 1:100
    if 1; iteri<60
        [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
        indSort = hybrid(:,3);
        X = X(indSort,:); 
        W = W(indSort,:,:); 
    end
    clf; WHPlot(W(:,:,:),H(:,tstart:end), X(:,tstart:end), 0,[],0)
    title(['Iteration ' num2str(iteri)]); drawnow
    cdata = print(gcf, '-RGBImage', '-r300');
    step(obj, cdata)
    [W, H]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .2, 'lambda', .01, 'maxiter', 1, 'showPlot', 0,...
        'W_init', W, 'H_init', H, 'SortFactors', 0);
end
saveas(gcf, [savepath(1:end-4) '_post.pdf']);
clf; WHPlot(W(:,:,:),H(:,tstart:end), X(:,tstart:end), 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end),0)
saveas(gcf, [savepath(1:end-4) '_postPlusSpec.pdf']);

release(obj)



%%  SPECTROGRAMS
% add feelab_common_code to path for spectrogramELM
% load an example of labeled data
clear all
rng(1000); % set seed for exact reproducability
pathname = '\\feevault\data0\AcqGui\TotsSomeDataFromCrcns';
datafolder = fullfile(pathname, 'to3567\2013-03-19\HVCp-1'); 
analysisfile =  fullfile(pathname, 'Analysis_files\HVCp\to3567_2013-03-19_HVCp-1');
SoundFiles = dir(fullfile(datafolder, '*chan0.mat')); 
load(analysisfile)
SONG = [];
for fi = [1:2];%length(SoundFiles)
    clf
    % plot the spectrogram
    A = load(fullfile(datafolder, SoundFiles(fi).name))
    [S,Time,F] = spectrogramELM(A.sound,A.fs, .005, 1,[500 8000]); 
    Spec = log(S); Spec = Spec-prctile(Spec(:),70); Spec(Spec<0)=0; 
    clims = [0 prctile(Spec(:),99.5)]; 
    cmap = 1/256*flipud([158,1,66;213,62,79;244,109,67;253,174,97;254,224,139;255,255,191;230,245,152;171,221,164;102,194,165;50,136,189;94,79,162;1 1 1]);
% cmap = flipud(gray); 
    colormap(cmap)
    imagesc(Spec, 'xdata', Time, 'ydata', F/1000,clims); set(gca, 'ydir', 'normal')
    % parse Tots' analysis file
    SegmentTimes = dbase.SegmentTimes{fi};
    SegmentNames = dbase.SegmentTitles{fi}
    for si = 1:length(SegmentNames)
        if numel(SegmentNames{si}) == 0
            SegmentNames{si} = ''; 
        end
    end
    
    SegmentTimes = SegmentTimes(dbase.SegmentIsSelected{fi}==1,:); 
    SegmentNames(dbase.SegmentIsSelected{fi}==0)=[]; 
    
    Sylls = unique(SegmentNames); 
    sColors = [[.5 .5 .5]; hsv(length(Sylls) - 1)]; 
%     
    hold on
    TopOfSpec = max(ylim);
    for si = 1:size(SegmentTimes,1)
        syllMask = strncmp(SegmentNames{si}, Sylls, 2);
        patch(SegmentTimes(si, [1, 2, 2, 1, 1])/dbase.Fs, TopOfSpec + TopOfSpec*.1 * [0, 0, 1, 1, 0], sColors(syllMask, :), 'Edgecolor','none')
        text(mean(SegmentTimes(si, :))/dbase.Fs,TopOfSpec*1.1, SegmentNames{si}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontsize', 8)
    end
    ylim([0 TopOfSpec*1.1 ])
%     title(num2str(fi)); pause(1)
    
    % zero out part of spec that Tots didn't label
    Spec(:,(Time<SegmentTimes(1)/dbase.Fs)|(Time>SegmentTimes(end)/dbase.Fs)) = [];
    SONG = [SONG Spec]; 
    % do seqNMF
%     RawSong = Spec/max(Spec(:)); 

end

% shg; imagesc(flipud(SONG))
%%
seed = 169%
rng(seed);
figure(seed); 
X = flipud(SONG); 
K = 8; 
L = 40; 
[W, H, cost]= seqNMF(X,'K',K,'L',L,...
     'lambda', .0002,'lambdaL1W', 0, 'lambdaL1H',.1,...  %, 'maxiter', 1, 'showPlot', 0,'W_init', W, 'H_init', H,...
      'SortFactors', 0, 'maxiter', 100);
SimpleWHPlotSpectrograms(W,H,X)
title(seed)
Wold = W; 
Hold = H; 

%%
% for K = 5:10
    seed = 169
savepath = fullfile('C:\Users\emackev\Downloads\tmpPlots\', ['seqNMFDemoSpec' num2str(seed) '_' datestr(now, 'mmddyyyyMMSS'),'.avi']);
obj = vision.VideoFileWriter(savepath, 'AudioInputPort', 0);%,  'fps', 20);
obj.FrameRate = 5; 

X = flipud(SONG); 
K = 8; 
L = 40; 
[N,T] = size(X);
tstart = 410; 
clf; SimpleWHPlotSpectrograms(rand(N,K,L),rand(K,T-tstart+1),X(:,tstart:end),0); 
rng(seed);
set(gcf, 'color', [0 0 0])
saveas(gcf, [savepath(1:end-4) '_pre.pdf']);
cdata = print(gcf, '-RGBImage', '-r300');
[w, h, ~] = size(cdata);
cdata = reshape(cdata, w*h,3);
cdata(mean(cdata,2)==255,:)=1;
cdata = reshape(cdata,w,h,3); 
step(obj, cdata)
[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',L,...
            'lambda', .0002,'lambdaL1W', 0, 'lambdaL1H', .1, ...
            'maxiter', 1, 'showPlot', 0); 

for iteri = 1:1000
    
    if iteri<20 || mod(iteri,10)==0
        clf; SimpleWHPlotSpectrograms(W(:,:,:),H(:,tstart:end), [], 0)
        set(gcf, 'color', [0 0 0])
        title(['Iteration ' num2str(iteri)]); drawnow
        cdata = print(gcf, '-RGBImage', '-r300');
        [w, h, ~] = size(cdata);
        cdata = reshape(cdata, w*h,3);
        cdata(mean(cdata,2)==255,:)=1;
        cdata = reshape(cdata,w,h,3); 
        step(obj, cdata);
    end
    [W, H]= seqNMF(X,'K',K,'L',L,...
            'lambda', .0002,'lambdaL1W', 0, 'lambdaL1H', .1,...
            'maxiter', 1, 'showPlot', 0,...
            'W_init', W, 'H_init', H, 'SortFactors', 0);
end
saveas(gcf, [savepath(1:end-4) '_post.pdf']);
SimpleWHPlotSpectrograms(W(:,:,:),H(:,tstart:end), X(:,tstart:end), 0)
saveas(gcf, [savepath(1:end-4) '_post_rawdata.pdf']);
release(obj)
% end
%%
% 
% % change order of factors to match tots if needed
% 
% W(:,2:3,:) = Wold(:,[3 2],:);
% H(2:3,:) = Hold([3 2],:);
% 
% % plot 
% % spectrogram subplot
% clf
% h(1) = subplot('position', [.1 .5 .6 .35])
% imagesc(Spec, 'xdata', Time, 'ydata', F/1000); set(gca, 'ydir', 'normal')
% hold on
% TopOfSpec = max(ylim);
% sColors = [[.5 .5 .5]; cbrewer('qual', 'Dark2', 3)];%; hsv(length(Sylls) - 1)]; 
% sColors((end+1):(K+1),:) = .5;
% for si = 1:size(SegmentTimes,1)
%     syllMask = strncmp(SegmentNames{si}, Sylls, 2);
%     patch(SegmentTimes(si, [1, 2, 2, 1, 1])/dbase.Fs, TopOfSpec + TopOfSpec*.1 * [0, 0, 1, 1, 0], sColors(syllMask, :), 'Edgecolor','none')
%     text(mean(SegmentTimes(si, :))/dbase.Fs,TopOfSpec*1.1, SegmentNames{si}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontsize', 8)
% end
% ylim([min(ylim) TopOfSpec*1.3])
% set(gca, 'color', 'none', 'tickdir', 'out', 'xtick', [], 'ytick', [])
% cmap = jet; 
% cmap(1,:) = 0; 
% colormap(cmap)
% axis off
% % reconstruction subplot
% h(2) = subplot('position', [.1 .1 .6 .35]);cla
% imagesc(helper.reconstruct(W,H), 'xdata', Time, 'ydata', F/1000); set(gca, 'ydir', 'normal');axis tight
% hold all
% TopOfSpec = max(ylim);
% Kfull = max(find(sum(H,2)>0)); 
% dn = TopOfSpec*.1/Kfull;
% % for ki = 1:K
% %     bump = squeeze(sum(W(:,ki,:),1)); 
% %     bb = min(find(bump>max(bump(:))/2)); 
% %     be = max(find(bump>max(bump(:))/2)); 
% %     bump = 0*bump; 
% %     bump(bb:be)=1;
% %     Bump{ki} = bump; 
% %     Hconv(ki,:) = conv(H(ki,:), [zeros(Lsong,1); bump], 'same')+ ...
% %         conv(H(ki,:), [zeros(Lsong,1); 0*max(squeeze(sum(W(:,ki,:),1)))*ones(Lsong,1)], 'same');
% % end
% Hplot = H/max(H(:)',[],2)*dn*2; 
% for ki = Kfull:-1:1
%     Xs = [1 1:size(H,2) size(H,2)]*.005; % spec is sampled at .005s bins
%     Ys = TopOfSpec+[dn*ki (dn*ki + Hplot(ki,:)) dn*ki]-dn/2;
%     patch(Xs,Ys, sColors(ki+1,:), 'edgecolor', 'none')
% end
% set(gca, 'color', 'none', 'tickdir', 'out', 'ytick', [])
% axis tight
% 
% linkaxes(h, 'x')
% xlim([4 6.3])
% 
% axis off
% %     xlabel('Time(s)');
% drawnow; shg
% % reconstruction subplot
% 
% for ki = 1:K
%         axStart = .71 + (ki-1)*.2/K; 
%         axW = .2/K*.85; 
%         axL(ki) = subplot('Position', [axStart .1 axW .35]);
%         WsToPlot = squeeze(W(:,ki,:));
%         hold on
%         imagesc(WsToPlot, 'ydata', F/1000, [0 max(W(:))])
%     %     title('w')
%         Xs = [1 1:Lsong Lsong];
%         Ys = TopOfSpec+[0 2*[dn dn dn dn] zeros(1,Lsong-4) 0];
%         patch(Xs,Ys, sColors(ki+1, :), 'Edgecolor','none')
%         axis tight; axis off
%         set(gca, 'color', 'none')
% end
% papersize = [5 2]; 
% set(gcf, 'color', 'w', 'papersize', papersize, 'paperposition', [0 0 papersize])
% linkaxes([axL h(2)], 'y')
