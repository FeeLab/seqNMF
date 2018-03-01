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
for fi = [1:2 9];%length(SoundFiles)
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
xlim([4.38   6.09])
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'SpectrogramFigs'; 
papersize = [5 1];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
if readytosave
%     set(gcf, 'color', 'none')
    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
    saveas(gcf, fullfile(savedir, [saveTitle '_TotsLabeling.pdf'])); 
end
readytosave = 0

% shg; imagesc(flipud(SONG))
%%
%% Procedure for choosing lambda
nLambdas = 30; % increase if you're patient
X = flipud(SONG); 
K = 8; 
L = 40; 
lambdas = sort([logspace(-2,-5,nLambdas)], 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 

for li = 1:length(lambdas)
    
    [W, H]= seqNMF(X,'K',K,'L',L,...
     'lambda', lambdas(li),'lambdaL1W', 0, 'lambdaL1H',0,...  %, 'maxiter', 1, 'showPlot', 0,'W_init', W, 'H_init', H,...
      'SortFactors', 1, 'maxiter', 300, 'showplot', 0);
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    SimpleWHPlot(W,H); shg; drawnow
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end
%% plot costs as a function of lambda
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization); 
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization-minRs)/(maxRs-minRs); 
Cs = filtfilt(b,a,cost); 
minCs =  prctile(cost,10); maxCs =  prctile(cost,90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (cost -minCs)/(maxCs-minCs); 

clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'tickdir','out','ticklength', [0.025, 0.025])
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'Spectrograms'; 
papersize = [4 3];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_lambdas.pdf'])); 

readytosave = 0
%%
for seed = [202]%202 205 % 202, lambda .0003, K8, L40, maxiter1000 works well
    rng(seed);
X = flipud(SONG); 
K = 8; 
L = 40; 
[W, H, cost]= seqNMF(X,'K',K,'L',L,...
     'lambda', .0003,'lambdaL1W', 0, 'lambdaL1H',0,...  %, 'maxiter', 1, 'showPlot', 0,'W_init', W, 'H_init', H,...
      'SortFactors', 1, 'maxiter', 1000, 'showplot', 0);
figure; SimpleWHPlotSpectrograms(W,H(:,760:end),[],0)
title(seed); 
drawnow
end
Wold = W; 
Hold = H; 
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'SpectrogramFigs'; 
papersize = [5 2];
SimpleWHPlotSpectrograms(W,H(:,760:end),[],0)

set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
    set(gcf, 'color', 'none')
    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
%     saveas(gcf, fullfile(savedir, [saveTitle '_seqNMFLabeling.pdf'])); 


SimpleWHPlotSpectrograms(W,H(:,760:end),X(:,760:end),0)
    set(gcf, 'color', 'none')
    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
%     saveas(gcf, fullfile(savedir, [saveTitle '_seqNMFLabeling_raw.pdf'])); 


readytosave = 0
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
