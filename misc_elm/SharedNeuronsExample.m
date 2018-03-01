
seed = 1235; 
rng(seed); 

[Batch,~,~] = generate_data(5000,[10,10]',[10,10]',0.0,[0,0]',[1,1]',0,1,0,0,0);
% generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,shared,diff,stretch,bin,seed)

L = 250; K = 4;
figure(1)
[W,H,temp,Loadings,Power] = seqNMF(Batch, ...     
    'K', K, 'L', L, 'lambda',.0005, ...        
    'showPlot', 1, 'maxiter', 100, 'shift', 1,...
    'lambdaOrthoW', 0,'lambdaOrthoH',0);
SimpleWHPlot(W,H,[],0)
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'SharedNeurons'; 
papersize = [8 6];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_noshared.pdf'])); 

[W,H,temp,Loadings,Power] = seqNMF(Batch, ...     
    'K', K, 'L', L, 'lambda',.00000, ...        
    'showPlot', 1, 'maxiter', 100, 'shift', 1,...
    'lambdaOrthoW', 0,'lambdaOrthoH',0);
SimpleWHPlot(W,H,[],0)
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'SharedNeurons'; 
papersize = [8 6];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_cnmf.pdf'])); 

[Batch,~,~] = generate_data(5000,[10,10]',[10,10]',0.0,[0,0]',[1,1]',1,1,0,0,0);
[W,H,temp,Loadings,Power] = seqNMF(Batch, ...     
    'K', K, 'L', L, 'lambda',.0005, ...        
    'showPlot', 1, 'maxiter', 100, 'shift', 1,...
    'lambdaOrthoW', 0,'lambdaOrthoH',0);
SimpleWHPlot(W,H,[],0)
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'SharedNeurons'; 
papersize = [8 6];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_sharedDiff.pdf'])); 
    
[Batch,~,~] = generate_data(5000,[10,10]',[10,10]',0.0,[0,0]',[1,1]',1,0,0,0,0);
[W,H,temp,Loadings,Power] = seqNMF(Batch, ...     
    'K', K, 'L', L, 'lambda',.0005, ...        
    'showPlot', 1, 'maxiter', 100, 'shift', 1,...
    'lambdaOrthoW', 0,'lambdaOrthoH',0);
SimpleWHPlot(W,H,[],0)
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'SharedNeurons'; 
papersize = [8 6];
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
    saveas(gcf, fullfile(savedir, [saveTitle '_sharedSame.pdf'])); 
%% Example parts-based vs events-based factorizations (W or H soft orthogonality)

% change these parameters to switch between parts-based and events-based
lambdaOrthoH = 0; % favor events-based (these can take any value, don't need to be zero and one)
lambdaOrthoW = 1; % favor parts-based

% rng(236); % fixed rng seed for reproduceability
X = trainNEURAL;
K = 3;
L = 2/3; % units of seconds
Lneural = ceil(L*VIDEOfs);
Lsong = ceil(L*SONGfs);
shg
display('For parts-based, set lambdaOrthoW=1 and lambdaOrthoH=0; for events-based, vice versa')
[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambda', .00, 'maxiter', 100, 'showPlot', 1,...
            'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW); 
        
p = .05; % desired p value for factors

display('Testing significance of factors on held-out data')
[pvals,is_significant] = test_significance(testNEURAL,W,p);

W = W(:,is_significant,:); 
H = H(is_significant,:); 

% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 180; % plot data starting at this timebin
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end),0)
title(['lambdaOrthoW=' num2str(lambdaOrthoW) ', lambdaOrthoH=' num2str(lambdaOrthoH)])


% savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
% saveTitle = 'PartsBasedOrthW1'; 
% papersize = [8 6];
% set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
%     saveas(gcf, fullfile(savedir, [saveTitle '_raw.pdf'])); 



% figure; WHPlot(W(indSort,:,:),H(:,tstart:end), ...
%     helper.reconstruct(W(indSort,:,:),H(:,tstart:end)),...
%     0,trainSONG(:,floor(tstart*SONGfs/VIDEOfs):end),0)
% title('SeqNMF reconstruction')

% set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
%     saveas(gcf, fullfile(savedir, [saveTitle '_reconstruction.pdf'])); 