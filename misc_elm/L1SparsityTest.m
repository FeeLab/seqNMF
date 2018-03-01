seed = 107; 
rng(seed); 
L = 20;
N = 6; 
T = 3000; 
W = zeros(N,1,L); 
for ni = 0:N-1
    W(ni+1,1,(mod(1:L,N+1)-1)==ni)=1;
end
H = (rand(1,T)<.005);
H(3*L) = 1;
imagesc(squeeze(W));shg
SimpleWHPlot(W,H); 
%%
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
savedir = 'C:\Users\emackev\Dropbox (MIT)\SeqNMF\Figures';
saveTitle = 'L1Example'; 
papersize = [4 2];


[W,H] = seqNMF(helper.reconstruct(W,H)+.0*rand(N,T), 'K', 3, 'L', 2*L,...
    'lambdaL1H', 10, 'lambda', .001, 'maxiter', 150);
SimpleWHPlot(W,H,[],0)
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
saveas(gcf, fullfile(savedir, [saveTitle '_L1H1.pdf'])); 

[W,H] = seqNMF(helper.reconstruct(W,H)+.0*rand(N,T), 'K', 3, 'L', 2*L,...
    'lambdaL1W', 100, 'lambda', .001, 'maxiter', 150);
SimpleWHPlot(W,H,[],0)
set(gcf, 'papersize', papersize, 'paperposition', [0 0 papersize]);%, 'color', 'none')
saveas(gcf, fullfile(savedir, [saveTitle '_L1W1.pdf'])); 

