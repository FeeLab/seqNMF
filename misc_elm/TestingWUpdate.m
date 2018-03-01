seed = 123; 
rng(seed); 
[Batch,~,~] = generate_data(3000,[20,20]',[2,2]',0.0,[0,0]',[1,1]',1,0,0,0,0);



%%
% Batch = generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,0,0,0,0,0);
%
X = X(11:end,:);
X(end-10:3:end,:) = [];
imagesc(X)
X = X([11:end],:);
X(end-10:2:end,:) = [];
nShareds = repmat(1:19,1,3); 
Wcost = zeros(1,length(nShareds)); 
Hcost = zeros(1,length(nShareds)); 
L = 100; K = 4;
for nSharedi = 1:length(nShareds); 
    
    nShared = nShareds(nSharedi);
    nSpecific = 20-nShared;
    Rsh = randperm(20); 
    Rsp = randperm(20);
    [Batch,~,~] = generate_data(6000,[20,20]',[2,2]',0.0,[0,0]',[1,1]',1,0,0,0,0);
    X = Batch([40+sort(Rsh(1:nShared)) 20+sort(Rsp(1:nSpecific))],:);
    
    figure(1)
    [W,H,temp,Loadings,Power] = seqNMF(X, ...     
        'K', K, 'L', L, 'lambda',.0001, ...        
        'showPlot', 0, 'maxiter', 200, 'shift', 1,...
        'lambdaOrthoW', 0,'lambdaOrthoH',1);
    figure(1); SimpleWHPlot(W,H); 
    Wflat = squeeze(sum(W,3)/L); 
    Wcost(nSharedi) = sum(sum((Wflat'*Wflat).*~eye(K))); 
    Hcost(nSharedi) = sum(sum((conv2(H,ones(1,2*L-1), 'same')*H').*~eye(K))); 
    figure(2); clf; hold on
    jitter = .5*(rand(size(nShareds))-.5); 
    scatter(nShareds+jitter, Wcost,'k', 'markerfacecolor', 'flat')
    scatter(nShareds, Hcost/100,'m', 'markerfacecolor', 'flat')
    axis tight
    legend('Wcost', 'Hcost'); drawnow
end
    %%
lambdas = logspace(-5,-2,20); 
Wcost = zeros(1,length(lambdas)); 
Hcost = zeros(1,length(lambdas)); 
RecCost = zeros(1,length(lambdas)); 
seqNMFCost = zeros(1,length(lambdas)); 
K = 5; 
lambda = 0; 
for li = 1:50; %length(lambdas)
%     lambda = lambdas(li);
[W,H,temp,Loadings,Power] = seqNMF(X, ...     
    'K', K, 'L', 100, 'lambda',lambda, ...        
    'showPlot', 0, 'maxiter', 100, 'shift', 1,...
    'lambdaOrthoW', 0,'lambdaOrthoH',0);
figure(1); SimpleWHPlot(W,H); title(num2str(lambda)); shg 
 [RecCost(li),seqNMFCost(li),~] = helper.get_seqNMF_cost(X,W,H);
    Wflat = squeeze(sum(W,3)); 
    Wcost(li) = sum(sum((Wflat'*Wflat).*~eye(K))); 
    Hcost(li) = sum(sum((H*H').*~eye(K))); 
     figure(2); clf; hold on
    subplot(2,1,1);hold on
    scatter(seqNMFCost,Wcost,'k', 'markerfacecolor', 'flat')
    ylabel('Wcost'); xlabel('seqNMFCost')
    subplot(2,1,2);hold on
    scatter(seqNMFCost, Hcost,'k','markerfacecolor', 'flat')
    ylabel('Hcost'); xlabel('seqNMFCost')
%     scatter(lambdas,(RecCost),'b', 'markerfacecolor', 'flat')
%     plot(lambdas,smooth((RecCost),3), 'b'); axis tight
%     title('RecCost')
%     set(gca, 'xscale', 'log')
%     subplot(4,1,2);hold on
%     scatter(lambdas,(seqNMFCost),'r', 'markerfacecolor', 'flat')
%     plot(lambdas,smooth((seqNMFCost),3), 'r'); axis tight
%     set(gca, 'xscale', 'log')
%     title('seqNMFCost')
%     subplot(4,1,3);hold on
%     scatter(lambdas,(Wcost),'m', 'markerfacecolor', 'flat')
%     plot(lambdas,smooth((Wcost),3), 'm'); axis tight
%     set(gca, 'xscale', 'log')
%     title('Wcost')
%     subplot(4,1,4);hold on
%     scatter(lambdas,(Hcost),'g', 'markerfacecolor', 'flat')
%     plot(lambdas,smooth((Hcost),3), 'g'); axis tight
%     set(gca, 'xscale', 'log')
%     title('Hcost')
    drawnow
end
%%
lambdas = logspace(-5,-1,50); 

Wcost = zeros(1,length(lambdas)); 
Hcost = zeros(1,length(lambdas)); 
RecCost = zeros(1,length(lambdas)); 
seqNMFCost = zeros(1,length(lambdas)); 
K = 2; 
L = 100;
for li = 1:length(lambdas)
    lambda = lambdas(li);
    
    [W,H,temp,Loadings,Power] = seqNMF(X, ...     
    'K', K, 'L', L, 'lambda',lambda, ...        
    'showPlot', 0, 'maxiter', 100, 'shift', 1, 'useWupdate', 1);
    figure(1); SimpleWHPlot(W,H);shg

    [RecCost(li),seqNMFCost(li),~] = helper.get_seqNMF_cost(X,W,H);
    Wflat = squeeze(sum(W,3)); 
    Wcost(li) = sum(sum((Wflat'*Wflat).*~eye(K))); 
    Hcost(li) = sum(sum((H*H').*~eye(K))); 
    figure(2); clf; hold on
    subplot(4,1,1);hold on
    scatter(lambdas,(RecCost),'b', 'markerfacecolor', 'flat')
    plot(lambdas,smooth((RecCost),3), 'b'); axis tight
    title('RecCost')
    set(gca, 'xscale', 'log')
    subplot(4,1,2);hold on
    scatter(lambdas,(seqNMFCost),'r', 'markerfacecolor', 'flat')
    plot(lambdas,smooth((seqNMFCost),3), 'r'); axis tight
    set(gca, 'xscale', 'log')
    title('seqNMFCost')
    subplot(4,1,3);hold on
    scatter(lambdas,(Wcost),'m', 'markerfacecolor', 'flat')
    plot(lambdas,smooth((Wcost),3), 'm'); axis tight
    set(gca, 'xscale', 'log')
    title('Wcost')
    subplot(4,1,4);hold on
    scatter(lambdas,(Hcost),'g', 'markerfacecolor', 'flat')
    plot(lambdas,smooth((Hcost),3), 'g'); axis tight
    set(gca, 'xscale', 'log')
    title('Hcost')
    suptitle('Wupdate')
end
%
Wcost = zeros(1,length(lambdas)); 
Hcost = zeros(1,length(lambdas)); 
RecCost = zeros(1,length(lambdas)); 
seqNMFCost = zeros(1,length(lambdas)); 
lambdas = logspace(-5,-1,50); 
K = 2; 
L = 100;
for li = 1:length(lambdas)
    lambda = lambdas(li);
    [W,H,temp,Loadings,Power] = seqNMF(X, ...     
    'K', K, 'L', L, 'lambda',lambda, ...        
    'showPlot', 0, 'maxiter', 100, 'shift', 1, 'useWupdate', 0);
    figure(3); SimpleWHPlot(W,H);shg
    [RecCost(li),seqNMFCost(li),~] = helper.get_seqNMF_cost(X,W,H);
    Wflat = squeeze(sum(W,3)); 
    Wcost(li) = sum(sum((Wflat'*Wflat).*~eye(K))); 
    Hcost(li) = sum(sum((H*H').*~eye(K))); 
    figure(2);
    subplot(4,1,1);hold on
    scatter(lambdas,(RecCost),'cdata', .7*[1 1 1], 'markerfacecolor', 'flat')
    plot(lambdas,smooth((RecCost),3), 'k'); axis tight
    title('RecCost')
    set(gca, 'xscale', 'log')
    subplot(4,1,2);hold on
    scatter(lambdas,(seqNMFCost),'cdata', .7*[1 1 1],'markerfacecolor', 'flat')
    plot(lambdas,smooth((seqNMFCost),3), 'k'); axis tight
    set(gca, 'xscale', 'log')
    title('seqNMFCost')
    subplot(4,1,3);hold on
    scatter(lambdas,(Wcost),'cdata', .7*[1 1 1],'markerfacecolor', 'flat')
    plot(lambdas,smooth((Wcost),3), 'k'); axis tight
    set(gca, 'xscale', 'log')
    title('Wcost')
    subplot(4,1,4);hold on
    scatter(lambdas,(Hcost),'cdata', .7*[1 1 1], 'markerfacecolor', 'flat')
    plot(lambdas,smooth((Hcost),3), 'k')
    set(gca, 'xscale', 'log')
    title('Hcost')
    suptitle('black and gray are no W update')
end
%%
figure(2)

%seed = 1;
%rng(seed)
K = 2
L = 100;
Hcost = []; 
Wcost = []; 
ErrCost = [];
crossCost = []; 


for iteri = 1:50
    [W,H,temp,Loadings,Power] = seqNMF(X, ...     
        'K', K, 'L', L, 'lambda',0.000, ...        
        'showPlot',0, 'maxiter', 100, 'shift', 1, ...
        'lambdaOrthoW', .1);
    figure(3); SimpleWHPlot(W,H);
    Hcost(iteri) = sum(sum((H*H').*~eye(K)));
    Wflat = sum(W,3); 
    Wcost(iteri) = sum(sum((Wflat'*Wflat).*~eye(K)));
    [N,K,L] = size(W);
    [~,T] = size(H);
    [ErrCost(iteri),crossCost(iteri),~] = helper.get_seqNMF_cost(X,W,H);
%     figure(4); 
%     scatter(Hcost,Wcost, 'cdata', [1 0 0]); xlabel('Hcost');ylabel('Wcost')
end
folder = 'C:\Users\emackev\Downloads\tmpPlots\TestingDiffOrtho'; 
% save(fullfile(folder, 'orthWpoint1')); 
%%
folder = 'C:\Users\emackev\Downloads\tmpPlots\TestingDiffOrtho'; 
DIR = dir(fullfile(folder, '*.mat')); 
for di = 1:length(DIR)
    load(fullfile(folder, DIR(di).name))
    figure(4); 
    subplot(2,2,di)
    scatter(Wcost./(crossCost+eps),Hcost./(crossCost+eps), 'cdata', [0 0 0], 'markerfacecolor', 'flat'); 
    xlabel('Wcost/crossCost');ylabel('Hcost/crossCost')
    title(DIR(di).name); 
    xlim([0 .1])
    ylim([0 2.5e-5])
end