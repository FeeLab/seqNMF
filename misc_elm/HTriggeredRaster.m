function HTriggeredRaster(H,V,L); 
nPks = 5; % number of examples to plot
pkHt = .1; % min peak height compared to maximum
H(sum(H,2)==0,:)=[];
[N T] = size(V);
[K,~] = size(H);

clf
for fi = 1:K
    axStart = .1 + (fi-1)*.8/K; 
    axW = .8/K*.85; 
    axL(fi) = subplot('Position', [axStart .1 axW .8]);

    [pks, locs] = findpeaks(H(fi,:),'MinPeakHeight',max(H(fi,:))*pkHt);
    
    [~, ind] = sort(pks, 'descend'); 
    ind = ind(1:min(nPks,length(ind))); 
    pks = pks(ind); locs = locs(ind);


    R = []; 
    twin = [1:L];
    c = 1;
    for pi = 1:length(pks)
        R(c,:,:) = V(:,min(locs(pi)+twin,size(V,2)));
        c = c+1;
    end
    
    hold on
    for ni = 1:N
        plot([0 L+1], ni*length(pks)*[1 1], 'color', [1 .8 .8], 'linewidth', .1)
    end
    
    Rmat = reshape(R,N*length(pks), length(twin));
    [Xmat,Ymat] = meshgrid(1:size(Rmat,2),1:size(Rmat,1));
    indBurst = find(Rmat>0);
    s = scatter(Xmat(indBurst)+.9*(-.5+rand(size(indBurst))),Ymat(indBurst), ...
       15, 0*[1 1 1], '.', 'markerfacecolor', 'flat', 'MarkerfaceAlpha',5/nPks, 'MarkerfaceAlpha',5/nPks);

    set(gca, 'ydir', 'reverse')
    axis tight; axis off
end

