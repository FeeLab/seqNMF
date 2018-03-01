function HTriggeredSpec(H,Spec,Hfs,Specfs,Lspec); 
%%
H(sum(H,2)==0,:)=[];
nPks = 5; % number of examples to plot

[K,~] = size(H);
clf
for fi = 1:K
    axStart = .1 + (fi-1)*.8/K; 
    axW = .8/K*.85; 
    axL(fi) = subplot('Position', [axStart .1 axW .8]);

    
    [m,indH] =  sort(H(fi,:), 'descend');
    plotS = [];
    for ei = 1:nPks
        tSample = indH(ei); 
        if (round(tSample*Specfs/Hfs) +  Lspec)<=size(Spec,2)
            plotS = [Spec(:,round(tSample*Specfs/Hfs):(round(tSample*Specfs/Hfs) +  Lspec)); plotS]; 
        end
    end
    imagesc(plotS)
    cmap  = 1/256*flipud([158,1,66;213,62,79;244,109,67;253,174,97;254,224,139;255,255,191;230,245,152;171,221,164;102,194,165;50,136,189;94,79,162;1 1 1]);
    cmap(1,:) = 0;
    colormap(cmap)
     set(gca, 'ydir', 'normal')
     hold on; plot(repmat([0 Lspec+2],nPks+1,1)',.5+ repmat(0:nPks,2,1)*size(Spec,1),'r')

    axis tight; axis off
end

