function [max_factor, L_sort, max_sort, hybrid] = ClusterByFactor(W,nclust)
% 
% DESCRIPTION:
% ------------------------------------------------------------------------
% Returns a set outputs which characterize the prefered offset and factors
% and can be further used as indices for plotting etc.
% 
% ------------------------------------------------------------------------
% USAGE: 
% [max_factor, L_sort, max_sort, hybrid] = ClusterByFactor(W,1)
% indSort = hybrid(:,3);
% figure; SimpleWHPlot(W(indSort,:,:),H, helper.reconstruct(W(indSort,:,:),H),0); 
%
% ------------------------------------------------------------------------
% INPUTS
% W:      W is a NxKxL matrix which is the output of convolutional NMF
%         where N is the number of neurons, K is the number of factors in
%         NMF and L is the maximum lag.
%
% ------------------------------------------------------------------------
% OUTPUTS
% max_factor:  max_factor is an cell array of length K containing where
%              each cell contains an NxL matrix which is a logical containing the
%              preferred lag for each neuron for factor K
%                 
%      
% L_sort:      L sort is a KxN matrix where the Kth row is the index which
%              sorts max_factor{k} by lag.
%              USAGE: max_factor{3}(L_sort(3,:),:))
% 
% max_sort:    max_sort is similar to L_sort except the Kth row sorts W(:,K,:)
%              by neurons which are most active in factor K
%
% hybrid:      Hybrid is a Nx3 matrix where the first column is the best
%              lag of a neuron, the second column is the best factor of a
%              neuron and the third column is the neuron index
%
% TODO: Option to sort hybrid by activitation and not lag


[N,K,L] = size(W);
max_factor = cell(1,K);
for ii = 1:K
    data = reshape(W(:,ii,:),[N,L]);
    for jj = 1:size(data,1)
        max_factor{ii}(jj,:) = data(jj,:) == max(data(jj,:));
        if sum(max_factor{ii}(jj,:)) > 1
            max_factor{ii}(jj,:) = 0;
        end
    end
    clust(:,ii) = max(data,[],2);
    [tmp,~] = find(max_factor{ii});
    L_sort(ii,:) = 0;
    L_sort(ii,1:length(tmp)) = tmp;
    [~,xx] = sort(max(data,[],2));
    max_sort(ii,:) = flip(xx);
end

%%
[idx,~] = kmeans(clust,nclust);
N = size(W,1);
pos(:,3) = 1:N;
for ii = 1:N
    nni = reshape(W(ii,:,:),[K,L]);
    try
    [pos(ii,2),pos(ii,1)] = find(nni == max(nni(:),[],1));
    catch
       pos(ii,1:2) = [L,K+1] ;
    end
end

%%
hybrid = pos;
[~,yy] = sort(hybrid(:,2));
hybrid = hybrid(yy,:);
temp = [];
for ii = 1:K+1
    hy = hybrid(hybrid(:,2) == ii,:);
    [~,yy]= sort(hy(:,1),'ascend');
    temp = vertcat(temp,hy(yy,:));
end
hybrid = temp;
hybrid(:,4) = idx(hybrid(:,3));
%% What I was doing before
% hybrid = pos;
% [~,yy] = sort(hybrid(:,2));
% hybrid = hybrid(yy,:);
% hy1 = hybrid(hybrid(:,2) == 1,:);
% hy2 = hybrid(hybrid(:,2) == 2,:);
% [~,yy]= sort(hy1(:,1),'ascend');
% hy1 = hy1(yy,:);
% [~,yy]= sort(hy2(:,1),'ascend');
% hy2 = hy2(yy,:);
% hybrid = vertcat(hy1,hy2);
[~,ii] = sort(hybrid(:,4));
hybrid = hybrid(ii,:);
end