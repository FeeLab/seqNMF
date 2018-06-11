function [pvals,is_significant] = test_significance(TestData,W,p,nnull)
%
% USAGE: 
%
% [pvals,is_significant] = test_significance(TestData,W,0.01)
%
% ------------------------------------------------------------------------
% DESCRIPTION:
%
% Tests each factor in W for significance using a held out test dataset at
% a p-value of p using Bonferroni correction. 
%  
% ------------------------------------------------------------------------
%
% INPUTS:
%
% Name              Default                 Description
% TestData                                  Held out data matrix (NxT) 
% W                                         NxKxL tensor containing factor exemplars
% p                 0.05                    Desired p-value to test
% nnull             ceil(K/p)*2             Number of null datasets to use
%
% ------------------------------------------------------------------------
% OUTPUTS:
%
% pvals                      A vector (1xK) containing the p-value of each factor
% is_significant             A boolean vector (1xK) which is 1 if a factor
%                            is significant using the specified pvalue with correction
%
% ------------------------------------------------------------------------
% CREDITS:
%   Emily Mackevicius and Andrew Bahle, 2/1/2018
%
%   Please cite our paper: 
%       XXXXXXXXXXXXXXXXXXXXX
% Remove factors where there is obviously no sequence
% That is, W is empty, or has one neuron with >99.9% of the power

indempty = sum(sum(W>0,1),3)==0; % W is literally empty
Wflat = sum(W,3); 
indempty = indempty | (max(Wflat,[],1).^2> .999*sum(Wflat.^2,1)); % or one neuron has >99.9% of the power
W(:,indempty,:) = []; % Delete factors that meet the above critera

[N,K,L] = size(W);
[~,T] = size(TestData);

if nargin < 3
    p = 0.05;
end

if nargin < 4
    nnull = ceil(K/p)*2;
end

% make nnull shifted datasets 
skewnull = zeros(K,nnull);

X = TestData;

for n = 1:nnull
    % Make a null dataset
    Wnull = zeros(N,K,L);
    for k = 1:K
        for ni = 1:N
            %Wnull(ni,k,:) = circshift(W(ni,k,:),randi(L));
            Wnull(ni,k,:) = circshift(W(ni,k,:),[0,0,randi(L)]);
        end
    end
    %figure, imagesc(squeeze(Wnull(:,1,:))), colormap(gca, flip(gray)),axis off
    % Calculate WTX
    WTX = zeros(K, T);
    for l = 1 : L
        %X_shifted = circshift(X,-l+1,2);       
        X_shifted = circshift(X,[0,-l+1]);       
        WTX = WTX + Wnull(:, :, l)' * X_shifted;
    end   
    %figure, histogram(WTX(1,:),0:1:100,'FaceColor','k')
    %set(gca,'yscale','log'), ylim([0,10e2])
    % Get skewness of each
    skewnull(:,n) = skewness(WTX,1,2);
end


WTX = zeros(K, T);
for l = 1 : L
    %X_shifted = circshift(X,-l+1,2);       
    X_shifted = circshift(X,[0,-l+1]);       
    WTX = WTX + W(:, :, l)' * X_shifted;
end   
skew = skewness(WTX,1,2);
for k = 1:K
    % Assign pvals from skewness
    pvals(k) = (1+sum(skewnull(k,:)>skew(k)))/nnull;
end
allpvals(indempty) = Inf; 
allpvals(~indempty) = pvals; 
pvals = allpvals;
is_significant = (pvals <= p/K);
end