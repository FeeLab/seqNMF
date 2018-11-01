function diss = DISSX(H1,W1,H2,W2)
% ------------------------------------------------------------------------
% USAGE:  
% ------------------------------------------------------------------------
% INPUTS
% H1 and W1 are the output of a convNMF fit with K factors
% H2 and W2 are the output of a different convNMF with the same K.
% ------------------------------------------------------------------------
% OUTPUTS
% diss:     Diss is a measure of the disssimilarity of two factorizations.
% ------------------------------------------------------------------------
% Emily Mackevicius and Andrew Bahle
% adapted from Wu et al 2016

    [K,~] = size(H1);
    for i = 1:K
        for j = 1:K
        Xhat1 = helper.reconstruct(W1(:,i,:),H1(i,:));
        Xhat2 = helper.reconstruct(W2(:,j,:),H2(j,:));
        C(i,j) = (Xhat1(:)'*Xhat2(:))/((sqrt(Xhat1(:)'*Xhat1(:))*sqrt(Xhat2(:)'*Xhat2(:)))+eps);
        end
    end
    maxrow = max(C,[],1); 
    maxcol = max(C,[],2); 
    maxrow(isnan(maxrow)) = 0;
    maxcol(isnan(maxcol)) = 0;
    diss = 1/2/K*(2*K - sum(maxrow) -sum(maxcol)); 

end