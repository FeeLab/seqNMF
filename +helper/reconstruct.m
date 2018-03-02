function X_hat = reconstruct(W,H)
% ------------------------------------------------------------------------
% USAGE: X_hat = helper.reconstruct(W,H)
% ------------------------------------------------------------------------
% INPUTS
% W:      W is a NxKxL tensor which gives the neuron basis
%         functions which are used for the reconstructions. The L'th NxK slice
%         of W is the neural basis set for a lag of L.
%
% H:      H is a KxT matrix which gives timecourses for each factor 
% ------------------------------------------------------------------------
% OUTPUTS
% X_hat:  The reconstruction X_hat = W (*) H; 
% ------------------------------------------------------------------------
% Emily Mackevicius and Andrew Bahle

[N,K,L] = size(W);
[~,T] = size(H);

% zeropad by L
H = [zeros(K,L),H,zeros(K,L)];
T = T+2*L;
X_hat = zeros(N,T);

for tau = 1:L % go through every offset from 1:tau
    %X_hat = X_hat + W(:, :, tau) * circshift(H,tau-1,2);
    X_hat = X_hat + W(:, :, tau) * circshift(H,[0,tau-1]);
end

% undo zer0padding
X_hat = X_hat(:,(L+1):(end-L)); 

end