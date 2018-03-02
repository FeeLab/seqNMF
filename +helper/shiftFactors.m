function [W,H] = shiftFactors(W,H)
    % shift factors by center of mass
    
    % get size of W and H
    [N,K,L] = size(W);
    [~,T] = size(H);
    
    if L>1; % if L=1, no room to shift
    
    center = max(floor(L/2),1); % middle bin
    
    % pad with zeros, for use with circshift. data is zeropadded within seqNMF, so don't need to pad H
    Wpad = cat(3,zeros(N,K,L),W,zeros(N,K,L));
    
    for k = 1:K
        % compute center of mass
        temp = sum(squeeze(W(:,k,:)),1);
        cmass = max(floor(sum(temp.*(1:length(temp)))/sum(temp)),1);          
        %Wpad(:,k,:) = circshift(squeeze(Wpad(:,k,:)),center-cmass,2);
        %Changing for compatibility with 2016, this is the same
        Wpad(:,k,:) = circshift(squeeze(Wpad(:,k,:)),[0,center-cmass]); 
        %H(k,:) = circshift(H(k,:),cmass-center); 
        H(k,:) = circshift(H(k,:),[0,cmass-center]); 
    end
 

    % undo zero pad
    W = Wpad(:,:,(L+1):(end-L));

    end