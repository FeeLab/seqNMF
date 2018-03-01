function [data,W,H,V_hat] = generate_data(T,Nneurons,Dt,NeuronNoise,SeqNoiseTime,SeqNoiseNeuron,shared,diff,stretch,bin,seed)
%rng(2001)
if seed == 0
    rng shuffle
else
    rng(seed)
end
%stretch = 0;
additional_neurons = 0;
non_sparse = 1;
if shared
    Nneurons = [Nneurons;Nneurons(1)];
    Dt = [Dt;Dt(1)];
    SeqNoiseTime = [SeqNoiseTime;SeqNoiseTime(1)];
    SeqNoiseNeuron = [SeqNoiseNeuron;SeqNoiseNeuron(1)];
end

%% Parameters
% V = data.sequences(10000,randi(10,2,1)+5,randi(4,2,1),.01,rand(2),ones(2,1)*0.95);
% Nneurons = [5,15,10,8]; % the number of neurons in each sequence
% Dt = [2,1,3,3]; % the number of time steps between each neuron in the sequence
% Pseq = []; % the probability of the sequence occuring
% NeuronNoise = 0.01; % the noise in a neurons firing rate
% SeqNoiseTime = [0.2,0.2,0.1,0.1]; % the noise in the sequence aka jitter (p of each neuron jittered 1 dt)
% SeqNoiseNeuron = [0.95,0.95,0.95,0.95]; % the probability that a neuron participates in a given seq
% T = 1000;
% Share = []; % the propotion of the chain that is shared in other sequences
%% Calculate useful things
N = sum(Nneurons)+additional_neurons; % Total number of neurons
nseq = length(Nneurons); % The number of sequences
lseq = Dt.*Nneurons; % the length of each sequences

j = 1;
neurons = {};
for ii = 1:length(Nneurons)
    neurons{ii} = j:j+Nneurons(ii)-1;
    j = j+Nneurons(ii);
end

%% MAKE H's
xx = zeros(1,T);
H = zeros(nseq,T);
%randomly distribute seq starting points preventing the same sequence from initiation during itself
% for ii = 1:length(lseq)
%     pos(ii,:) = cumsum(randi(lseq(ii)+50,50,1)+lseq(ii)); %sp sets how often the seq happen
%     temp = pos(ii,:);
%     H(ii,temp(temp<T))= 1;
% %     H(ii,logical(xx)) = 0;
% %     for jj = 1:length(pos(ii,:))
% %         xx(pos(ii,jj):pos(ii,jj)+lseq(ii)) = 1;
% %     end
%       
% end
nn = nseq*1000; % make smaller
if stretch > 0
    stretches = randi(stretch,nn,1);
else
    stretches = zeros(nn,1);
end
%temp = cumsum(randi(450,nn,1)+max(lseq)+stretch); = 1:nn
temp = [];
for j = 1:nn
    % need to fix this to work for stretches = 0
    %temp = [temp;randi(100,1,1)+max(lseq)+stretches(ii)];
    temp = [temp;randi(100,1,1)+(max(lseq)/max(Dt)*(max(Dt)+stretches(ii)))];
end
temp = cumsum(temp);

%indx = randi(nseq,nn,1);
indx = ones(nn/nseq,1)';
for ii = 2:nseq
    indx = [indx,ii*ones(nn/nseq,1)'];
    %indx = [ones(nn/2,1);2*ones(nn/2,1)]';
end
indx = indx(randperm(nn));
for ii = 1:nseq
    H(ii,temp((indx == ii))) = 1;
    Hs{ii} = stretches((indx == ii));
end
if shared
     H(end,:) = sum(H(1:end-1,:));
end
H = H(:,1:T);


%% Make Data using noise parameters in the reconstruction
W = zeros(N,nseq,max(lseq)+150);
[N,K,L] = size(W);
if shared && diff
    index(1,:) = 1:Nneurons(1);
    for ii = 2:K-1
        index(ii,:) = randperm(Nneurons(1)); 
    end
end

%leng = max(Dt)+ max(lseq) + (stretch)*max(lseq); 
leng = (max(lseq)/max(Dt)*(max(Dt)+stretch));
L = leng+150;
H(:,T-(2*(max(lseq)/max(Dt)*(max(Dt)+(stretch)))):T) = 0;
%H(:,end-300:end) = 0;
[~,T] = size(H);
V_hat = zeros(N,T);%+L-1);

%Dont forget these things!
%NeuronNoise = [0.1,0.05,0.2,0.13]; % the noise in a neurons firing rate
%SeqNoiseTime = [0.1,0.2,0.1,0.1]; % the noise in the sequence aka jitter (p of each neuron jittered 1 dt)
%SeqNoiseNeuron = [0.9,0.9,0.9,0.9]; % the probability that a neuron participates in a given seq
if shared 
    Ktemp = K-1;
else 
    Ktemp = K;
end
for ii = 1:Ktemp % go through each factor
    ind = find(H(ii,:));
    for jj = 1:sum(H(ii,:)) % go through each iteration of the sequence
        tempH = zeros(1,size(H,2));
        tempH(ind(jj)) = 1;
                   
        if stretch > 0 % change the dt for each instance
            Dt_temp = Dt(ii)+Hs{ii}(jj);%+(randi(stretch))
            %*(-1+(2*(rand(1)>0.5))); If you want compression as well
            
            tempW = zeros(N,leng+150);          
            temp = eye(length(neurons{ii}));
            %rng(ii)
            %temp = rand(size(temp))>0.7;
            if size(temp,2) < leng
                temp2 = zeros(length(neurons{ii}),Dt_temp*Nneurons(ii));
                temp2(:,1:Dt_temp:Dt_temp*Nneurons(ii)) = temp;    
                tempW(neurons{ii},50:49+size(temp2,2)) = temp2;  
            else
                tempW(neurons{ii},50:lseq(ii)+49) = temp; 
            end
        else 
            Dt_temp = Dt(ii);
            
            tempW = zeros(N,leng+150);          
            temp = eye(length(neurons{ii}));
            %rng(ii)
            %temp = rand(size(temp))>0.7;
            if size(temp,2) < leng
                temp2 = zeros(length(neurons{ii}),Dt_temp*Nneurons(ii));
                temp2(:,1:Dt_temp:Dt_temp*Nneurons(ii)) = temp;    
                tempW(neurons{ii},50:49+size(temp2,2)) = temp2;  
                if shared
                     temp = eye(length(neurons{end}));
                     if diff
                        temp = temp(index(ii,:),:);
                     end
                     temp2 = zeros(length(neurons{ii}),Dt_temp*Nneurons(ii));
                     temp2(:,1:Dt_temp:Dt_temp*Nneurons(ii)) = temp;
                     tempW(neurons{end},50:49+size(temp2,2)) = temp2; 
                end
            else
                tempW(neurons{ii},50:lseq(ii)+49) = temp; 
            end      
        end

        % neurons are jittered with some lambda
        %shifts = poissrnd(SeqNoiseTime(ii),N,1).*(1-2*(rand(N,1)>0.5));
        shifts = round(normrnd(0,SeqNoiseTime(ii),N,1));
        %shifts(abs(shifts) >5) = 0; %stop poiss from getting to big
        for idx = 1:size(W,1)
            tempW(idx,:) = circshift(tempW(idx,:),shifts(idx)*1);
        end
        
        % neurons participate with some p
        tempW(rand(N,1)>SeqNoiseNeuron(ii),:) = 0; 
        [tempW, tempH] = helper.shiftFactors(tempW, tempH);
        shift = circshift(1:length(tempH),floor((Dt(ii)*Nneurons(ii))- Dt_temp*Nneurons(ii))/2);
        tempH = tempH(:,shift);
        newData = conv2(tempH,tempW);
        %V_hat = V_hat + newData(:,ceil(size(tempW,2)/2):end - floor(size(tempW,2)/2));
        V_hat = V_hat + newData(:,1:(end - size(tempW,2)+1));

    end
end


%% could add indepent noise later
V_hat = V_hat + (rand(size(V_hat))<NeuronNoise);
%V_hat  = V_hat./V_hat;
V_hat(isnan(V_hat)) = 0;
%%
%data = V_hat;
if ~bin
    filtbio = [zeros(1,10*10) exp(-(1:10*10)/10)]; 
    data = conv2(V_hat,filtbio,'same');
else
    data = V_hat;
end

W = [];
for ii = 1:K % go through each factor               
        Dt_temp = Dt(ii);
        tempW = zeros(N,leng+150);          
        temp = eye(length(neurons{ii}));  
        temp2 = zeros(length(neurons{ii}),Dt_temp*Nneurons(ii));
        temp2(:,1:Dt_temp:Dt_temp*Nneurons(ii)) = temp;    
        tempW(neurons{ii},50:49+size(temp2,2)) = temp2;    
        W(:,ii,:) = tempW;
end
[W, H] = helper.shiftFactors(W, H);

if ~bin
    for ii = 1:size(W,2)
        W(:,ii,:) = conv2(squeeze(W(:,ii,:)),filtbio,'same');
    end
end
rng shuffle
end
