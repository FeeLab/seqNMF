function score = similarity(W1,H1,W2,H2)
%%
% W1 = W1.*(W1 > 0.001);
% H1 = H1.*(H1 > 0.001);

X1 = {};
X2 = {};

for ii = 1:size(W1,2)
    X1{ii} = helper.reconstruct(W1(:,ii,:),H1(ii,:));
end

for ii = 1:size(W2,2)
    X2{ii} = helper.reconstruct(W2(:,ii,:),H2(ii,:));
end
%%
for ii = 1:size(W1,2)
    for jj = 1:size(W2,2)
        t1 = X1{ii};
        t2 = X2{jj};
        temp = (t1(:)'*t2(:))/((sqrt(t1(:)'*t1(:))*sqrt(t2(:)'*t2(:)))+eps);
        temp(isnan(temp)) = 0;
        S(ii,jj) = temp;%(1,2);
    end
end
% S(isnan(S)) = 0;
% S(S<0) = eps;
%%
temp = S;
num = 0;
for ii = 1:size(W1,2)
    maximum = temp(temp == max(temp(:)));
    maximum = maximum(1);
    if isempty(maximum)
        maximum = 0;
        num = num + maximum;
    else
        num = num + maximum;
        [r,c]= find(temp == max(temp(:)));
        temp(r,:) = 0;
        temp(:,c) = 0;
    end
end
score = num/(size(W1,2)+eps);
% score(isnan(score)) = 0;
end