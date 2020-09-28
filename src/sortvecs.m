function [W2,U2] = sortvecs(W,U)
W2 = W(:,1);
U2 = U(1,:);
W = W(:,2:end);
U = U(2:end,:);
for j = 2:(size(W,2)+1);
    [~,ind] = max(W2(:,j-1).'*W);
    W2 = [W2,W(:,ind)];
    U2 = [U2;U(ind,:)];
    W = [W(:,1:(ind-1)),W(:,(ind+1):end)];
    U = [U(1:(ind-1),:);U((ind+1):end,:)];
end
