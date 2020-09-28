% [H,S] = fact_bayes_resub(Y,b,nu,p,L,maxit,delta,showme) returns a
% convolutive NMF factorization of the square spectrogram Y into H and S,
% where H is associated to the RIR and S is the dereverberated spectrogram.
% 
% if length(B) = 1, then b_k = B*sqrt(sigma2_k).
% if lenght(B) = 2, then b_k is considered to be a random variable with a
% gamma distribution with shape = B(1) and scale = B(2); 
%
% if length(NU) = 1, then nu_k = NU.
% if lenght(NU) = 2, then nu_k is considered to be a random variable with a
% gamma distribution with shape = NU(1) and rate = NU(2); 

function [W,U,H,S] = fact_beta(Y,p)

%% Initializing
ep = 1e-10; %------------------ To avoid dividing by 0
[K,N] = size(Y);    
S = Y; %----------------------- Y is used as a first estimation of Sç

L = spdiags(ones(p.M,1)*[1 -1],[0,1],p.M-1,p.M);
sigma2 = diag(Y*Y');
W = random('exp',1,K,p.J);
W = W*diag(1./sum(W,1));
% save borrar W
U = random('exp',mean(mean(Y)),p.J,N);

%% STAGE I
beta = p.beta1;
gamma = 1/((2-beta)*(beta<1)+(beta>1)*(beta<=2)+(beta-1)*(beta>2));
H = [ones(K,1),zeros(K,p.M-1)];
W0 = 1000*ones(size(W));
for iter = 1:p.it1
    
    % Loop break ----------------------------------------------------------
    if norm(W-W0)<K*1e-6;
        break
    else
        W0 = W;
    end
    % ---------------------------------------------------------------------
        
    % W updating ----------------------------------------------------------
    % update the representation X
    X = W*U;
    % update the auxiliary matrices for diagonalized convolution
    A1 = Y./max(X.^(2-beta),ep);
    if beta > 1
        A2 = X.^(beta-1);
    else
        A2 = 1./max(X.^(1-beta),ep);
    end;
    % compute full bot and top
    top = (A1*U.').^gamma;
    bot = (A2*U.').^gamma;
    % update W
    W = W.*max(top,ep)./(bot+ep);
    if or(iter == 20, iter == 50)
        W = W + random('exp',1,K,p.J);
    end
    W = W*diag(1./sum(W,1));
    % ---------------------------------------------------------------------
        
    % U updating ----------------------------------------------------------
    % update the auxiliary matrices for diagonalized convolution
    A1 = Y./max(X.^(2-beta),ep);
    if beta > 1
        A2 = X.^(beta-1);
    else
        A2 = 1./max(X.^(1-beta),ep);
    end
    % compute full bot and top
    top = (W.'*A1).^gamma;
    bot = (W.'*A2).^gamma;
    % Update U
    U = U.*(max(top,ep)./(bot+ep));
    % ---------------------------------------------------------------------
    
end

%% Stage II

H = exp(-linspace(0,1,p.M));     
H = ones(K,1)*(H/norm(H,1)); 
penalizer = p.ldaU*mean(mean(Y))*(1+10./(ones(p.J,p.J)*U));
beta = p.beta2;
gamma = 1/((2-beta)*(beta<1)+(beta>1)*(beta<=2)+(beta-1)*(beta>2));

for iter = 1:p.it2
    Sp = S;
    
    % X updating ----------------------------------------------------------
    for k = 1:K;
        X(k,:) = filter(H(k,:),1,S(k,:));
    end
    % update the auxiliary matrices for diagonalized convolution
    A1 = Y./max(X.^(2-beta),ep);
    if beta > 1
        A2 = X.^(beta-1);
    else
        A2 = 1./max(X.^(1-beta),ep);
    end
    % ---------------------------------------------------------------------
        
    % S updating ----------------------------------------------------------
    % Discrete convolution
    top = zeros(K,N);
    bot = zeros(K,N);
    for m = 1:p.M;
        top = top + spdiags(H(:,m),0,K,K)*[A1(:,m:end),zeros(K,m-1)];
        bot = bot + spdiags(H(:,m),0,K,K)*[A2(:,m:end),zeros(K,m-1)];
    end
    % Multiplicative update
    top = W.'*top-penalizer;
    bot = W.'*bot;
    U = U.*(max(top.^gamma,ep)./(bot.^gamma+ep));
    S = W*U;
    % ---------------------------------------------------------------------
      
    % H updating ----------------------------------------------------------
    % Auxiliaries
    A1 = Y./max(X.^((2-beta)*(beta<2)),ep);
    A2 = X.^((beta-1)*(beta>1));
    % Discrete convolution
    top = zeros(K,p.M);
    bot = zeros(K,p.M);
    for t = 1:N;
        if (t+p.M-1)<= N
            top = top + spdiags(S(:,t),0,K,K)*A1(:,t:t+p.M-1);
            bot = bot + spdiags(S(:,t),0,K,K)*A2(:,t:t+p.M-1);
        else
            top = top + spdiags(S(:,t),0,K,K)*[A1(:,t:end),zeros(K,t-N+p.M-1)];
            bot = bot + spdiags(S(:,t),0,K,K)*[A2(:,t:end),zeros(K,t-N+p.M-1)];
        end
    end    
    for k = 1:K
        Hk = H(k,:).';
        uno = spdiags(bot(k,:).',0,p.M,p.M);
        dos = 2*p.ldaH*sigma2(k)*spdiags(Hk,0,p.M,p.M)*(L'*L);
        tres = spdiags(Hk,0,p.M,p.M)*top(k,:).';
        Hk = (uno + dos) \ tres;
        Hk = Hk/max(Hk);
        H(k,:) = Hk.';
    end
    % ---------------------------------------------------------------------
        
    % Stopping criterion --------------------------------------------------
    if norm(Sp-S)<p.delta
        fprintf('1%d iterations',iter);
        break
    end
    % ---------------------------------------------------------------------
end

% Gain function -----------------------------------------------------------
S = W*U;
for k = 1:K;
    X(k,:) = filter(H(k,:),1,S(k,:));
end
G = S./(X+ep);
S = G.*Y;
% -------------------------------------------------------------------------

end


