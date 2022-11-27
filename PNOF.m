function [Ystar,iter, obj] = PNOF(KH,k,lambda,beta,Y,numclass)

num = size(KH, 1); %the number of samples
numker = size(KH, 3); %m represents the number of kernels
maxIter = 100; %the number of iterations
H = zeros(num,k,numker);
G = zeros(num,num,numker);
Z = zeros(num,num);
opt.disp = 0;
Hstar = zeros(num,k);
F = zeros(num,k,numker);
for p = 1:numker
    P(p) = 1;
end

flag = 1;
iter = 0;

while flag
    iter = iter +1;  
    %% the first step-- optimize H_i
    for p=1:numker  
        Gp = G(:,:,p);
        A = KH(:,:,p)+lambda*(Gp+Gp'-eye(num)-Gp'*Gp);
        [AP,~] = eigs(A, k, 'la', opt);
        H(:,:,p) = AP;
    end
    %%
    %the second step-- optimize S_i
    for p = 1:numker
        K =  H(:,:,p)*H(:,:,p)';
        tmp = (beta*Hstar*F(:,:,p)' + lambda*K)/(lambda*K + beta*eye(num));
        for ii= 1:num
            idx = 1:num;
            idx(ii) = [];
            G(ii,idx,p) = EProjSimplex_new(tmp(ii,idx));
        end
    end
    %%
    %the third step-- optimize Hstar
    T = zeros(num,k);
    for i = 1:numker
        Fp = F(:,:,i);
        Sp = G(:,:,i);
        temp = Sp*Fp;
        pp = P(i)/sum(P);
        T = T+pp*temp;
    end
    Hstar = max(T/numker,0);
    
    %%
    %the fourth step-- optimize Fi
    
    for p = 1:numker
        Sp = G(:,:,p);
        temp = Sp'*Hstar;
        [Ur, ~, Vr] = svds(temp,k);
        F(:,:,p) = Ur*Vr';
    end
    
    %%
    for i = 1:numker
       P(i) = 1/norm(G(:,:,i)-Hstar*F(:,:,i)','fro'); 
    end
    %%
    term1 = 0;
    term2 = 0;
    term3 = 0;
    for j = 1:numker
        term1 = term1+ trace(KH(:,:,j)-KH(:,:,j)*H(:,:,j)*H(:,:,j)');
        term2 = term2+ lambda*norm((H(:,:,j)-G(:,:,j)*H(:,:,j)),'fro')^2;
        term3 = term3 + beta*norm(G(:,:,j)-Hstar*F(:,:,j)','fro')^2;
    end
   
    obj(iter) = term1+term2+term3;
    if (iter>2) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-7 || iter>maxIter)
        flag = 0;
    end
end
[non ,Ystar] = max(Hstar,[],2);
% plot(obj,'r','LineWidth',1.5);
% grid on;