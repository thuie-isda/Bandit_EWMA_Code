clc;
clear;

n =10; % Image size n*n
T = 2000; % The Total time length
bd = 5; % How smooth the defect is (Order of Spline)
k0 = 20; % scale of the background
bd0 = 3;
nq0 = 2;
k1  = 2; % scale of the defect
ndefect = 1;  % How many defect in the system
%%delta =0; % Defect magnitude
sigma = 0.05; % The noise level
Tin = 1; % When the change occur
Y = zeros(n,n,T); % Generated Videos

repNum=200;
tau=0;
m =7; % budget of pixels that can be observed for each image
lambda=0.1;
limit=11.9;
m0=0;

loss = @(t) sqrt(max(0,log((1-(1-lambda).^t)./lambda)));
Test = zeros(repNum,T-m0);
count1 = zeros(repNum,n);
meanU = count1;

knots0 = [ones(1,bd) linspace(1,n,round(n/k0)) n * ones(1,bd0)];
nKnots0 = length(knots0) - (bd0-1);
kspline0 = spmak(knots0,eye(nKnots0));
kd = fnder(kspline0);
B0 = fourierbasis(nq0,n);
B00 = B0;
sigma0 = 0.1;
Sigma0 = eye(size(B00,2),size(B00,2))*sigma0;
covMat = eye(n)*sigma^2+sigma0.*B00*eye(size(B00,2))*B00'; 

    
indM = combnk(1:n,m);
[nM,~] = size(indM);
invSigM = cell(1,nM);

for i = 1:nM
invSigM{i} = inv(covMat(indM(i,:),indM(i,:)));
end

deltaset=0.1:0.1:1;
ARL_oc=zeros(1,size(deltaset,2));
sdrl_oc=zeros(1,size(deltaset,2));
parpool('local',6);

for dd=1:size(deltaset,2)
    delta=deltaset(dd);
parfor ss=1:repNum
    %generate the background
    knots0 = [ones(1,bd) linspace(1,n,round(n/k0)) n * ones(1,bd0)];
    nKnots0 = length(knots0) - (bd0-1);
    kspline0 = spmak(knots0,eye(nKnots0));
    kd = fnder(kspline0);
    B0=spval(kspline0,1:n)';
    runlen = n;
%    B0 = B0(:,4:end-2);
 %   B0 = B0(:,5:end-2);
     B0 = fourierbasis(nq0,n);
    B00 = B0;
   % B00 = kron(B0,B0);
    mu0 = zeros(1,size(B00,2));
    sigma0 = 0.1;
    Sigma0 = eye(size(B00,2),size(B00,2))*sigma0;
    theta0 = zeros(T,size(B00,2));
    % sample IC data
    Y = [];
    for i = 1:T
    theta0(i,:) = mvnrnd(mu0,Sigma0,1);
    Y(:,i) = B00*theta0(i,:)'+ sigma*randn(n,1);
    end

    %Generate Defect by Spline Bases (We can change this later if we need to generate other shape of defect
    knots = [ones(1,bd) linspace(1,n,round(n/k1)) n * ones(1,bd)];
    nKnots = length(knots) - (bd-1);
    kspline = spmak(knots,eye(nKnots));
    kd = fnder(kspline);
    B=spval(kspline,1:n)';
    B = B(:,3:end-2);
    kcoef = size(B,2);
    Z = zeros(kcoef,1);


    for i = 1:ndefect
        Z(randsample(1:kcoef,1)) = 1;
    end
    %imagesc(Z)

    D = B*Z;
    %imagesc(D)
    %colorbar;

    for i = (Tin+1):T
        Y(:,i) = Y(:,i)+delta*D;
    end
 %   Y1 = reshape(Y,n*n,T);
     Y1 = Y;
%     B = kron(B,B);
 % length of markov chain
    nq = size(B,2); 
    w=0.1*ones(nq,1);
    sigma9=3*ones(nq,1);
    dd
    ss
    %covMat = cov(Y1');

[count1,Q,meanU,runlen,Test] = CMAB(Y1',m0,m,covMat,loss,lambda,limit,indM,invSigM,tau);
TTest(ss,:) = Test;
RRunlength(ss) = runlen;
end
ARL_oc(dd)= mean(RRunlength);
sdrl_oc(dd)=std(RRunlength);
end
delete(gcp('nocreate'))
save('GenerateVector_CMAB.mat')


% xlswrite('..\..\result_update.xlsx',ARL_oc','Sheet3','F53:F72');
% xlswrite('..\..\result_update.xlsx',sdrl_oc','Sheet3','G53:G72');
% limit =49.5
% for rep = 1:repNum
%     ind =  find(Test(rep,2:end)>limit,1);
%     if isempty(ind)
%         runlen1(rep) = T;
%     else
%         runlen1(rep) = ind;
%     end
% end
% mean(runlen1)
% std(runlen1)
% 
% 
% 
% 
% parfor rep =  1:repNum
%     %generate the background
% 
%     rep
%      [count1,Q,meanU,runlen(rep),Test(rep,:)] = MABind(Y1',m0,m,covMat,loss,lambda,limit,tau);
%      plot(Test(rep,:))
% end
% arl_mab = mean(runlen)
% stdrl_mab = std(runlen)

% for rep = 1:repNum
%     ind =  find(Test(rep,tau+1:end)>limit,1);
%     if isempty(ind)
%         runlen1(rep) = k0-m0-tau;
%     else
%         runlen1(rep) = ind;
%     end
% end
% mean(runlen1)
% std(runlen1)
% runlen(find(runlen<=tau)) = NaN;
% nanmean(runlen,2)
% 
% limit = 12.15
% for rep = 1:repNum
%     ind =  find(Test(rep,:)>limit,1);
%     if isempty(ind)
%         runlen1(rep) = k0-m0;
%     else
%         runlen1(rep) = ind;
%     end
% end
% mean(runlen1)
% std(runlen1)
% 
% for i = 1:p
% figure('position',[40,100,1500,1000]);
%     plot(WU(:,i),'r');
% hold on;
% plot(WL(:,i),'k');
% plot(W0(:,i));
% end
% 
% % power
% ll = distparams1{1}*inv(covMat)*distparams1{1}';
% ncx = ncx2cdf(limit,p,ll);
% OCARL = 1/(1-ncx);