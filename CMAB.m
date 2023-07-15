function [count1,Q,meanU1,runlen,Test,regret] = CMAB(data,m0,m,covMat,loss,lambda,limit,indM,invSigM,tau)
[kmax,dim] = size(data);
%meanU0 = mean(data(1:m0,:));
meanU0=zeros(1,dim);
%covU0 = cov(data(1:m0,:));
invcovMat = inv(covMat);
%inv5covMat = covMat^-0.5;
%loss0 = shiftmean*invcovMat*shiftmean';
%regret = zeros(1,kmax-m0);
%data = data*inv5covMat;
%Q = zeros(kmax-m0,size(indM,1));
Q = zeros(1,size(indM,1));
% W0 = zeros(kmax-m0,dim);
% WL = W0;
% WU = W0;
count = zeros(kmax-m0,dim);
count(1,:) = 1;
runlen = kmax-m0-tau;
%indicate = zeros(kmax-m0,dim);
index = zeros(m,kmax-m0);
meanU = zeros(kmax-m0,dim);
sumU= zeros(kmax-m0,dim);
Omega = cell(kmax-m0,1);
ww = 1-lambda;
%ww1 = zeros(1,kmax-m0);
Test = zeros(1,kmax-m0);
for i = 1:kmax-m0
    if i == 1
       Omega{i} = invcovMat+inv(covMat);
       sumU(i,:) = meanU0 + data(m0+1,:);
       meanU(i,:) = sumU(i,:)*inv(Omega{i});
       ww(i) = 1;
    end
    if i >1
         ww(1:i-1) =ww(1:i-1).*(1-lambda);
         ww(i) = 1;
        tmp = randsample(size(indM,1),size(indM,1),false);
        tempQ = Q(i-1,tmp);
        [dd,ccc] = sort(tempQ,'descend');
        index(:,i) = indM(tmp(ccc(1)),:)';
   %     regret(i) = loss0-shiftmean(index(:,i))*inv(covMat(index(:,i),index(:,i)))*shiftmean(index(:,i))';
        count(i,:) = count(i-1,:);
        count(i,index(:,i)) =  count(i,index(:,i))+ 1;            
        obsM = data(i+m0,index(:,i));
        non_ind = 1:dim;
        non_ind(index(:,i)) = [];
        dataM = zeros(1,dim);
        dataM(index(:,i)) = obsM;
        Z = eye(dim);
        Z(non_ind,:) = [];
        tmp =sort(index(:,i));
        Omega{i} = Omega{i-1}*(1-lambda) +Z'*inv(covMat(tmp,tmp))*Z;
        sumU(i,:) = sumU(i-1,:).*(1-lambda) + dataM(tmp)*inv(covMat(tmp,tmp))*Z;
  %      invOmega{i} = inv(Omega{i});
        meanU(i,:)= sumU(i,:)*inv(Omega{i}+1e-4*eye(dim));

%  Test(i) = (1-(1-lambda)^(sqrt(i))).*sumU(i,:)*inv(Omega{i})*sumU(i,:)';   
 % Test(i) = meanU(i,:)*Omega{i}*meanU(i,:)';
    end
 %   Test(i) = 1./sum(ww(1:i).^2)*(sum(ww(1:i)))*meanU(i,:)*Omega{i}*meanU(i,:)';
  Test(i) = meanU(i,:)*Omega{i}*meanU(i,:)';
    if i > tau
    if Test(i) > limit
        runlen = i-tau;
        break;
    end
    end
       Q(i,:) = findM1(meanU(i,:),indM,invSigM,sqrt(diag(inv(Omega{i}+1e-4*eye(dim)))),loss(i));
%      QU(i,:) = abs(meanU(i,:)) + sqrt(1./2./count(i,:)*log(alpha*dim*i^4));
%      QU(i,:) = Q(i,:) + sqrt(diag(inv(Omega{i})).*log(max(1.1,loss(i))))';
%      QL(i,:) = Q(i,:) - sqrt(diag(inv(Omega{i})).*log(max(1.1,loss(i))))';          
end
count1 = count(runlen+tau,:);
meanU1 = meanU(runlen+tau,:);
% if runlen<kmax-m0
% diag = find(Test==1);
% else
%     diag = 0;
% end
end