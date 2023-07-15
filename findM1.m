function [Q,Q1]  = findM1(meanUt,indM,invSigM,stdM,losst)
[nM,m] = size(indM);
Q = zeros(1,nM);
Q1 = Q;
for i = 1:nM
%Q1(i) = 0;
Q1(i) = meanUt(indM(i,:))*invSigM{i}*meanUt(indM(i,:))';
 tmp = zeros(m,m);
 tmpPhi = invSigM{i};
 for j = 1:m
     for k = 1:m
         tmp(j,k) = sqrt(losst).*abs(tmpPhi(j,k)*meanUt(indM(i,j))*stdM(indM(i,k)))+sqrt(losst).*abs(tmpPhi(j,k)*meanUt(indM(i,k))*stdM(indM(i,j)))+...
         losst.*abs(tmpPhi(j,k)*stdM(indM(i,k))*stdM(indM(i,j)));
     end
 end
 Q(i) = Q1(i)+sum(tmp(:));
%Q(i)=Q1(i);
end
end