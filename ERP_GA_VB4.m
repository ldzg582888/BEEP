function [C,Z1,Z2,Tao,A1,A2,errflag]=ERP_GA_VB4(X1,X2,deltaL)
% solve the problem in VB framework
% p(X1|C A1 Z1 Tao Psi)p(C)p(Z1)p(A1|alpha1)p(alpha1);
% p(X2|C A2 Z2 Tao Psi)p(C)p(Z2)p(A2|alpha2)p(alpha2);

% X1: EEG data from experimental condition1; format: chan*len*trial
% X2: EEG data from experimental condition2; format: chan*len*trial 
% deltaL: stop criterion (less than 10^-6 is recommanded)
% If you get any question, please contact with Chaohua Wu
% (xs.wuchaohua@gmail.com)
[chan1,len1,trial1] = size(X1);
[chan2,len2,trial2] = size(X2);
if (chan1~=chan2)||(len1~=len2)
    error('dimension of X1 and X2 disagree!');
end

nc = chan1;
errflag = 0;
nzflag = ones(1,nc); % a flag to label whether the amplitude factor is zero
%%%% initialize %%%%%


%%%% phi1&phi2
phi1 = zeros(chan1); % spontaneous EEG covariance matrix in condition 1

X1_mean = mean(X1,3);
for it = 1:trial1
    temp = X1(:,:,it)-X1_mean;
    phi1 = phi1 + temp*temp';
    
end
phi1 = phi1/(trial1*len1);

phi2 = zeros(chan2); % spontaneous EEG covariance matrix in condition 2
X2_mean = mean(X2,3);
for it = 1:trial2
    temp = X2(:,:,it)-X2_mean;
    phi2 = phi2 + temp*temp';
    
end
phi2 = phi2/(trial2*len2);

%%%% Z1&Z2
[U,S,V] = svd(X1_mean,'econ');
Z1 = V(:,1:nc)'; % ERP time course
Z2 = Z1;

%%%% A
A1 = eye(nc); % amplitude factors
A2 = eye(nc);

%%%% tao
Tao = zeros(1,nc); % latency shift
search_inv = 30; % latency shift interval


EC = U; % EC: posterior expectation of spatial filter
% ECcov = zeros(chan1,chan1,nc);

EZ1 = Z1; % EZ: posterior expectation of ERP time course 
EZ1sqr = Z1.^2;
EZ2 = Z2;
EZ2sqr = Z2.^2;

EA1 = A1;  % EA: posterior expectation of amplitude factors 
EA1sqr = A1.^2;
EA2 = A2;
EA2sqr = A2.^2;

Ealpha1 = ones(1,nc); % alpha: parameter in left truncted normal distribution 
%Ealpha2 = ones(1,nc);

Einvphi1 = inv(phi1); % posterior expectation of the inverse of spontaneous EEG covariance matrix
Einvphi2 = inv(phi2);
lowerbound0 =0;
dZ = 1;
dl = 1;
i =0;
%for j=1:60
while (((i<650)||(dZ>deltaL)) && (i<1500))%%%% iteration at least 650 times, at most 1500 times

     i = i+1;   
  
%%%% update C

    [EC,ECcov,ECsqr] = updateC(X1_mean,X2_mean, EC, EZ1, EZ2, EA1, EA2, Einvphi1, Einvphi2, EZ1sqr, EZ2sqr, EA1sqr, EA2sqr,nc,trial1,nzflag);
    % ECcov: posterior expectation of cov(C)
    % ECsqr: posterior expectation of C*C^T
     pEZ1 = EZ1;pEZ2 = EZ2;   
     [EZ1,EZ1sqr,EZ2,EZ2sqr] = updateZ(X1_mean,X2_mean,EA1,EA2,EA1sqr,EA2sqr,EC,ECcov,Tao,EZ1,EZ2,Einvphi1,Einvphi2,nc,trial1,nzflag);
     dZ = (norm(EZ1-pEZ1,'fro')/norm(pEZ1,'fro'))^2+(norm(EZ2-pEZ2,'fro')/norm(pEZ2,'fro'))^2;
     % update rate of ERP time course 
     Tao = updateTao(X2_mean,EC,EA2,EZ2,EZ1,nc,Einvphi2,nzflag);
     
    [EA1,EA1sqr,SigmaA1,MuA1,nzflag1] = updateA(X1_mean,EC, ECcov,EA1, Einvphi1, EZ1, EZ1sqr, Ealpha1,nc,trial1,nzflag);
     
    [EA2,EA2sqr,SigmaA2,MuA2,nzflag2] = updateA(X2_mean,EC, ECcov,EA2, Einvphi2, EZ2, EZ2sqr, Ealpha1,nc,trial1,nzflag);
    nzflag = nzflag1& nzflag2;
    
  
     Ealpha1 = updatealpha(EA1sqr,EA2sqr,nc,nzflag);

     Einvphi1 = updateinvphi(X1,EC,ECcov,ECsqr,EA1,EA1sqr,EZ1,EZ1sqr,nc);

     Einvphi2 = updateinvphi(X2,EC,ECcov,ECsqr,EA2,EA2sqr,EZ2,EZ2sqr,nc);
    
    

      lower_bound = lb(X1,X2,EC, ECsqr, ECcov, EA1, EA2, EA1sqr, EA2sqr, EZ1, EZ2, EZ1sqr, EZ2sqr, Ealpha1, Einvphi1, Einvphi2, nc, SigmaA1,MuA1,SigmaA2,MuA2,nzflag);
      % VB lower bound
       dl = (lower_bound - lowerbound0)/(abs(lowerbound0)+10^-30);
       if ((dl <0)&&(i>1))
         errflag =1;
       end
       lowerbound0 = lower_bound;

    disp(['iteration ',num2str(i),' deltaL = ',num2str(dl), ' dZ = ', num2str(dZ)]);
end
% normalize
normEC = sqrt(sum(EC.^2));
normEZ1 = sqrt(sum(EZ1.^2,2));
normEZ2 = sqrt(sum(EZ2.^2,2));
for i = 1:nc
  if nzflag(i) == 1
      EC(:,i) = EC(:,i)./normEC(i);
      EZ1(i,:) = EZ1(i,:)./normEZ1(i);
      EZ2(i,:) = EZ2(i,:)./normEZ2(i);
      EA1(i,i) = EA1(i,i)*normEC(i)*normEZ1(i);
      EA2(i,i) = EA2(i,i)*normEC(i)*normEZ2(i);
  end
end
%end
C = EC;
Z1 = EZ1;
Z2 = EZ2;
A1 = EA1;
A2 = EA2;
return

%%%% updata C function
function [EC,ECcov,ECsqr] = updateC(X1mean,X2mean, EC, EZ1, EZ2, EA1, EA2, Einvphi1, Einvphi2, EZ1sqr, EZ2sqr, EA1sqr, EA2sqr,nc,trial,nzflag)
[chan1,len1] = size(X1mean);   
invphi1 = Einvphi1;
invphi2 = Einvphi2;
tempEC = zeros(chan1,nc);
tempECcov = zeros(chan1,chan1,nc);
tempECsqr = zeros(chan1,chan1,nc);

for i = 1:nc
    if nzflag(i)==1
        R1 = trial*(X1mean-EC*EA1*EZ1+EC(:,i)*EA1(i,i)*EZ1(i,:));
        R1 = EA1(i,i)*EZ1(i,:)*R1';
     %   R1 = R1./len1;
        R2 = trial*(X2mean-EC*EA2*EZ2+EC(:,i)*EA2(i,i)*EZ2(i,:));
        R2 = EA2(i,i)*EZ2(i,:)*R2';
    %    R2 = R2./len1;

        xi1 = trial*EA1sqr(i,i)*sum(EZ1sqr(i,:));
        xi2 = trial*EA2sqr(i,i)*sum(EZ2sqr(i,:));

        sigmaC = (xi1*invphi1+xi2*invphi2+eye(nc))^-1;    
        muC = ((R1*invphi1+R2*invphi2)*sigmaC)';

        tempEC(:,i) = muC;

        tempECcov(:,:,i) = sigmaC; 
        tempECsqr(:,:,i) = sigmaC+muC*muC';

    end
    
end
EC = tempEC;
ECcov = tempECcov;
ECsqr = tempECsqr;
return

%%%% update   A
function [EA1,EA1sqr,SigmaA1,MuA1,nzflag] = updateA(X1mean,EC, ECcov, EA1, Einvphi1, EZ1, EZ1sqr, Ealpha1,nc,trial,nzflag)

[chan1,len1] = size(X1mean); 
invphi1 = Einvphi1;
tempEA1 = zeros(nc);
tempEA1sqr = zeros(nc); 
SigmaA1 = zeros(1,nc);
MuA1 = zeros(1,nc);
tempnzflag = ones(1,nc);
powz = zeros(sum(nzflag == 1),1);

for i = 1:nc
  if nzflag(i) == 1
 
      R1 = trial*(X1mean-EC*EA1*EZ1+EC(:,i)*EA1(i,i)*EZ1(i,:));
      p1 = (EZ1(i,:)*R1')*invphi1*EC(:,i);

      sigmaA1 = (Ealpha1(i)+trial*sum(EZ1sqr(i,:))*(trace(invphi1*ECcov(:,:,i))+EC(:,i)'*invphi1*EC(:,i)))^-1; 
      muA1 = p1*sigmaA1;
      [tempe,tempv]= myltnormal(muA1,sigmaA1,0);
      tempEA1(i,i) = tempe;
      tempEA1sqr(i,i) = tempv+tempe^2;
      powz(i)= norm(EC(:,i)*tempEA1(i,i)*EZ1(i,:),'fro')^2;
%       if tempe <5*10^-5
%           tempnzflag(i) = 0;
%       end
      SigmaA1(i) = sigmaA1;
      MuA1(i) = muA1;
  end
  

end
powmax = max(powz);
for i = 1:nc
  if nzflag(i) == 1
      if ((powz(i)/powmax)<10^-3)
          tempnzflag(i) = 0;
          tempEA1(i,i) = 0;
          tempEA1sqr(i,i) = 0;
      end
  end
      
end
nzflag = nzflag&tempnzflag;
EA1 = tempEA1;
EA1sqr = tempEA1sqr;
return
%%
%%%% update Z1
function [EZ1,EZ1sqr,EZ2,EZ2sqr] = updateZ(X1mean,X2mean,EA1,EA2,EA1sqr, EA2sqr, EC,ECcov,Tao,EZ1,EZ2,Einvphi1,Einvphi2,nc,trial,nzflag)
        [chan1,len1] = size(X1mean); 
        tempEZ1 = zeros(nc,len1);
        tempEZ1sqr = zeros(nc,len1);
        tempEZ2 = zeros(nc,len1);
        tempEZ2sqr = zeros(nc,len1);
        invphi1 = Einvphi1;
        invphi2 = Einvphi2;

        for i = 1:nc
            if nzflag(i) == 1
            R1 = trial*(X1mean-EC*EA1*EZ1+EC(:,i)*EA1(i,i)*EZ1(i,:));
            R2 = trial*(X2mean-EC*EA2*EZ2+EC(:,i)*EA2(i,i)*EZ2(i,:));
            R2 = circshift(R2,[0,-Tao(i)]);

            lamda1 = EC(:,i)*EA1(i,i);
            lamda2 = EC(:,i)*EA2(i,i);

            p1 = trial*EA1sqr(i,i)*(trace(invphi1*ECcov(:,:,i))+EC(:,i)'*invphi1*EC(:,i));
            p2 = trial*EA2sqr(i,i)*(trace(invphi2*ECcov(:,:,i))+EC(:,i)'*invphi2*EC(:,i));

            sigmaZ1 = (p1+p2+2)^-1;
            muZ1 = sigmaZ1*(lamda1'*invphi1*R1+lamda2'*invphi2*R2);

            tempEZ1(i,:) = muZ1;
            tempEZ1sqr(i,:) = sigmaZ1+muZ1.^2;
            tempEZ2(i,:) = circshift(muZ1,[0,Tao(i)]);
            tempEZ2sqr(i,:) = circshift(tempEZ1sqr(i,:),[0,Tao(i)]);
            end
                

        end

        EZ1 = tempEZ1;
        EZ1sqr = tempEZ1sqr;

        EZ2 = tempEZ2;
        EZ2sqr = tempEZ2sqr;
return

%%%% update alpha
function Ealpha1 = updatealpha(EA1sqr,EA2sqr,nc,nzflag)
Ealpha1 = zeros(nc,1);
 for i = 1:nc
     if nzflag(i) == 1
         
          Ealpha1(i) = 2/(EA1sqr(i,i)+EA2sqr(i,i));
   
     end
 end
% Ealpha1 = nc./diag(EA1);
return

%%%% in question

function Einvphi1 = updateinvphi(X1,EC,ECcov,ECsqr,EA1,EA1sqr,EZ1,EZ1sqr,nc)
[chan1,len1,trial1] = size(X1);
W = zeros(chan1);
ECAZ = EC*EA1*EZ1;
S = zeros(chan1);
for i = 1:nc
    for j = 1:nc
        if i == j
            S = S+sum(EZ1sqr(i,:))*EA1sqr(i,i)*ECsqr(:,:,i);
        else 
            S = S+EZ1(i,:)*EZ1(j,:)'*(EA1(i,i)*EA1(j,j))*(EC(:,i)*EC(:,j)');
        end
    end
end

for i = 1:trial1
    W = W + X1(:,:,i)*X1(:,:,i)'-ECAZ*X1(:,:,i)'-X1(:,:,i)*ECAZ'+S;
end
Einvphi1 = (nc+len1*trial1)*inv(W);
return


%%%% update Tao
function Tao = updateTao(X2mean,EC,EA2,EZ2,EZ1,nc,Einvphi2,nzflag)

	searchintv = [-30:30];
	searchlen = length(searchintv);
	tg = zeros(1,searchlen);
	tempTao = zeros(1,nc);
 %   EZ2sqr = zeros(size(EZ1sqr));
    %%%% ����ƶ�component 1

	for i = 1:nc
        if nzflag(i) == 1
            R2 = X2mean-EC*EA2*EZ2+EC(:,i)*EA2(i,i)*EZ2(i,:);
            for j = 1:searchlen
                tempEZ = circshift(EZ1(i,:),[0,searchintv(j)]);
                tg(j) = tempEZ*R2'*Einvphi2*EC(:,i)*EA2(i,i);
            end

            [Y,I] = max(tg);
            tempTao(i) = searchintv(I);
        end
	end
	Tao = tempTao;
% 	for i = 1:nc
% 		EZ2(i,:) = circshift(EZ1(i,:),[0,Tao(i)]);
%         EZ2sqr(i,:) = circshift(EZ1sqr(i,:),[0,Tao(i)]);
% 	end
return  

    
function lower_bound = lb(X1,X2,EC, ECsqr, ECcov, EA1, EA2, EA1sqr, EA2sqr, EZ1, EZ2, EZ1sqr, EZ2sqr, Ealpha1, Einvphi1, Einvphi2, nc, SigmaA1,MuA1,SigmaA2,MuA2,nzflag)
%%%% term1
[chan1,len1,trial1] = size(X1);
ECAZ1 = EC*EA1*EZ1;
ECAZ2 = EC*EA2*EZ2;
% 
%%% S = VAR(C*A*Z)
S1 = zeros(chan1);
ww1 = zeros(chan1);
for i = 1:nc
    for j = 1:nc
        if i == j
            S1 = S1+sum(EZ1sqr(i,:))*EA1sqr(i,i)*ECsqr(:,:,i);
        else 
            S1 = S1+EZ1(i,:)*EZ1(j,:)'*(EA1(i,i)*EA1(j,j))*(EC(:,i)*EC(:,j)');
        end
    end
end
for i = 1:trial1
    ww1 = ww1 + X1(:,:,i)*X1(:,:,i)'-ECAZ1*X1(:,:,i)'-X1(:,:,i)*ECAZ1'+S1;
end

S2 = zeros(chan1);
ww2 = zeros(chan1);
for i = 1:nc
    for j = 1:nc
        if i == j
            S2 = S2+sum(EZ2sqr(i,:))*EA2sqr(i,i)*ECsqr(:,:,i);
        else 
            S2 = S2+EZ2(i,:)*EZ2(j,:)'*(EA2(i,i)*EA2(j,j))*(EC(:,i)*EC(:,j)');
        end
    end
end
for i = 1:trial1
    ww2 = ww2 + X2(:,:,i)*X2(:,:,i)'-ECAZ2*X2(:,:,i)'-X2(:,:,i)*ECAZ2'+S2;
end

w1 = inv(Einvphi1/(nc+len1*trial1));
w2 = inv(Einvphi2/(nc+len1*trial1));
term1 = -0.5*trial1*len1*log(det(w1))-0.5*trace(ww1*Einvphi1);
term2 = -0.5*trial1*len1*log(det(w2))-0.5*trace(ww2*Einvphi2);


%%%% term3
term3 = 0;
for i = 1:nc
    if nzflag(i) == 1
    term3 = term3+trace(ECsqr(:,:,i));
    end
end
term3 = -0.5*term3;

%%%% term4
term4 = 0;
for i = 1:nc
    if nzflag(i) == 1
       term4 = term4 -0.5*sum(EZ1sqr(i,:))-0.5*sum(EZ2sqr(i,:));
    end
end
%%%% term5
term5 = 0;
dea1sqr = diag(EA1sqr);
dea2sqr = diag(EA2sqr);

for i = 1:nc
    if nzflag(i) == 1
    term5 = term5 -0.5*(Ealpha1(i)*dea1sqr(i)+Ealpha1(i)*dea2sqr(i)+2*log(dea1sqr(i)+dea2sqr(i)));

    end
end
%%%% term6
term6 = 0;
for i = 1:nc
    if nzflag(i) == 1
        term6 =term6+log(dea1sqr(i)+dea2sqr(i));
    end
end
term7 = 0.5*log(det(w1))+0.5*log(det(w2));
%%%% log C
term8 = 0;
for i=1:nc
    if nzflag(i) == 1
    term8 = term8 -0.5*log(10^-150+det(ECcov(:,:,i)));
    end
end
term8 = -term8;
%%%% log Z
term9 = 0;
for i = 1:nc
    if nzflag(i) == 1
    term9 = term9 -0.5*sum(log(10^-150+EZ1sqr(i,:)-EZ1(i,:).^2))-0.5*sum(10^-150+log(EZ2sqr(i,:)-EZ2(i,:).^2));
    end
end
term9 = -term9;


%%%% log A
term10 = 0;
for i = 1:nc
    if nzflag(i) == 1
    g = 0.5*erfc(-MuA1(i)/sqrt(2*SigmaA1(i)));
    term10 = term10 -log(g)-0.5*log(SigmaA1(i))-0.5*((dea1sqr(i)-EA1(i,i)^2)/SigmaA1(i)+(EA1(i,i)-MuA1(i))^2/SigmaA1(i));
    end
end

for i = 1:nc
    if nzflag(i) == 1
    g = 0.5*erfc(-MuA2(i)/sqrt(2*SigmaA2(i)));
    term10 = term10 -log(g)-0.5*log(SigmaA2(i))-0.5*((dea2sqr(i)-EA2(i,i)^2)/SigmaA2(i)+(EA2(i,i)-MuA2(i))^2/SigmaA2(i));
    end
end
term10 = -term10;
%%%% log alpha
term11 = 0;
for i = 1:nc
    if nzflag(i) == 1
    term11 = term11+log(dea1sqr(i)+dea2sqr(i));
    end
end
term11 = -term11;

%%%% log phi
term12 = 0.5*(chan1+1)*(log(det(w1))+log(det(w2)));
term12 = -term12;
    
lower_bound = term1+term2+term3+term4+term5+term6+term7+term8+term9+term10+term11+term12;

return


