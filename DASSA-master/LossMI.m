function[DeltaI]= LossMI(Prob_1tild,Prob_2_1tild,nx)
n1 = size(Prob_1tild,1);
if nargin > 2
    n1 = nx;
end
npairs = n1 * (n1-1) / 2;
DeltaI = zeros(3,npairs);
tmp = 1:n1;
Index = {tmp tmp};
[I1,I2] = ndgrid(Index{:});
I = [I1(:),I2(:)];
repeated=(I1<=I2);
I(repeated,:)=[];
DeltaI(2:3,:) = I';

%from here it's the new version of vectorized implementation.
P1=Prob_1tild(DeltaI(2,:))';
P2=Prob_1tild(DeltaI(3,:))';
Pstar=P1+P2;
Pi1=P1./Pstar;
Pi2=P2./Pstar;
Py1=Prob_2_1tild(:,DeltaI(2,:));
Py2=Prob_2_1tild(:,DeltaI(3,:));
%Ptild=Py1*diag(Pi1)+Py2*diag(Pi2);
%size(Py1)
%size(Pi1)
Ptild=Py1.*repmat(Pi1,size(Py1,1),1)+Py2.*repmat(Pi2,size(Py2,1),1);

l=(full(Ptild==0));
Ptild(l)=0.001;
temp1=Py1./Ptild;
temp2=Py2./Ptild;
clear Ptild;
temp1(temp1==0)=1;
temp2(temp2==0)=1;
KL1=sum(Py1.*log(temp1));
clear Py1;
KL2=sum(Py2.*log(temp2));
clear Py2;
Djs=Pi1.*KL1+Pi2.*KL2;

DeltaI(1,:)=Pstar.*Djs;

%the new implementation ends here

%Old implementation starts here
%for i = 1 : npairs
%    I1 = DeltaI(2,i);
%    I2 = DeltaI(3,i);
%    Pstar = Prob_1tild(I1) + Prob_1tild(I2);
%    DeltaI(1,i) = Pstar .* Djs(Prob_2_1tild(:,I1),Prob_2_1tild(:,I2),Pstar,Prob_1tild(I1),Prob_1tild(I2));
%end
%old implementation till here.
return

function DJS = Djs(P1,P2,Pstar,Px1,Px2) 
Pi_1 = Px1/Pstar;
Pi_2 = Px2/Pstar;
Pbar = (Pi_1 .* P1) + (Pi_2 .* P2);
DJS = Pi_1*Dkl(P1,Pbar) + Pi_2*Dkl(P2,Pbar);
return
