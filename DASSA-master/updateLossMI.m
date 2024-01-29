function[DeltaI]= updateLossMI(Prob_1tild,Prob_2_1tild,DeltaI,I1,I2)
oldI = max(I1,I2);
Removed=(DeltaI(2,:)==oldI | DeltaI(3,:)==oldI);
DeltaI(:,Removed) = [];

DeltaI(2,:) = DeltaI(2,:) - (DeltaI(2,:) > oldI);
DeltaI(3,:) = DeltaI(3,:) - (DeltaI(3,:) > oldI);

newI = min(I1,I2);
Index = find(DeltaI(2,:)==newI | DeltaI(3,:)==newI);
% npairs = size(Index,2);
% Extra_Delta= [zeros(1,npairs) ; ones(1,npairs).*newI ; Index];

%the new implementation
if ~isempty(Index)
	DeltaI_2=DeltaI(:,Index);
	P1=Prob_1tild(DeltaI_2(2,:))';
	P2=Prob_1tild(DeltaI_2(3,:))';
	Pstar=P1+P2;
	Pi1=P1./Pstar;
	Pi2=P2./Pstar;
	Py1=Prob_2_1tild(:,DeltaI_2(2,:));
	Py2=Prob_2_1tild(:,DeltaI_2(3,:));
	%Ptild=Py1*diag(Pi1)+Py2*diag(Pi2);
	%size(Py1)
	%size(Pi1)
	Ptild=Py1.*repmat(Pi1,size(Py1,1),1)+Py2.*repmat(Pi2,size(Py2,1),1);

	l=Ptild==0;
	Ptild(l)=0.001;
	temp1=Py1./Ptild;
	temp2=Py2./Ptild;
	temp1(temp1==0)=1;
	temp2(temp2==0)=1;
	KL1=sum(Py1.*log(temp1));
	KL2=sum(Py2.*log(temp2));
	Djs=Pi1.*KL1+Pi2.*KL2;

	DeltaI_2(1,:)=Pstar.*Djs;
	leng=size(DeltaI_2,2);
	DeltaI(:,Index)=DeltaI_2(:,1:leng);
end
%new implementation ends.
%old implementation 
%for i = Index
%    I1 = DeltaI(2,i);
%    I2 = DeltaI(3,i);
%    Pstar = Prob_1tild(I1) + Prob_1tild(I2);
%    DeltaI(1,i) = Pstar .* Djs(Prob_2_1tild(:,I1),Prob_2_1tild(:,I2),Pstar,Prob_1tild(I1),Prob_1tild(I2));
%end
%old implementation ends.
% DeltaI = [DeltaI Extra_Delta];
return

function DJS = Djs(P1,P2,Pstar,Px1,Px2) 
Pi_1 = Px1/Pstar;
Pi_2 = Px2/Pstar;
Pbar = (Pi_1 .* P1) + (Pi_2 .* P2);
DJS = Pi_1*Dkl(P1,Pbar) + Pi_2*Dkl(P2,Pbar);
return
