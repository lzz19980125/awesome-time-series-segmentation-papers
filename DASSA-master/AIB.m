function AIB(Address, Saveadd, MinW, N_mdl, CalcMode)
clc
%% Input Parameteres
% Address = '../../Journals/BMSJ/BMSJ1890_IB.txt';
% Saveadd = '../../Journals/BMSJ/IB/';
% MinW = '604800';
% CalcMode = '1';
% N_mdl = '10';

Wmin = str2num(MinW);
N_mdl = str2num(N_mdl);
CalcMode = str2num(CalcMode);
Wstep = Wmin;
L1 = 1;
L2 = 1;
AIBMode = 1;
%%
if ~exist(Saveadd,'dir')
	mkdir(Saveadd);
end
%% Definition
'Definition...'
[X,XID,Y,YID,CatIndxtilde,WinIndxtilde,Timestamp, BeginTime,EndTime] = Definition(Wmin,Wstep,Address);
nx = size(CatIndxtilde,2);

%% Computing D
'Computing D...'
D = ComputeD(XID,Y,YID,Timestamp);
%% First Step
'Calculating Xtilde..'
[Xtilde,CatIndxtilde,Prob_X,Prob_Xtild,Prob_Xtilde_X] = CalculateXtilde(CalcMode,L1,XID,Y,YID,WinIndxtilde,CatIndxtilde,Timestamp,Saveadd,D,nx,AIBMode,N_mdl);
%% Second Step
[Ytilde,WinIndxtilde,Prob_Xtilde_Y,Prob_Ytild,Prob_Ytilde_Y,Prob_Y] = CalculateYtilde(L2,Xtilde,Y,YID,WinIndxtilde,CatIndxtilde,Timestamp,Saveadd,D);
%% Saving Final resultis
Prob_X_YYtilde = 0;
size(Prob_Xtilde_Y)
% save(strcat(Saveadd,'Output'));
Saving(Saveadd,(Xtilde-1),[X,Timestamp],(Ytilde-1),[(YID-1),Y],Prob_Xtilde_Y,Prob_X_YYtilde,[BeginTime EndTime]);
'done!'
exit

%% save results
function Saving(Saveadd,Xtilde,X,Ytilde,Y,Prob_Xtilde_Y,Prob_X_YYtilde,Time)
'Saving...'
%tic
Address = strcat(Saveadd,'Xtilde.txt');
dlmwrite(Address, Xtilde,'\t')

Address = strcat(Saveadd,'X.txt');
dlmwrite(Address, X,'\t')

Address = strcat(Saveadd,'Ytilde.txt');
dlmwrite(Address, Ytilde,'\t')

Address = strcat(Saveadd,'Y.txt');
dlmwrite(Address, Y,'\t')
Address = strcat(Saveadd,'P_Xtilde_Y.txt');
dlmwrite(Address, Prob_Xtilde_Y,'\t')

Address = strcat(Saveadd,'BeginEndTime.txt');
dlmwrite(Address, Time,'\t')
return
%% Other functions

function[X,XID,Y,YID,CatIndx,WinIndxtilde,Timestamp, Begintime, Endtime] = Definition(Wmin,Wstep,Address)
%% Construct X 
[X,XID,Timestamp,CatIndx,  Begintime, Endtime] = ConstX(Address);
Wmax = Endtime - Begintime +0.01;
% %% Construct X_tilde as X
% [CatIndxtilde] = ConstXtilde(XID,CatIndx);
%% Construct Y
[Y,YID] = ConstY(Begintime,Endtime,Wmin,Wmax,Wstep);
%% Construct Y_tilde as Y
[WinIndxtilde] = ConstYtilde(YID);
return

function [X,XID,Timestamp,CatIndx, Begintime, Endtime] = ConstX(Address)
%% Import Data as X and Timstamp
Data = load(Address);
X = Data(:,1:(end-1));
Timestamp = Data(:,end);
Begintime = min(Timestamp);
Endtime = max(Timestamp);
%% Find the Categories in X
XID = zeros(size(X,1),1);
ID = zeros(size(X,1),1); % The array of not categoriezed rows
Catnumber = 1;
CatIndx = [];
while any(ID==0)
    NotVisited = find(ID == 0);
    Index = find(ismember(X,X(NotVisited(1),:),'rows'));
    XID(Index) = Catnumber;
    ID(Index) = 1;
    CatIndx = [CatIndx,{Index}];
    Catnumber = Catnumber + 1;
end
return

function [Y,YID] = ConstY(BeginTime,EndTime,Wmin,Wmax,Wstep)
% DefineY is a function to construct the Y list, which is the set of all possible time windows
% ConstY has the 4D output Y. 
% each Dimension represent, 1- window size  2- Window number 3- Start time 4- End time
Y = [];
YID = [];
Wsize = Wmin;
Overlap = Wmin;
% if EndTime > BeginTime
%     Duration = EndTime-BeginTime;
% else
%     Duration = 0;
% end
while (Wsize <= Wmax)
%         NW= ceil(Duration / Wsize);
        StartID = size(Y,1);
        [L,ID] = GetWNumber(Overlap,Wsize,BeginTime,EndTime,StartID);
        Y = [Y;L];
        YID = [YID;ID];
        Wsize = Wsize + Wstep;
end
return

function [WinIndxtilde] = ConstYtilde(YID)
% Ytilde = YID;
WinIndxtilde = num2cell(YID');
return

function [L,ID] = GetWNumber(Overlap,Wsize,BeginTime,EndTime,StartID)
% GetWNumber is a function to concatenate the window size and window number
% the format of each element in L is : 'window size'+ 'N' + 'Window number'
% the output is an string list
    ID = [];
    L = [];
    Begin = BeginTime;
    Endtmp = Begin +  Wsize;
    End = EndTime;
    i = 0;
    while (Endtmp <= EndTime)
        End = Endtmp;
        i = i + 1;
        ID = [ID;StartID + i];
        L = [L;[Wsize,i,Begin,End]];
        Begin = Begin + Overlap;
        Endtmp = Begin + Wsize;
    end 
    if End < EndTime
       i = i +1;
       ID = [ID;StartID + i];
       L = [L;[Wsize,i,Begin,EndTime]];
    end
return

function D = ComputeD(XID,Y,YID,Timestamp)
Smallest = min(Y(:,1));
SmallIndx = find(YID(Y(:,1)==Smallest));
SmallY = Y(SmallIndx,:);

XTime = repmat(Timestamp,1,size(SmallY,1));
YBegin = repmat(SmallY(:,end-1)',size(XID,1),1);
YEnd = repmat(SmallY(:,end)',size(XID,1),1);

[~,Index] = max(((XTime>= YBegin) & (XTime < YEnd)),[],2);
Win = SmallIndx(Index).*(XTime(:,1)~=YEnd(1,end)) + (SmallIndx(end)).*(XTime(:,1)==YEnd(1,end));
D = [XID,Win];
return

%%
function [Xtilde,CatIndxtilde,Prob_X,Prob_Xtild,Prob_Xtilde_X] = CalculateXtilde(CalcMode,L1,XID,Y,YID,WinIndxtilde,CatIndxtilde,Timestamp,Saveadd,D,nx,AIBMode,N_mdl)
%% Initialization
'X Initialization...'
% Xtilde is constracted as X
Xtilde = XID;
[~,Prob_X,Prob_Y_X] = InitialProb(XID,Y,YID,Timestamp,CatIndxtilde,WinIndxtilde);
[Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde] = InitialCPT(Prob_X,Prob_Y_X);
% DeltaI : 3 x n matrix,  DeltaI(1,x) is the DelataI when DeltaI(2,x) and
% DeltaI(3,x) is merged
if AIBMode == 1
[DeltaI]= LossMI(Prob_Xtild,Prob_Y_Xtild,nx);
N = size(CatIndxtilde,2);
%% Updating
'X Updating...'
if CalcMode == 1
    NoMDL = N;
else 
    NoMDL = N_mdl
end
% Merging will be continued untile 'NoMDL' clusters remained. 
% Then MDL will be runnig inorder to find the best number of clusters
NoMDL = NoMDL - 2;
if N == size(Prob_X,1)
    NoMDL = NoMDL + 1;
end
%% Smallest
Smallest = min(Y(:,1));
SmallIndx = find(YID(Y(:,1)==Smallest));
T = size(SmallIndx,1);
% NewProb = Prob_Y_Xtild(SmallIndx,:);%./repmat(sum(Prob_Y_Xtild(SmallIndx,:),1),T,1);
% [Cost(end),model_cost(end),data_cost(end)]= MDL(N,N,T,NewProb,Prob_Xtild,Prob_X_Xtilde,Prob_Xtilde_X, D);

[I1,I2,DeltaI,CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde] = WhioutMDL_x(N,NoMDL,DeltaI,CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde);
if CalcMode==1
    [Cost,model_cost,data_cost,L_XBest]=MDL_Cost_x(NoMDL,DeltaI,CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde,I1,I2,T,Saveadd,D,SmallIndx,N,Prob_X);
    [~,~,~,CatIndxtilde,Xtilde,Prob_Xtild,~,Prob_Xtilde_X,~] = WhioutMDL_x(NoMDL+1,L_XBest-1,DeltaI,CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde);
end
end
return

function [CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde] = Merge_x(CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde,I1,I2)
IMax = max(I1,I2);
IMin = min(I1,I2);
%% Updating CatIndxtilde
Cat1 = cell2mat(CatIndxtilde(I1));
Cat2 = cell2mat(CatIndxtilde(I2));
CatNew = [Cat1;Cat2];
CatIndxtilde(IMax) = [];
CatIndxtilde(IMin) = {CatNew};
%% Updating Xtilde
Xtilde(CatNew) = IMin;
newX = (Xtilde>IMax);
Xtilde(newX) = Xtilde(newX)-1;
%% Updating Prob_Xtild and Prob_Y_Xtild
P1 = Prob_Xtild(I1);
P2 = Prob_Xtild(I2);
Pnew = P1 + P2;
Pyxnew = ((P1/Pnew).*Prob_Y_Xtild(:,I1)) + ((P2/Pnew).*Prob_Y_Xtild(:,I2));
Pxxnew = ((P1/Pnew).*Prob_X_Xtilde(:,I1)) + ((P2/Pnew).*Prob_X_Xtilde(:,I2));
% updating Prob_Xtild
Prob_Xtild(IMax) = [];
Prob_Xtild(IMin) = Pnew;
% updating Prob_Y_Xtild
Prob_Y_Xtild(:,IMax) = [];
Prob_Y_Xtild(:,IMin) = Pyxnew;
% updating Prob_X_Xtilde
Prob_X_Xtilde(:,IMax) = [];
Prob_X_Xtilde(:,IMin) = Pxxnew;
%% Updating Prob_Xtilde_X
PNew = (Prob_Xtilde_X(I1,:) | Prob_Xtilde_X(I2,:));
Prob_Xtilde_X(IMax,:) = [];
Prob_Xtilde_X(IMin,:) = PNew;
return

function[Cost,model_cost,data_cost,L_XBest] = MDL_Cost_x(NoMDL,DeltaI,CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde,I1,I2,T,Saveadd,D,SmallIndx,N,Prob_X)
Cost = zeros(NoMDL,1);
model_cost = zeros(NoMDL,1);
data_cost = zeros(NoMDL,1);

[CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde] = Merge_x(CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde,I1,I2);
[DeltaI]= updateLossMI(Prob_Xtild,Prob_Y_Xtild,DeltaI,I1,I2);
NewProb = Prob_Y_Xtild(SmallIndx,:);%./repmat(sum(Prob_Y_Xtild(SmallIndx,:),1),T,1);
%[Cost(NoMDL),model_cost(NoMDL),data_cost(NoMDL)]= MDL2(NoMDL,N,T,NewProb,Prob_Xtild,Prob_X_Xtilde,Prob_Xtilde_X, D,Prob_X);
[Cost(NoMDL),model_cost(NoMDL),data_cost(NoMDL)]= MDL(NoMDL,N,T,NewProb,Prob_Xtild,Prob_X_Xtilde,Prob_Xtilde_X, D);
%%

for i = NoMDL-1 :-1:1
    %i
    if ~isempty(DeltaI)
        [~,MinInx] = min(DeltaI(1,:));
        I1 = DeltaI(2,MinInx);
        I2 = DeltaI(3,MinInx);
    else
        I1 = 2;
        I2 = 1;
    end
    [CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde] = Merge_x(CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde,I1,I2);
    [DeltaI]= updateLossMI(Prob_Xtild,Prob_Y_Xtild,DeltaI,I1,I2);
    NewProb = Prob_Y_Xtild(SmallIndx,:);%./repmat(sum(Prob_Y_Xtild(SmallIndx,:),1),T,1);
    %[Cost(i),model_cost(i),data_cost(i)] = MDL2(i,N,T,NewProb,Prob_Xtild,Prob_X_Xtilde,Prob_Xtilde_X, D,Prob_X);
    [Cost(i),model_cost(i),data_cost(i)] = MDL(i,N,T,NewProb,Prob_Xtild,Prob_X_Xtilde,Prob_Xtilde_X, D);
%     if (NoMDL>=(Step+i)) && (Cost(i) >= Cost(i+Step))
%         NumberOfLabels = (i + Step)
%         break;
%     end
end
if N == size(Prob_X,1)
   Cost(1)=[];
   model_cost(1)=[];
   data_cost(1)=[];
end
% save(strcat(Saveadd,'CostX.mat'));
figmode = 0;
if figmode == 1
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MDL Cost 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = Cost;
X = 1:(size(Y,1)*size(Y,2));
[m,i] = min(Y);
minX = X(i);
Margin = 2.5;
SX = max(min(X),minX-Margin);
EX = min(max(X),minX+Margin);
SY = min(Y(SX:EX));
EY = max(Y(SX:EX));
h=figure(1);
plot(Cost,'LineWidth',2);set(gca,'FontSize',20);
xlabel('Number of Clusters (l)','FontSize',20)
ylabel('MDL Cost','FontSize',20)
grid on;
magnifyOnFigure(...
        h,...
        'units', 'pixels',...
        'initialPositionSecondaryAxes', [max(X)-130 max(Y)-150 100],...
        'initialPositionMagnifier',     [160 0 (EX-SX)*2 (EY-SY)],...    
        'mode', 'interactive',...    
        'displayLinkStyle', 'straight',...        
        'edgeWidth', 2,...
        'edgeColor', 'black'); 
saveas(h,strcat(Saveadd,'Cost.fig'));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model Cost 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure(2);
plot(model_cost,'LineWidth',2);
xlabel('Number of Clusters (l)','FontSize',20)
ylabel('Model Cost','FontSize',20)
grid on;    
saveas(h,strcat(Saveadd,'model_cost.fig'));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Data Cost 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure(3);
plot(data_cost,'LineWidth',2);
xlabel('Number of Clusters (l)','FontSize',20)
ylabel('Data Cost','FontSize',20)
grid on;
saveas(h,strcat(Saveadd,'data_cost.fig'));
end

[~,L_XBest] = min(Cost);
return

function [I1,I2,DeltaI,CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde] = WhioutMDL_x(StartInx,EndInx,DeltaI,CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde)
for i = StartInx-1 :-1:(EndInx+1)
    %i
    if ~isempty(DeltaI)
        [~,MinInx] = min(DeltaI(1,:));
    %     Plot = [Plot , Min];
        I1 = DeltaI(2,MinInx);
        I2 = DeltaI(3,MinInx);
    else
        I1 = 2;
        I2 = 1;
    end
    [CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde] = Merge_x(CatIndxtilde,Xtilde,Prob_Xtild,Prob_Y_Xtild,Prob_Xtilde_X,Prob_X_Xtilde,I1,I2);
    [DeltaI]= updateLossMI(Prob_Xtild,Prob_Y_Xtild,DeltaI,I1,I2);
end
%%
if ~isempty(DeltaI)
    [~,MinInx] = min(DeltaI(1,:));
    I1 = DeltaI(2,MinInx);
    I2 = DeltaI(3,MinInx);
else
    I1 = 2;
    I2 = 1;
end
return


function [Prob_XY,Prob_X,Prob_Y_X] = InitialProb(XID,Y,YID,Timestamp,CatIndxtilde,WinIndxtilde)
Xsize = size(CatIndxtilde,2) + 1;
Prob_XY = zeros(Xsize,size(WinIndxtilde,2)); 

  
for i = 1 : size(CatIndxtilde,2)
    Index = CatIndxtilde{i};
    XTime = repmat(Timestamp(Index),1,size(WinIndxtilde,2) );
    YBegin = repmat(Y(:,end-1)',size(Index,1),1);
    YEnd = repmat(Y(:,end)',size(Index,1),1); 
    Num = sum((XTime>= YBegin) & (XTime < YEnd),1)  + (YEnd(1,:)==YEnd(1,end)).*sum(XTime(:,1)==YEnd(1,end));
    Prob_XY(i,:) = Num;
end
Prob_XY(end,:) = (sum(Prob_XY,1)==0);
Prob_XY = Prob_XY./sum(sum(Prob_XY));
if all(Prob_XY(end,:)==0)
    Prob_XY(end,:) = [];
end
Prob_X = sum(Prob_XY,2);
Denominator = repmat(Prob_X,1,size(Prob_XY,2));
Prob_Y_X = (Prob_XY./(Denominator + (Denominator==0)))';
% [n1,n2] = size(Prob_XY);
% Prob_Y_X = Prob_Y_X + ((Denominator==0).*(ones(n1,n2))./n2)';
return
%% 
function [Ytilde,WinIndxtilde,Prob_X_Y,Prob_Ytild,Prob_Ytilde_Y,Prob_Y] = CalculateYtilde(L2,XID,Y,YID,WinIndxtilde,CatIndxtilde,Timestamp,Saveadd,D)
%% Initialization
'Y Initialization...'
% Ytilde is constracted as Y
Ytilde = YID;
[~,Prob_Y,Prob_X_Y] = InitialProbY(XID,Y,YID,Timestamp,CatIndxtilde,WinIndxtilde);
[Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde] = InitialCPT(Prob_Y,Prob_X_Y);
Condition = 0;
if Condition == 1
% DeltaI : 3 x n matrix,  DeltaI(1,x) is the DelataI when DeltaI(2,x) and
% DeltaI(3,x) merge
[DeltaI]= LossMI(Prob_Ytild,Prob_X_Ytild);
%% Updating
'Y Updating...'
NoMDL = 20; % Merging will be continued untile 'NoMDL' clusters remained. 
% Then MDL will be runnig inorder to find the best number of clusters
NoMDL = NoMDL + 1;
N = size(Prob_X_Y,1);
T = size(WinIndxtilde,2);

[I1,I2,DeltaI,WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde] = WhioutMDL(T,NoMDL,DeltaI,WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde);
[Cost,model_cost,data_cost,L_XBest]=MDL_Cost_y(NoMDL,DeltaI,WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde,I1,I2,T,Saveadd,D,N,Prob_Y);
[~,~,~,WinIndxtilde,Ytilde,Prob_Ytild,~,Prob_Ytilde_Y,~] = WhioutMDL_y(NoMDL+1,L_XBest-1,DeltaI,WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde);
end
return

function [WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde] = Merge_y(WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde,I1,I2)
IMax = max(I1,I2);
IMin = min(I1,I2);
%% Updating CatIndxtilde
Win1 = cell2mat(WinIndxtilde(I1));
Win2 = cell2mat(WinIndxtilde(I2));
WinNew = [Win1;Win2];
WinIndxtilde(IMax) = [];
WinIndxtilde(IMin) = {WinNew};
%% Updating Xtilde
Ytilde(WinNew) = IMin;
newY = (Ytilde>IMax);
Ytilde(newY) = Ytilde(newY)-1;
%% Updating Prob_Xtild and Prob_Y_Xtild
P1 = Prob_Ytild(I1);
P2 = Prob_Ytild(I2);
Pnew = P1 + P2;
Pyxnew = ((P1/Pnew).*Prob_X_Ytild(:,I1)) + ((P2/Pnew).*Prob_X_Ytild(:,I2));
Pyynew = ((P1/Pnew).*Prob_Y_Ytilde(:,I1)) + ((P2/Pnew).*Prob_Y_Ytilde(:,I2));
% updating Prob_Xtild
Prob_Ytild(IMax) = [];
Prob_Ytild(IMin) = Pnew;
% updating Prob_X_Ytild
Prob_X_Ytild(:,IMax) = [];
Prob_X_Ytild(:,IMin) = Pyxnew;
% updating Prob_Y_Ytilde
Prob_Y_Ytilde(:,IMax) = [];
Prob_Y_Ytilde(:,IMin) = Pyynew;
%% Updating Prob_Xtilde_X
PNew = (Prob_Ytilde_Y(I1,:) | Prob_Ytilde_Y(I2,:));
Prob_Ytilde_Y(IMax,:) = [];
Prob_Ytilde_Y(IMin,:) = PNew;
return

function[Cost,model_cost,data_cost,L_XBest] = MDL_Cost_y(NoMDL,DeltaI,WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde,I1,I2,T,Saveadd,D,N,Prob_Y)
Cost = zeros(NoMDL,1);
model_cost = zeros(NoMDL,1);
data_cost = zeros(NoMDL,1);

[WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde] = Merge_y(WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde,I1,I2);
[DeltaI]= updateLossMI(Prob_Ytild,Prob_X_Ytild,DeltaI,I1,I2);
[Cost(NoMDL),model_cost(NoMDL),data_cost(NoMDL)]= MDL_y(NoMDL,T,N,Prob_X_Ytild,Prob_Ytild,Prob_Y_Ytilde,Prob_Ytilde_Y, D);
%%

for i = NoMDL-1 :-1:1
    %i
    if ~isempty(DeltaI)
        [~,MinInx] = min(DeltaI(1,:));
        I1 = DeltaI(2,MinInx);
        I2 = DeltaI(3,MinInx);
    else
        I1 = 2;
        I2 = 1;
    end

    [WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde] = Merge_y(WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde,I1,I2);
    [DeltaI]= updateLossMI(Prob_Ytild,Prob_X_Ytild,DeltaI,I1,I2);
    [Cost(i),model_cost(i),data_cost(i)] = MDL_y(i,T,N,Prob_X_Ytild,Prob_Ytild,Prob_Y_Ytilde,Prob_Ytilde_Y, D);
end
if N == size(Prob_Y,1)
   Cost(1)=[];
   model_cost(1)=[];
   data_cost(1)=[];
end
save(strcat(Saveadd,'CostY.mat'));
% figure(4)
% plot(Cost,'r')
% legend('Cost')
% figure(5)
% plot(model_cost,'g');
% legend('model_cost')
% figure(6)
% plot(data_cost,'b');
% legend('data_cost')
[~,L_XBest] = min(Cost);
return

function [I1,I2,DeltaI,WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde] = WhioutMDL_y(StartInx,EndInx,DeltaI,WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde)
for i = StartInx-1 :-1:(EndInx+1)
    %i
    if ~isempty(DeltaI)
        [~,MinInx] = min(DeltaI(1,:));
    %     Plot = [Plot , Min];
        I1 = DeltaI(2,MinInx);
        I2 = DeltaI(3,MinInx);
    else
        I1 = 2;
        I2 = 1;
    end
    [WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde] = Merge_y(WinIndxtilde,Ytilde,Prob_Ytild,Prob_X_Ytild,Prob_Ytilde_Y,Prob_Y_Ytilde,I1,I2);
    [DeltaI]= updateLossMI(Prob_Ytild,Prob_X_Ytild,DeltaI,I1,I2);
end
%%
if ~isempty(DeltaI)
    [~,MinInx] = min(DeltaI(1,:));
    I1 = DeltaI(2,MinInx);
    I2 = DeltaI(3,MinInx);
else
    I1 = 2;
    I2 = 1;
end
return

function [Prob_YX,Prob_Y,Prob_X_Y] = InitialProbY(XID,Y,YID,Timestamp,CatIndxtilde,WinIndxtilde)
Xsize = size(CatIndxtilde,2)+1;
Prob_YX = zeros(size(WinIndxtilde,2),Xsize);

for i = 1 : size(CatIndxtilde,2)
    Index = CatIndxtilde{i};
    XTime = repmat(Timestamp(Index),1,size(WinIndxtilde,2) );
    YBegin = repmat(Y(:,end-1)',size(Index,1),1);
    YEnd = repmat(Y(:,end)',size(Index,1),1); 
    Num = sum((XTime>= YBegin) & (XTime < YEnd),1)  + (YEnd(1,:)==YEnd(1,end)).*sum(XTime(:,1)==YEnd(1,end));
    Prob_YX(:,i) = Num';
end
Prob_YX(:,end) = (sum(Prob_YX,2)==0);
Prob_YX = Prob_YX./sum(sum(Prob_YX));
if all(Prob_YX(:,end)==0)
    Prob_YX(:,end) = [];
end
Prob_Y = sum(Prob_YX,2);
Denominator = repmat(Prob_Y,1,size(Prob_YX,2));
Prob_X_Y = (Prob_YX./(Denominator + (Denominator==0)))';
[n1,n2] = size(Prob_YX);
Prob_X_Y = Prob_X_Y + ((Denominator==0).*(ones(n1,n2))./n2)';
return

function [Prob_1tild,Prob_2_1tild,Prob_1tilde_1,Prob_1_1tilde] = InitialCPT(Prob_1,Prob_2_1)
%% Computing Prob_Xtild
Prob_1tild = Prob_1;
%% Computing Prob_YXtild
[n1,n2] = size(Prob_2_1);
P1 = repmat(Prob_1',n1,1);
P1tild = repmat(Prob_1tild',n1,1);
numinator = (P1.*Prob_2_1); 

Prob_2_1tild = (P1tild ~= 0).* (numinator./(P1tild + ((P1tild == 0))));
%% Computing Prob_XtildeX
Prob_1tilde_1 = eye(size(Prob_1,1));

[n1,n2] = size(Prob_1tilde_1);
% P1 = repmat(Prob_1',n1,1);
% P1tild = repmat(Prob_1tild,1,n2);
% Prob_1_1tilde = ((Prob_1tilde_1 .* P1) ./ P1tild)';
Prob_1_1tilde = Prob_1tilde_1;
return
