function [ total_cost,model_cost,data_cost ] = MDL2( l,N,T,y_xt,xt,x_xt,xt_x, D,x)
%MDL Summary of this function goes here
%   Detailed explanation goes here

%model cost
%model_cost=log_s(l)+log_s(N)+log_s(T)+N*log_2(l)+entropy(xt);
model_cost=L_N(l)+L_N(N)+L_N(T)+N*log_2(l)+entropy(x);
for i=1:l
    %model_cost=model_cost+entropy(y_xt(:,i))+entropy(x_xt(:,i));
    model_cost=model_cost+entropy(y_xt(:,i));
end

%data description cost
data_cost=0;
[xx,y]=size(D);
[~,I]=max(xt_x,[],1);
xt2=I(D(:,1));
x2=D(:,1);
y2=D(:,2);
 
m1=(xt2-1)'.*size(x_xt,1)+x2;
m2=(xt2-1)'.*size(y_xt,1)+y2;
%Ans1 = x_xt(m1);
%if size(Ans1,1)<size(Ans1,2)
%    Ans1 = Ans1';
%end
%Ans2 = y_xt(m2);
%if size(Ans2,1)<size(Ans2,2)
%    Ans2 = Ans2';
%end
%Ans3 = xt(xt2);
%if size(Ans3,1)<size(Ans3,2)
%    Ans3 = Ans3';
%end
%data_cost=sum(-log_2(Ans1.*Ans2.*Ans3));
for i=1:xx
    data_point=D(i,1);
    data_label=find(xt_x(:,data_point)==1);
    time_window=D(i,2);
    %data_cost=data_cost-log_2(x_xt(data_point,data_label)*y_xt(time_window,data_label)*xt(data_label));
    data_cost=data_cost-log_2(y_xt(time_window,data_label)*x(data_point));
end
total_cost=model_cost+data_cost;
end
