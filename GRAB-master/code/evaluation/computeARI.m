% function:
% compute the adjusted rand index between two clustering results
function ri = computeARI(c1, c2)
%% Check input
if length(c1)~=numel(c1)
	c1 = c1(:);
end
if length(c2)~=numel(c2)
    c2 = c2(:);
end

%% Preliminary computations
N = length(c1);
[~, ~, c1] = unique(c1);
N1 = max(c1);
[~, ~, c2] = unique(c2);
N2 = max(c2);

for i=1:1:N1
    for j=1:1:N2
        G1 = find(c1==i);
        G2 = find(c2==j);
        n(i,j) = length(intersect(G1,G2));
    end
end

%% Calculate the adjusted rand index
ssm = 0;
sm1 = 0;
sm2 = 0;
for i=1:1:N1
    for j=1:1:N2
        ssm = ssm + nchoosek2(n(i,j),2);
    end
end
temp = sum(n,2);
for i=1:1:N1
    sm1 = sm1 + nchoosek2(temp(i),2);
end
temp = sum(n,1);
for i=1:1:N2
    sm2 = sm2 + nchoosek2(temp(i),2);
end
NN = ssm - sm1*sm2/nchoosek2(N,2);
DD = (sm1 + sm2)/2 - sm1*sm2/nchoosek2(N,2);
ri = NN/DD;

%% Definition of n choose k
function c = nchoosek2(a,b)
if a>1
    c = nchoosek(a,b);
else
    c = 0;
end