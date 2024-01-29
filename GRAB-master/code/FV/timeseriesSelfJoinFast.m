% The following codes are provided in the following URL
% https://sites.google.com/site/onlinesemanticsegmentation/

% [MatrixProfile, MatrixProfileIndex] = Time_series_Self_Join_Fast(A, SubsequenceLength)
% Output:
%     MatrixProfile: matrix porfile of the self-join (vector)
%     MatrixProfileIndex: matrix porfile index of the self-join (vector)
% Input:
%     A: input time series (vector)
%     SubsequenceLength: interested subsequence length (scalar)
%

function [MatrixProfile, MPindex] = timeseriesSelfJoinFast(A, SubsequenceLength)
%% set trivial match exclusion zone
exclusionZone = round(SubsequenceLength/4);

%% check input
if SubsequenceLength > length(A)/2
    error('Error: Time series is too short relative to desired subsequence length');
end
if SubsequenceLength < 4
    error('Error: Subsequence length must be at least 4');
end
if length(A) == size(A, 2)
   A = A'; 
end

%% initialization
MatrixProfileLength = length(A) - SubsequenceLength + 1;
MatrixProfile = zeros(MatrixProfileLength, 1);
MPindex = zeros(MatrixProfileLength, 1);
[X, n, sumx2, sumx, meanx, sigmax2, sigmax] = ...
    fastfindNNPre(A, SubsequenceLength);

%% compute the matrix profile
pickedIdx = randperm(MatrixProfileLength);
for i = 1:MatrixProfileLength
    % compute the distance profile
    idx = pickedIdx(i);
    subsequence = A(idx:idx+SubsequenceLength-1);
    distanceProfile = fastfindNN(X, subsequence, n, SubsequenceLength, ...
        sumx2, sumx, meanx, sigmax2, sigmax);
    distanceProfile = abs(distanceProfile);
    
    % apply exclusion zone
    exclusionZoneStart = max(1, idx-exclusionZone);
    exclusionZoneEnd = min(MatrixProfileLength, idx+exclusionZone);
    distanceProfile(exclusionZoneStart:exclusionZoneEnd) = inf;
    
    % figure out and store the neareest neighbor
    if i == 1
        MatrixProfile = distanceProfile;
        MPindex(:) = idx;
        [MatrixProfile(idx), MPindex(idx)] = min(distanceProfile);
    else
        updatePos = distanceProfile < MatrixProfile;
        MPindex(updatePos) = idx;
        MatrixProfile(updatePos) = distanceProfile(updatePos);
        [MatrixProfile(idx), MPindex(idx)] = min(distanceProfile);
    end
end

%% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [X, n, sumx2, sumx, meanx, sigmax2, sigmax] = fastfindNNPre(x, m)
n = length(x);
x(n+1:2*n) = 0;
X = fft(x);
cum_sumx = cumsum(x);
cum_sumx2 =  cumsum(x.^2);
sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];
sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];
meanx = sumx./m;
sigmax2 = (sumx2./m)-(meanx.^2);
sigmax = sqrt(sigmax2);

function dist = fastfindNN(X, y, n, m, sumx2, sumx, meanx, sigmax2, sigmax)
%x is the data, y is the query
y = (y-mean(y))./std(y,1);                      %Normalize the query
y = y(end:-1:1);                                %Reverse the query
y(m+1:2*n) = 0;

%The main trick of getting dot products in O(n log n) time
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

%compute y stats -- O(n)
sumy = sum(y);
sumy2 = sum(y.^2);

%computing the distances -- O(n) time
dist = (sumx2 - 2*sumx.*meanx + m*(meanx.^2))./sigmax2 - 2*(z(m:n) - sumy.*meanx)./sigmax + sumy2;
%dist = 1-dist./(2*m);
dist = sqrt(dist);
