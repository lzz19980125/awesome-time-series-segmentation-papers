% function:
% Find natural structures in time series
% 
% Input parameters:
% A - the time series
% SubsequenceLength - the pilot subsequence length
% topN - the number of most similar subsequences to a pilot
%
% Output parameters:
% newCommunity - segments and their state memberships
function [newCommunity] = GRAB(A,SubsequenceLength,topN)
%% Check input
if length(A) == size(A, 2)
   A = A';
end

%% Initialization
MatrixProfileLength = length(A)-SubsequenceLength+1;
topN = round(MatrixProfileLength/topN);
Eps = round(MatrixProfileLength/60);
MinPts = 20;
numPilots = 100;

[X, n, sumx2, sumx, meanx, sigmax2, sigmax] = ...
    fastfindNNPre(A, SubsequenceLength);

segments = cell(numPilots, 1);
prob = ones(MatrixProfileLength, 1);

%% Local fragment generation
i = 0;
cntGoodPilots = 0;
cntPoorPilots = 0;
while cntGoodPilots ~= numPilots
    % determine the pilot in the current round
    i = i + 1;
    cumProb = cumsum(prob);   
    idx = rand(1,1) * cumProb(length(cumProb));
    if (idx <= cumProb(1))
        idx = 1;
    else
        for j = 2:MatrixProfileLength
            if (cumProb(j-1) < idx && idx <= cumProb(j))
                idx = j;
                break;
            end
        end
    end
    
    % compute the z-normalized Euclidean distances
    % The following four lines are the code provided in the following URL
    % https://sites.google.com/site/onlinesemanticsegmentation/
    subsequence = A(idx:idx+SubsequenceLength-1);
    distanceProfile = fastfindNN(X, subsequence, n, SubsequenceLength, ...
        sumx2, sumx, meanx, sigmax2, sigmax);
    distanceProfile = abs(distanceProfile);
   
    % cluster topN subsequences by adjacency
    [~, index] = sort(distanceProfile);
    [k, C] = DBSCAN(index(1:topN), Eps, MinPts);
    
    prob(idx) = 0;
    
    numClustered = 0;
    for j = 1:k
        numClustered = numClustered + length(C{j});
    end
        
    % pilot filtering
    if numClustered >= topN*0.7
        cntGoodPilots = cntGoodPilots + 1; 
        clear recordCovered;      
        [segments, flags] = updateSegments(cntGoodPilots, k, C, index, segments, SubsequenceLength, MatrixProfileLength);
        covered = find(flags==1);
        prob(covered) = prob(covered)/2; 
        cntPoorPilots = 0;
    else
        cntPoorPilots = cntPoorPilots + 1;
    end
    
    % parameter adjustment
    if (cntPoorPilots >= 20 || ((i == MatrixProfileLength) && (cntGoodPilots ~= numPilots)))
        MinPts = MinPts - 2;
        i = 0;
        cntGoodPilots = 0;
        cntPoorPilots = 0;
        clear segments;
        clear prob;
        segments = cell(numPilots, 1);
        prob = ones(MatrixProfileLength, 1);
    end
end

%% Global fragment generation
[newSegments, weights, missing] = processSegments(segments,length(A));

%% Graph partition
[COMTY,~] = cluster_jl(weights,1,1,0,0);
level = length(COMTY.COM);
community = cell(length(COMTY.SIZE{level}), 1);
counters = zeros(length(COMTY.SIZE{level}), 1);
for i = 1:length(COMTY.COM{level})
    id = COMTY.COM{level}(i);
    a = newSegments(i*2-1);
    b = newSegments(i*2);
    counters(id, 1) = counters(id, 1)+1;
    community{id}{1}{counters(id, 1)}(1) = a;
    counters(id, 1) = counters(id, 1)+1;
    community{id}{1}{counters(id, 1)}(1) = b;
end

%% Post-processing
community = processMissing(community, missing, length(A), SubsequenceLength, A, topN, X, n, sumx2, sumx, meanx, sigmax2, sigmax);

% merge adjacent subsequences from the same state
for i = 1:length(community)
    for j = length(community{i}{1})-1:-2:3
        x1 = community{i}{1}{j}(1);
        y1 = community{i}{1}{j+1}(1);
        for k = j-2:-2:1
            x2 = community{i}{1}{k}(1);
            y2 = community{i}{1}{k+1}(1);
            if x1 - y2 == 1
                community{i}{1}{k+1}(1) = y1;
                community{i}{1}(j+1) = [];
                community{i}{1}(j) = [];
                break;
            end
            if x2 - y1 == 1
                community{i}{1}{k}(1) = x1;
                community{i}{1}(j+1) = [];
                community{i}{1}(j) = [];
                break;
            end
        end
    end
end

% rearrange state membership
for i = 1:length(community)
    minVal = inf;
    for j = 1:2:length(community{i}{1})
        x = community{i}{1}{j}(1);
        if x < minVal
            minVal = x;
        end
    end
    minStart(i) = minVal; 
end
[~,index] = sort(minStart);
for i = 1:length(community)
    newCommunity{i}(1) = community{index(i)}(1);
end

%% Local fragment generation with a single pilot
function [segments, flags] = updateSegments(round, k, C, index, segments, SubsequenceLength, MatrixProfileLength)
RS = cell(2*k,1);
flags = zeros(MatrixProfileLength, 1);
for i = 1:k
    minVal = min(index(C{i}));
    maxVal = max(index(C{i}))+SubsequenceLength-1;
    if maxVal > MatrixProfileLength
        maxVal = MatrixProfileLength;
    end   
    flags(minVal:maxVal) = 1;
    RS{2*i-1}(1) = minVal;
    RS{2*i}(1) = maxVal;
end
segments{round}{1} = RS;

%% Global fragment generation
function [newSegments, weights, missing] = processSegments(segments,totalLength)
% ePoints: endpoints of local fragments
counter = 1;
ePoints(counter) = 1;
for i = 1:length(segments)
    for j = 1:length(segments{i}{1})
        counter = counter + 1;
        ePoints(counter) = segments{i}{1}{j}(1);
    end
end
counter = counter + 1;
ePoints(counter) = totalLength+1;
ePoints = sort(ePoints);

% get raw global fragments by endpoints
counter = 0;
pre = ePoints(1);
for i = 2:length(ePoints)
    cur = ePoints(i); 
    if pre ~= cur
        counter = counter + 1;
        newSegments(counter) = pre;
        counter = counter + 1;
        newSegments(counter) = cur-1;
        pre = cur;
    end
end

% get missing subsequences
covered = zeros(totalLength, 1);
for i = 1:length(segments)
    for j = 1:2:length(segments{i}{1})
        x = segments{i}{1}{j}(1);
        y = segments{i}{1}{j+1}(1);
        covered(x:y) = 1;
    end
end
counter = 0;
for i = length(newSegments)-1:-2:1
    a = newSegments(i);
    b = newSegments(i+1);    
    if (b-a) >= 2
        numOfZeros = length(find(covered(a+1:b-1)==0));
        if numOfZeros ~= 0
            counter = counter + 1;
            missing(counter) = a;
            counter = counter + 1;
            missing(counter) = b;
            newSegments(i+1) = [];
            newSegments(i) = [];
        end
    end
end

% compute occurrence of every global fragment
occurrence = zeros(round(length(newSegments)/2),length(segments));
for i = 1:2:length(newSegments)-1
    x1 = newSegments(i);
    y1 = newSegments(i+1);
    for j = 1:length(segments)
        for k = 1:2:length(segments{j}{1})
            x = segments{j}{1}{k}(1);
            y = segments{j}{1}{k+1}(1);
            if (x-2 <= x1 && y+1 >= y1) || (x-1 <= x1 && y+2 >= y1)
                occurrence(round(i/2),j) = 1;
                break;
            end
        end
    end
end
   
% merge raw global fragments to avoid short fragments
targetNum = 50;
curNum = round(length(newSegments)/2);
coOccurNum = zeros(curNum-1,1);
mergeNum = zeros(curNum,1);
for i = 1:curNum-1
    coOccurNum(i) = sum(occurrence(i,:) & occurrence(i+1,:));
    x2 = newSegments(2*i);
    y1 = newSegments(2*i+1);
    if y1-x2 ~= 1
        coOccurNum(i) = -1;
    end    
end

while (curNum > targetNum)
    [~,index] = max(coOccurNum);
    newSegments(2*index) = newSegments(2*index+2);
    newSegments(2*index+2) = [];
    newSegments(2*index+1) = [];
    occurrence(index,:) = occurrence(index,:)|occurrence(index+1,:);
    occurrence(index+1,:) = [];
    mergeNum(index) = mergeNum(index)+mergeNum(index+1)+1;
    mergeNum(index+1) = [];
    curNum = curNum-1;
    
    if (index~=1)
        num = mergeNum(index-1)+mergeNum(index)+1;
        coOccurNum(index-1) = sum(occurrence(index-1,:)&occurrence(index,:))/num;
        x2 = newSegments(2*index-2);
        y1 = newSegments(2*index-1);
        if y1-x2 ~= 1
            coOccurNum(index-1) = -1;
        end 
    end
    
    if (index<length(coOccurNum))
        num = mergeNum(index)+mergeNum(index+1)+1;
        coOccurNum(index) = sum(occurrence(index,:)&occurrence(index+1,:))/num;
        x2 = newSegments(2*index);
        y1 = newSegments(2*index+1);
        if y1-x2 ~= 1
            coOccurNum(index) = -1;
        end 
        coOccurNum(index+1) = [];
    else
        coOccurNum(index) = [];
    end
end 
    
% compute the weights of the edges
weights = zeros(round(length(newSegments)/2),round(length(newSegments)/2));
for i = 1:2:length(newSegments)-3
    for j = i+2:2:length(newSegments)-1
        a = round(i/2);
        b = round(j/2);
        occurrenceA = sum(occurrence(a,:));
        occurrenceB = sum(occurrence(b,:));
        occurrenceC = sum(occurrence(a,:)&occurrence(b,:));
        weight = occurrenceC/(min(occurrenceA,occurrenceB));
        weights(a,b) = weight;
        weights(b,a) = weight;
    end
end   
        
%% Process missing subsequences
function [segments] = processMissing(segments, missing, totalLength, SubsequenceLength, A, topN, X, n, sumx2, sumx, meanx, sigmax2, sigmax)
for i = 1:2:length(missing)-1
    x = missing(i);
    y = missing(i+1);
    if y-x <= totalLength*0.05
        left = 0;
        right = 0;
        for j = 1:length(segments)
            for k = 1:2:length(segments{j}{1})
                if (segments{j}{1}{k+1}(1)+1 == x)
                    left = j;
                end
                if (segments{j}{1}{k}(1)-1 == y)
                    right = j;
                end
            end
        end
        if left == 0
            val = length(segments{right}{1});
            segments{right}{1}{val+1} = x;
            segments{right}{1}{val+2} = y;
        elseif right == 0
            val = length(segments{left}{1});
            segments{left}{1}{val+1} = x;
            segments{left}{1}{val+2} = y;
        elseif left == right
            val = length(segments{right}{1});
            segments{right}{1}{val+1} = x;
            segments{right}{1}{val+2} = y;
        else
            idx = round((x+y)/2);
            if ((idx+SubsequenceLength-1)>totalLength)
                idx = totalLength-SubsequenceLength+1;
            end
            subsequence = A(idx:idx+SubsequenceLength-1);
            distanceProfile = fastfindNN(X, subsequence, n, SubsequenceLength, ...
                sumx2, sumx, meanx, sigmax2, sigmax);
            distanceProfile = abs(distanceProfile);
            [~,index] = sort(distanceProfile);
            totalNum = 0;
            for k = 1:2:length(segments{left}{1})
                start = segments{left}{1}{k}(1);
                finish = segments{left}{1}{k+1}(1);
                for l = 1:topN
                    if index(l) >= start && index(l) <= finish
                        totalNum = totalNum + 1;
                    end
                end
            end
            max = totalNum;
            maxI = left;
            totalNum = 0;
            for k = 1:2:length(segments{right}{1})
                start = segments{right}{1}{k}(1);
                finish = segments{right}{1}{k+1}(1);
                for l = 1:topN
                    if index(l) >= start && index(l) <= finish
                        totalNum = totalNum + 1;
                    end
                end
            end
            if totalNum > max
                maxI = right;
            end
            val = length(segments{maxI}{1});
            segments{maxI}{1}{val+1} = x;
            segments{maxI}{1}{val+2} = y;
        end
    else
        idx = round((x+y)/2);
        if ((idx+SubsequenceLength-1) > totalLength)
            idx = totalLength-SubsequenceLength+1;
        end
        subsequence = A(idx:idx+SubsequenceLength-1);
        distanceProfile = fastfindNN(X, subsequence, n, SubsequenceLength, ...
            sumx2, sumx, meanx, sigmax2, sigmax);
        distanceProfile = abs(distanceProfile);
        [~,index] = sort(distanceProfile);
        max = 0;
        maxI = 0;
        for j = 1:length(segments)
            totalNum = 0;
            for k = 1:2:length(segments{j}{1})
                start = segments{j}{1}{k}(1);
                finish = segments{j}{1}{k+1}(1);
                for l = 1:topN
                    if index(l) >= start && index(l) <= finish
                        totalNum = totalNum + 1;
                    end
                end
            end
            if totalNum >= max
                max = totalNum;
                maxI = j;
            end
        end
        val = length(segments{maxI}{1});
        segments{maxI}{1}{val+1} = x;
        segments{maxI}{1}{val+2} = y;
    end
end
        
%% The following two functions are the code provided in the following URL
% https://sites.google.com/site/onlinesemanticsegmentation/
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