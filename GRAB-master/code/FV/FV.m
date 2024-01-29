% function:
% Find natural structures in time series based on the FLUSS algorithm
%
% Input parameter:
% ts - the time series
% SubsequenceLength - the length of subsequences
% zone - the exclusion zone size to avoid trivial minimum
% Output parameter:
% newCommunity - segments and their state memberships
function [newCommunity] = FV(ts, SubsequenceLength, zone)
%% get arc curve
[~, MPindex] = timeseriesSelfJoinFast(ts, SubsequenceLength);
[crosscount] = segmentTimeSeries(MPindex);
[crosscount] = normCrossCountAll(crosscount,SubsequenceLength);      
splitLoc = [];

%% find all valley points
while true
    [minVal, minIdx] = min(crosscount);
    if minVal == inf
        break;
    else
        splitLoc = [splitLoc, minIdx];
        exu_st = max(1, minIdx-zone);
        exu_ed = min(length(crosscount), minIdx+zone);
        crosscount(exu_st:exu_ed) = inf;
    end
end

%% get segments by valley points
splitLoc = sort(splitLoc);
cnt = 1;
for i = 1:length(splitLoc)
    if i == 1
        segments(cnt) = 1;
        cnt = cnt+1;
        segments(cnt) = splitLoc(i);
        cnt = cnt+1;
    else
        segments(cnt) = splitLoc(i-1)+1;
        cnt = cnt+1;
        segments(cnt) = splitLoc(i);
        cnt = cnt+1;
    end
end
segments(cnt) = splitLoc(i)+1;
cnt = cnt+1;
segments(cnt) = length(ts);

%% compute the weights of the edges
weights = zeros(length(splitLoc)+1,length(splitLoc)+1);
bound = length(ts)-SubsequenceLength+1;
for i = 1:2:length(segments)
    for j = 1:2:length(segments)
            x1 = segments(i);
            y1 = segments(i+1);
            x2 = segments(j);
            y2 = segments(j+1);
            if y1 > bound
                y1 = bound;
            end
            if y2 > bound
                y2 = bound;
            end
            if x1 < bound && x2 < bound
                weights(round(i/2),round(j/2)) = length(find(MPindex(x1:y1)>=x2 & MPindex(x1:y1)<=y2));
                weights(round(i/2),round(j/2)) = weights(round(i/2),round(j/2))/length(MPindex);
            end
    end
end

%% graph partition
[COMTY,~] = cluster_jl(weights,1,1,0,0);
level = length(COMTY.COM);
community = cell(length(COMTY.SIZE{level}), 1);
counters = zeros(length(COMTY.SIZE{level}), 1);
for i = 1:length(COMTY.COM{level})
    id = COMTY.COM{level}(i);
    a = segments(i*2-1);
    b = segments(i*2);
    counters(id, 1) = counters(id, 1)+1;
    community{id}{1}{counters(id, 1)}(1) = a;
    counters(id, 1) = counters(id, 1)+1;
    community{id}{1}{counters(id, 1)}(1) = b;
end

%% post-processing
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