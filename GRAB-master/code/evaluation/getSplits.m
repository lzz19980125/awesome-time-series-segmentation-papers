% function:
% Get change points between different states
%
% Input parameter:
% TS - the time series
% community - the segments and segments and their state memberships
% 
% Output parameter
% splitLoc - the change points
function [splitLoc] = getSplits(TS,community)
idx = 1;
isNull = 1;
for i = 1:length(community)
    for j = 1:length(community{i}{1})
        val = community{i}{1}{j}(1);
        if val ~= 1 && val ~= length(TS)
            flag = 0;
            if isNull == 0
                for k = 1:length(splitLoc)
                    if splitLoc(k)-val == 1
                        splitLoc(k) = val;
                        flag = 1;
                        break;
                    elseif splitLoc(k)-val == -1
                        flag = 1;
                        break;
                    end
                end
            end
            isNull = 0;
            if flag == 0
                splitLoc(idx) = val;
                idx = idx+1;
            end
        end
    end
end