% function:
% Compute clustering accuracy
% 
% Input parameter:
% GT - ground truth
% community - segments and their state memberships
%
% Output parameter
% accuracy - the clustering accuracy
function [accuracy] = computeClu(GT,community)
for i = 1:length(community)
    for j = 1:2:length(community{i}{1})
        for k = community{i}{1}{j}(1):1:community{i}{1}{j+1}(1)
            c1(k) = i;
        end
    end
end
c2 = GT;
accuracy = computeARI(c1, c2);
