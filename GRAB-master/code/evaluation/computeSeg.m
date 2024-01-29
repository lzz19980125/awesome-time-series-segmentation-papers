% function:
% Compute segmentation error
% 
% Input parameter:
% GT - ground truth
% splitLoc - the change points
%
% Output parameter
% error - the segmentation error
%% compute segmentation error
function [error] = computeSeg(GT,splitLoc)
splitLoc = sort(splitLoc);
idx = 1;
for i = 1:length(GT)-1
    if GT(i) ~= GT(i+1)
        gt(idx) = i;
        idx = idx+1;
    end
end

S = zeros(length(gt),1);
F = zeros(length(gt),1);
for i = 1:length(gt)
    if length(gt) == 1
        S(1) = 1;
        F(1) = length(GT);
        break;
    end
    if i == 1
        start = 1;
        finish = gt(i)+ceil((gt(i+1)-gt(i)-1)/2);
    elseif i == length(gt)
        start = gt(i)-floor((gt(i)-gt(i-1)-1)/2);
        finish = length(GT);
    else
        start = gt(i)-floor((gt(i)-gt(i-1)-1)/2);
        finish = gt(i)+ceil((gt(i+1)-gt(i)-1)/2);
    end
    S(i) = start;
    F(i) = finish;
end

error = 0;
for i = 1:length(gt)
    val = gt(i);
    start = S(i);
    finish = F(i);
    flag = 0;
    for j = 1:length(splitLoc)
        predict = splitLoc(j);
        if start <= predict && predict <= finish
            %se_p
            error = error+abs(predict-val);
            flag = 1;
        end
    end
    if flag == 0
        % se_r
        error = error+finish-start+1;
    end
end
error = error/length(GT);