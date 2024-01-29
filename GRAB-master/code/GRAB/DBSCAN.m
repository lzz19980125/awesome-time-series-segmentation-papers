% function:
% Return the density-based clustering result
%
% Input parameters:
% D - data points
% Eps - the parameter to determine the neighborhood of a point
% MinPts - the parameter to identify core points
%
% Output parameters:
% k - number of clusters
% C - cluster identifiers
function [k, C] = DBSCAN(D, Eps, MinPts)
%% Initialization
dLength = length(D);
O = zeros(dLength, 1);
C = cell(dLength);
d = zeros(dLength, dLength);

%% identify core points
for i = 1:dLength
    for j = dLength:-1:i+1
        if abs(D(i)-D(j)) <= Eps
            d(i, j) = 1;
            d(j, i) = d(i, j);
        end
    end
    if sum(d(i, :)) >= MinPts
            O(i) = 1;
    end
end

%% cluster data points by DBSCAN
k = 0;
Tau = ones(dLength, 1);
while sum(O) ~= 0
    Tau_old = Tau;
    j = 1;
    while O(j) == 0
        j = j+1;
    end
    Q = zeros(dLength, 1);
    Q(j) = 1;
    Tau(j) = 0;
    while sum(Q) ~= 0
        m = 1;
        while Q(m) == 0
            m = m+1;
        end
        Q(m) = 0;
        if sum(d(m,:)) >= MinPts
            for l = 1:dLength
                if d(m, l) == 1 && Tau(l) ~= 0
                    Q(l) = 1;
                    Tau(l) = 0;
                end
            end
        end
    end
    k = k+1;
    for i = 1:dLength
        if Tau(i) ~= 0
            Tau_old(i) = 0;
        end
        if Tau_old(i) ~= 0
            O(i) = 0;
        end
    end
    C{k} = Tau_old;
end
 
%% assign cluster identifiers
for i = 1:k
    C{i} = find(C{i}~=0);
end