% The following codes are provided in the following URL
% https://sites.google.com/site/onlinesemanticsegmentation/

% Segments the time series into several parts based on repeatability of a 
% pattern within each of the regions.
% 
% function [crosscount, splitLoc] = SegmentTimeSeries(slWindow, MPindex)
%
% Input parameters:
% slWindow - sliding window size, used to reduce spurious segmentation of
%           the regions which lengths are less than the sliding window size
% MPindex - matrix profile index for the time series to segment
%
% Output parameters:
% crosscount - number of crossings at each point
% splitLoc - split locations
%
function [crosscount] = segmentTimeSeries(MPindex)
profile_len = length(MPindex);
% find the number of arcs crossing over each index i
nnmark=zeros(1,profile_len);
for i=1:profile_len
    small=min(i,MPindex(i));
    large=max(i,MPindex(i));
    nnmark(small)=nnmark(small)+1;
    nnmark(large)=nnmark(large)-1;
end
% iterate over nnmark and cumulatively sum its values for each index i
crosscount = cumsum(nnmark); % number of arcs crossing over each index
end
   
   
