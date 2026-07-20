function [marker, markerSize] = getMarkersAndSize(iter)
    markerSize = 10;
    MarkerList = {'+','o','x','v','d','^','s','>','<'};
    marker = MarkerList(mod(iter-1,length(MarkerList))+1);
end