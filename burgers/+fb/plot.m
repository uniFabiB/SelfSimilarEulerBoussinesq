
function plot(varargin)
    lineWidth = 1;
    myLegend = strings(nargin,1);


    title = inputname(1);
    if nargin > 1
        for i = 2:nargin
            title = strcat(title , ", ", inputname(i)); % does not work if argument is not a single variable like a+b
        end
    end

    figure("Name",title);
    for i = 1:nargin
        if(i>1)
            hold on;
        end
        myLegend(i) = inputname(i); % does not work if argument is not a single variable like a+b
        f = varargin{i};
        [marker, markerSize] = fb.getMarkersAndSize(i);
        lineStyleString = marker+"-";
        plot(f, lineStyleString, 'MarkerSize', markerSize, 'LineWidth', lineWidth);
    end

    legend(myLegend, 'Location', 'southeast');
    hold off;
end