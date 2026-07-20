
function plot(varargin)
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
        Markers = {'+','o','*','x','v','d','^','s','>','<'};
        f = varargin{i};
        marker = Markers(mod(i-1,length(Markers))+1);
        lineStyleString = "-" + marker;
        plot(f,lineStyleString);
    end

    legend(myLegend, 'Location', 'southeast');
    hold off;
end