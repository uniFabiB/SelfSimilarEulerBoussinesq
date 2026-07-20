

function plotFunctionsAndDiff(varargin)
    myLegend = strings(nargin+ (nargin*(nargin-1))/2,1);


    title = inputname(1);
    if nargin > 1
        for i = 2:nargin
            title = strcat(title , ", ", inputname(i)); % does not work if argument is not a single variable like a+b
        end
    end

    figure("Name",title);
    for i = 1:nargin
        if(i>0)
            hold on;
        end
        myLegend(i) = inputname(i); % does not work if argument is not a single variable like a+b
        f = varargin{i};
        semilogy(abs(f),'x-');
    end

    legendIndex = nargin;



    for i = 1:nargin
        for j = i+1:nargin
            f = varargin{i};
            g = varargin{j};
            fName = inputname(i);
            gName = inputname(j);
            h=abs(f-g);
            %[Linfty, posMax]= max(h{1});
            [Linfty, posMax]= max(h);
            legendIndex = legendIndex + 1;
            myLegend(legendIndex) = fName + " - " + gName + " (max=" + Linfty + " at x="+ posMax + ")";
            semilogy(h+1e-19,'o-');
        end
    end

    legend(myLegend, 'Location', 'southeast');
    hold off;
end