function plotFunctionsAndDiff(varargin)
    lineWidth = 1;
    [~, markerSize] = fb.getMarkersAndSize(0);

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
        semilogy(abs(f), 'x-', 'MarkerSize', markerSize, 'LineWidth', lineWidth);
    end

    legendIndex = nargin;



    for i = 1:nargin
        for j = i+1:nargin
            f = varargin{i};
            g = varargin{j};
            fName = inputname(i);
            gName = inputname(j);

            if domain(f) == domain(g)
                h_temp=f-g;
            else
                dom_f = domain(f);
                dom_g = domain(g);
                dom_inter = [max(dom_f(1), dom_g(1)), min(dom_f(2), dom_g(2))];
                dom_union = [min(dom_f(1), dom_g(1)), max(dom_f(2), dom_g(2))];
                h_temp = chebfun(@(x) (x>dom_inter(1))*(x<dom_inter(2))*(f(x) - g(x)) + 0*x, dom_union, 'splitting', 'on');
            end
            h=abs(h_temp);
            
            %[Linfty, posMax]= max(h{1});
            [Linfty, posMax]= max(h);
            legendIndex = legendIndex + 1;
            myLegend(legendIndex) = fName + " - " + gName + " (max=" + Linfty + " at x="+ posMax + ")";
            semilogy(h+1e-19,'o-', 'MarkerSize', markerSize, 'LineWidth', lineWidth);
        end
    end

    legend(myLegend, 'Location', 'southeast');
    hold off;
end