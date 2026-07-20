
function [tau, dOut, dBar] = calcDampedStepSizeAndDirectionV0(myU, d, N_d_lin, f, fPrime, dom, newtonIter, tauOld, dOld, dBarOld)
    tau = tauOld;
    epsilon = 1e-6;
    tauMin = epsilon;
    endSearch = false;
    cFactor=1.0;
    dOut = d;
    initPrediction = 1;
    fprintf('calc damped newton step\n');
    fprintf('\t%8s   %8s   %8s\n','||d||','cFactor','tau');
    while (~ endSearch)


        if(newtonIter>1 && initPrediction)
            mu = norm(dOld)*norm(dBarOld)/(norm(dBarOld-d, 'fro')*norm(d))*tau;
            tau=min(1,mu);
            initPrediction = 0;
        end

        if(tau<tauMin)
            % tau too small, try full newton to try a new base point
            tau=1.0;
            dOut = d;
            dBar = d;
            resTrial = f(myU + tau*d);
            break;
        end

        uTrial = myU + tau*d;
        resTrial = f(uTr¸ial);



        N_temp = chebop(dom);
        N_temp.op = @(x,d) fPrime(myU,d);
        N_temp.init = chebfun(@(x) 0, dom);
        N_temp.bc = @(x,d) [d(-2); d(0)];
        %N_temp.lbc = @(d,x) d;
        %N_temp.rbc = @(d,x) d;
        dBar = N_temp\resTrial;

        %dBar = N_d_lin\0;
        dBar = -dBar;


        % contraction factor
        cFactor = norm(dBar)/norm(d);

        % correction factor
        muPrime = 1.0/2.0*norm(d)*tau^2/(norm(dBar - (1-tau)*d,'fro'));



        if(cFactor>1)
            tau = min(muPrime,1.0/2.0*tau);
            continue
        end

        tauPrime = min(muPrime,1.0);
        if(abs(tauPrime-1.0)<epsilon)
            tau = 1.0;
            %dOut = dBar;
            break;
        else
            if(tauPrime>4*tau)
                tau = tauPrime;
            else
                break
            end
        end
        
        fprintf('\t%8.4f   %8.4f   %8.4f\n', norm(dOut), cFactor, tau);
    end

    %norm(d)
    %norm(dBar)
    %norm(myU+tau*dOut)
    %save(strcat("myDampVars",string(newtonIter)), "myU", "d", "dBar", "fPrime", "resTrial")
    fprintf('\t%8.4f   %8.4f   %8.4f\n', norm(dOut), cFactor, tau);
end