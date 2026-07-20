
%% ===================== Safe getter (zero outside -N..N) =====================
function f = getVar(Xc, k, Ndom, dom)
    if k < -Ndom || k > Ndom
        f = chebfun(0, dom);
    else
        f = Xc{k + Ndom + 1};
    end
end