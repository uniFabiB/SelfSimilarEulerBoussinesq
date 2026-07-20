
%% ===================== Main residual builder =====================
function out = buildSystem(r, allVars, Ndom, lambda, dom)
    nModes = 2*Ndom + 1;
    Ac = allVars(1:nModes);              % A_{-N}...A_{N}
    Bc = allVars(nModes+1:2*nModes);     % B_{-N}...B_{N}
    Cc = allVars(2*nModes+1:3*nModes);   % C_{-N}...C_{N}

    idxOf = @(k) k + Ndom + 1;
    invr  = 1./r;
    invr2 = 1./r.^2;

    eqA  = cell(nModes,1);
    eqB  = cell(nModes,1);
    eqCa = cell(nModes,1);
    eqCb = cell(nModes,1);

    for m = -Ndom:Ndom
        %% ---- Equation 0: A_m  (A_m + r A_m' = 0) ----
        Am  = getVar(Ac, m, Ndom, dom);
        dAm = diff(Am);
        eqA{idxOf(m)} = Am + r.*dAm;
        
        %% ---- Equation 1: B_m ----
        Bm  = getVar(Bc, m,  Ndom, dom);
        dBm = diff(Bm);  d2Bm = diff(Bm,2);

        term1 = invr.*Bm + dBm;
        term2 = (lambda+1) * ( -invr.*Bm + dBm + r.*d2Bm );

        sumB = chebfun(0, dom);
        for k = -Ndom:Ndom
            Ak   = getVar(Ac, k, Ndom, dom);
            Bmk  = getVar(Bc, m-k, Ndom, dom);
            dBmk = diff(Bmk); d2Bmk = diff(Bmk,2);
            sumB = sumB + Ak.*( -invr2.*Bmk + invr.*dBmk + d2Bmk );
        end

        Cm1 = getVar(Cc, m-1, Ndom, dom);
        Cp1 = getVar(Cc, m+1, Ndom, dom);
        rhsB = 0.5*( diff(Cm1) + diff(Cp1) );

        eqB{idxOf(m)} = term1 + term2 + sumB - rhsB;

        %% ---- Equation 2 (eqC_a) ----
        dCm1 = diff(Cm1); dCp1 = diff(Cp1);
        ddCm1 = diff(Cm1,2); ddCp1 = diff(Cp1,2);

        LHS = dCm1 + dCp1;

        sum1 = chebfun(0, dom); sum2 = chebfun(0, dom);
        sum3 = chebfun(0, dom); sum4 = chebfun(0, dom);
        rSum1 = chebfun(0, dom); rSum2 = chebfun(0, dom); rSum3 = chebfun(0, dom);

        for k = -Ndom:Ndom
            Ak  = getVar(Ac, k, Ndom, dom); dAk = diff(Ak);
            Bk  = getVar(Bc, k, Ndom, dom); dBk = diff(Bk);

            Cm_k_3  = getVar(Cc, m-k-3, Ndom, dom); dCm_k_3  = diff(Cm_k_3);
            Cm_k_1p = getVar(Cc, m-k+1, Ndom, dom); dCm_k_1p = diff(Cm_k_1p);
            Cm_k_1m = getVar(Cc, m-k-1, Ndom, dom); dCm_k_1m = diff(Cm_k_1m);
            Cm_k_3p = getVar(Cc, m-k+3, Ndom, dom); dCm_k_3p = diff(Cm_k_3p);
            ddCm_k_1m = diff(Cm_k_1m,2); ddCm_k_1p = diff(Cm_k_1p,2);

            sum1 = sum1 + 0.25*dAk.*( dCm_k_3 + dCm_k_1p + dCm_k_1m + dCm_k_3p );
            sum2 = sum2 + (1/(8*1i))*( -dBk + invr.*Bk ).*( dCm_k_3 - dCm_k_1p + dCm_k_1m - dCm_k_3p );
            sum3 = sum3 + 0.5*Ak.*( ddCm_k_1m + ddCm_k_1p );
            sum4 = sum4 - (1/(2*1i))*( invr.*Bk ).*( dCm_k_1m - dCm_k_1p );

            rSum1 = rSum1 + 0.25*dAk.*( dCm_k_3 - dCm_k_1p - dCm_k_1m + dCm_k_3p );
            rSum2 = rSum2 - (1/(4*1i))*( dBk + invr.*Bk ).*( dCm_k_1m - dCm_k_1p );
            rSum3 = rSum3 + (1/(8*1i))*( invr.*Bk - dBk ).*( dCm_k_3 + dCm_k_1p - dCm_k_1m - dCm_k_3p );
        end

        LHS = LHS + sum1 + sum2 + 0.5*(lambda+1)*r.*( ddCm1 + ddCp1 ) + sum3 + sum4;
        RHS = rSum1 + rSum2 + rSum3;

        eqCa{idxOf(m)} = LHS - RHS;

        %% ---- Equation 3 (eqC_b) ----
        LHSb = (1/1i)*( dCm1 - dCp1 );

        s1 = chebfun(0, dom); s2 = chebfun(0, dom);
        s3 = chebfun(0, dom); s4 = chebfun(0, dom);
        r1 = chebfun(0, dom); r2 = chebfun(0, dom); r3 = chebfun(0, dom);

        for k = -Ndom:Ndom
            Ak  = getVar(Ac, k, Ndom, dom); dAk = diff(Ak);
            Bk  = getVar(Bc, k, Ndom, dom); dBk = diff(Bk);

            Cm_k_3  = getVar(Cc, m-k-3, Ndom, dom); dCm_k_3  = diff(Cm_k_3);
            Cm_k_1p = getVar(Cc, m-k+1, Ndom, dom); dCm_k_1p = diff(Cm_k_1p);
            Cm_k_1m = getVar(Cc, m-k-1, Ndom, dom); dCm_k_1m = diff(Cm_k_1m);
            Cm_k_3p = getVar(Cc, m-k+3, Ndom, dom); dCm_k_3p = diff(Cm_k_3p);
            ddCm_k_1m = diff(Cm_k_1m,2); ddCm_k_1p = diff(Cm_k_1p,2);

            s1 = s1 - (1/(4*1i))*dAk.*( dCm_k_3 + dCm_k_1p - dCm_k_1m - dCm_k_3p );
            s2 = s2 - 0.125*( dBk - invr.*Bk ).*( dCm_k_3 - dCm_k_1p - dCm_k_1m + dCm_k_3p );
            s3 = s3 + (1/(2*1i))*Ak.*( ddCm_k_1m - ddCm_k_1p );
            s4 = s4 + 0.5*( invr.*Bk ).*( dCm_k_1m + dCm_k_1p );

            r1 = r1 + (1i/4)*dAk.*( dCm_k_3 - dCm_k_1p + dCm_k_1m - dCm_k_3p );
            r2 = r2 + 0.25*( dBk + invr.*Bk ).*( dCm_k_1m + dCm_k_1p );
            r3 = r3 + 0.125*( invr.*Bk - dBk ).*( dCm_k_3 + dCm_k_1p + dCm_k_1m + dCm_k_3p );
        end

        LHSb = LHSb + s1 + s2 + (1/(2*1i))*(lambda+1)*r.*( ddCm1 - ddCp1 ) + s3 + s4;
        RHSb = r1 + r2 + r3;

        eqCb{idxOf(m)} = LHSb - RHSb;
    end

    out = [eqB{:}, eqCa{:}, eqCb{:}].';
    %out = out(:);
end