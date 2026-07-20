function out = genOp(r,A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,B6,C0,C1,C2,C3,C4,C5,C6)
    out = buildSystem(r, {A0,A1,A2,A3,A4,A5,A6,B0,B1,B2,B3,B4,B5,B6,C0,C1,C2,C3,C4,C5,C6}, 3, 0.5, [1, 5]);
end
