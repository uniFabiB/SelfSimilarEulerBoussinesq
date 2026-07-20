N = chebop(@(x,u,v) [diff(u) - u.*(1 - v);
                      diff(v) - v.*(u - 2)], [0 5]);

N.lbc = @(u,v) [u(0)-0.5; v(0)-0.3];

[u,v] = N\[0;0];