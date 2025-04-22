
function i = decomp_grad(M,Pb,x0) %decomposicao de gradientes conjugados
    r0 = M*x0-Pb;
    p0 = -r0;
    alfa = r0'*r0/(r0'*M*r0);
    x = x0+alfa*p0;
    r = r0+alfa*M*p0;
    p = p0;
    i = 1;

    while norm(x-x0,2) >= norm(x,2)*1e-14 && i < 10000
        x0 = x;
        p0 = p;
        a  = r'*r/(r0'*r0);
        p = -r+a*p0;
        alfa = r'*r/(p'*M*p);
        r0 = r;
        x = x0+alfa*p;
        r = r0+alfa*M*p;
        i = i+1;
    end

    return
end