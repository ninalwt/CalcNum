function x = decomp_LU(M,Pb) %decomposicao LU
    [INF,SUP] = lu(M);

    y = INF\Pb;
    x = SUP\y;

    return
end