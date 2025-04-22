function x = decomp_cholesky(M,Pb) %decomposicao de cholesky
    H = chol(M);
    
    y = H'\Pb;
    x = H\y;

    return
end