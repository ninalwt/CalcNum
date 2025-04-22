function i = decomp_jacobi(M,Pb,x0) %decomposicao Jacobi
    n = size(M, 1);
    D = diag(diag(M));
    C = eye(n) - D \ M;
    g = D \ Pb;
    i = 0;

    % Loop principal
    while (norm(Pb - M*x0) > 1e-4 && i < 10000)
        i = i+1;
        x0 = C*x0+g;
    end

    % Se atingir o número máximo de iterações, exibe mensagem de erro
    if (i == 10000)
        disp('Erro: o método (Jacobi) não converge.');
        return;
    end
end