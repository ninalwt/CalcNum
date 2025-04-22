function i = decomp_gauss_seidel(M,Pb,x0)
    INF = tril(M);
    SUP = triu(M,1);
    C = -INF \ SUP; % Calculando a matriz de iteração
    g = INF \ Pb; % Calculando o vetor de iteração
    i = 0;

    while (norm(Pb-M*x0) > 1e-4 && i < 10000)
        i = i+1;
        x0 = C*x0+g;
    end

    % Se atingir o número máximo de iterações, exibe mensagem de erro
    if (i == 10000)
        disp('Erro: o metodo (Seidel) nao converge.');
        return;
    end
end