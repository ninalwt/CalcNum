arestas = load('manh.el') + 1; %Le o arquivo e soma um aos indices do grafo, pois eles precisam ser maiores que 0
G = graph(arestas(:, 1), arestas(:, 2)); %Grafo das arestas (arestas(:, 1) = primeira coluna)
selecionado = 15;tamanho = numnodes(G);vertices_aleatorios = randperm(tamanho, selecionado); %Escolhe um subgrafo aleatorio de 15 verticesvalores_aleatorios = rand(1, selecionado) * selecionado; %Atribui valores aleatorios aos vertices
coordenadas = load('manh.xy');figure(1)scatter(coordenadas(:,1),coordenadas(:,2),6,'.');title('Ilha de Manhattan')
A = adjacency(G);D = diag(degree(G));L = D - A;
P = sparse(vertices_aleatorios, vertices_aleatorios, 1e7, tamanho, tamanho);
b = zeros(tamanho, 1);b(vertices_aleatorios) = valores_aleatorios;
M = L + P;Pb = P * b;
tic, x1 = decomposicao_LU(M,Pb); tempo_lu = toc;tic, x2 = decomposicao_cholesky(M,Pb); tempo_ch = toc;
x0 = zeros(n,1);
tic, i1 = decomposicao_jacobi(M,Pb,x0); tempo_ja = toc;tic, i2 = decomposicao_gaussseidel(M,Pb,x0); tempo_gs = toc;tic, i3 = decomposicao_gradientes(M,Pb,x0); tempo_gr = toc;
disp(['Tempo:\n ' + 'LU -> ' + tempo_lu + '.\n'+ 'Cholesky -> ' + tempo_ch + '.\n '+ 'Jacobi -> ' + tempo_ja + '.\n'+ 'Gaussseidel -> ' + tempo_gs + '.\n '+ 'Gradientes -> ' + tempo_gr + '.\n ']);
function x = decomposicao_LU(M,Pb)    [INF,SUP] = lu(M);
    y = INF\Pb;    x = SUP\y;
    returnend
function x = decomposicao_cholesky(M,Pb)    H = chol(M);        y = H'\Pb;    x = H\y;
    returnend
function i = decomposicao_jacobi(M,Pb,x0)    n = size(M, 1);    D = diag(diag(M));    C = eye(n) - D \ M;    g = D \ Pb;    i = 0;
    % Loop principal    while (norm(Pb - M*x0) > 1e-4 && i < 10000)        i = i+1;        x0 = C*x0+g;    end
    % Se atingir o número máximo de iterações, exibe mensagem de erro    if (i == 10000)        disp('Erro: o método (Jacobi) não converge.');        return;    endend
function i = decomposicao_gaussseidel(M,Pb,x0)    INF = tril(M);    SUP = triu(M,1);    C = -INF \ SUP; % Calculando a matriz de iteração    g = INF \ Pb; % Calculando o vetor de iteração    i = 0;
    while (norm(Pb-M*x0) > 1e-4 && i < 10000)        i = i+1;        x0 = C*x0+g;    end
    % Se atingir o número máximo de iterações, exibe mensagem de erro    if (i == 10000)        disp('Erro: o metodo (Seidel) nao converge.');        return;    endend
function i = decomposicao_gradientes(M,Pb,x0)    r0 = M*x0-Pb;    p0 = -r0;    alfa = r0'*r0/(r0'*M*r0);    x = x0+alfa*p0;    r = r0+alfa*M*p0;    p = p0;    i = 1;
    while norm(x-x0,2) >= norm(x,2)*1e-14 && i < 10000        x0 = x;        p0 = p;        a  = r'*r/(r0'*r0);        p = -r+a*p0;        alfa = r'*r/(p'*M*p);        r0 = r;        x = x0+alfa*p;        r = r0+alfa*M*p;        i = i+1;    end
    returnend
