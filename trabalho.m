arestas = load('manh.el') + 1; %Le o arquivo e soma um aos indices do grafo, pois eles precisam ser maiores que 0

G = graph(arestas(:, 1), arestas(:, 2)); %Grafo das arestas (arestas(:, 1) = primeira coluna)

selecionado = 15;
tamanho = numnodes(G);
vertices_aleatorios = randperm(tamanho, selecionado); %Escolhe um subgrafo aleatorio de 15 vertices
valores_aleatorios = rand(1, selecionado) * selecionado; %Atribui valores aleatorios aos vertices

coordenadas = load('manh.xy');
figure(1)
scatter(coordenadas(:,1),coordenadas(:,2),6,'.');
title('Ilha de Manhattan')

A = adjacency(G);
D = diag(degree(G));
L = D - A;

P = sparse(vertices_aleatorios, vertices_aleatorios, 1e7, tamanho, tamanho);

b = zeros(tamanho, 1);
b(vertices_aleatorios) = valores_aleatorios;

M = L + P;
Pb = P * b;

tic, x1 = decomp_LU(M,Pb); tempo_lu = toc;
tic, x2 = decomp_cholesky(M,Pb); tempo_ch = toc;

x0 = zeros(n,1);

tic, i1 = decomp_jacobi(M,Pb,x0); tempo_ja = toc;
tic, i2 = decomp_gauss_seidel(M,Pb,x0); tempo_gs = toc;
tic, i3 = decomp_grad(M,Pb,x0); tempo_gr = toc;

disp(['Tempo:\n ' 
+ 'LU -> ' + tempo_lu + '.\n'
+ 'Cholesky -> ' + tempo_ch + '.\n '
+ 'Jacobi -> ' + tempo_ja + '.\n'
+ 'Gaussseidel -> ' + tempo_gs + '.\n '
+ 'Gradientes -> ' + tempo_gr + '.\n '
]);
