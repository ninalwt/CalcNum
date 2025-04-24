% LEITURA DOS DADOS
disp('oi');
arestas = load('manh.el') + 1;
coordenadas = load('manh.xy');

G = graph(arestas(:,1), arestas(:,2)); % Grafo completo
A = adjacency(G);

% EXTRAINDO A MAIOR COMPONENTE CONEXA
bins = conncomp(G); % Identifica componentes conexas
contagens = histcounts(bins, 1:max(bins)+1);
[~, maior_bin] = max(contagens); % Pega o maior grupo
idx_maior = find(bins == maior_bin); % Índices dos vértices dessa componente

G_sub = subgraph(G, idx_maior); % Subgrafo com a maior componente
A_sub = adjacency(G_sub);
n = numnodes(G_sub);

% PLOT DA COMPONENTE CONEXA
figure(1)
coords_sub = coordenadas(idx_maior,:);
scatter(coords_sub(:,1), coords_sub(:,2), 6, '.');
    
% MATRIZ LAPLACIANA
D = diag(degree(G_sub));
L = D - A_sub;

% FONTES DE CALOR ALEATÓRIAS EM 1% DOS NÓS
num_fontes = ceil(0.01 * n);
idx_fontes = randperm(n, num_fontes);
valores_fontes = rand(num_fontes,1) * 10;

P = sparse(idx_fontes, idx_fontes, 1e7, n, n); % Penalização grande
b = zeros(n, 1);
b(idx_fontes) = valores_fontes;

M = L + P;
Pb = P * b;

% RESOLUÇÃO DO SISTEMA LINEAR
x0 = zeros(n,1);

% LU
tic, x_lu = decomposicao_LU(M,Pb); tempo_lu = toc;

% Cholesky
tic, x_ch = decomposicao_cholesky(M,Pb); tempo_ch = toc;

% Jacobi
tic, [x_ja, i_ja] = decomposicao_jacobi(M,Pb,x0); tempo_ja = toc;

% Gauss-Seidel
tic, [x_gs, i_gs] = decomposicao_gauss_seidel(M,Pb,x0); tempo_gs = toc;

% Gradientes Conjugados
tic, [x_gr, i_gr] = decomposicao_gradientes(M,Pb,x0); tempo_gr = toc;

% PLOTS DAS SOLUÇÕES
figure(2)
x = coords_sub(:,1);
y = coords_sub(:,2);
scatter3(x, y, x_lu, 2, x_lu)
colorbar
colormap jet
title('Distribuição de Temperatura (Solução via LU)')

figure(3)
x = coords_sub(:,1);
y = coords_sub(:,2);
scatter3(x, y, x_ch, 2, x_ch)
colorbar
colormap jet
title('Distribuição de Temperatura (Solução via Cholesky)')

figure(4)
x = coords_sub(:,1);
y = coords_sub(:,2);
scatter3(x, y, x_ja, 2, x_ja)
colorbar
colormap jet
title('Distribuição de Temperatura (Solução via Jacobi)')

figure(5)
x = coords_sub(:,1);
y = coords_sub(:,2);
scatter3(x, y, x_gs, 2, x_gs)
colorbar
colormap jet
title('Distribuição de Temperatura (Solução objetivo principal é observar o tempo de solução do sistema construído para diferentesvia Gauss-Seidel)')

figure(6)
x = coords_sub(:,1);
y = coords_sub(:,2);
scatter3(x, y, x_gr, 2, x_gr)
colorbar
colormap jet
title('Distribuição de Temperatura (Solução via Gradientes Conjugados)')



% TEMPOS DE EXECUÇÃO
fprintf('\nTempos de Execução:\n');
fprintf('LU: %.4f s\n', tempo_lu);
fprintf('Cholesky: %.4f s\n', tempo_ch);
fprintf('Jacobi: %.4f s (%d iterações)\n', tempo_ja, i_ja);
fprintf('Gauss-Seidel: %.4f s (%d iterações)\n', tempo_gs, i_gs);
fprintf('Gradientes: %.4f s (%d iterações)\n', tempo_gr, i_gr);

