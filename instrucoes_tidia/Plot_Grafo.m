%% Exemplo de como plotar um grafo a partir da matriz de adjacencia.
%% Para isso, modifiquei o código 'Teste_Split_Edges.m' (em anexo) que o professor fez.
%% O professor fez o código com base na figura 'Grafo.jpeg' (em anexo).
%% Obs: A função usada para plotar o grafo funciona apenas para o Matlab,
%% ou seja, não funciona para o Octave.

%% E-mail para tirar dúvidas: aquisson@usp.br.

nv = 7; %% numero de vertices.
%% Construindo a matriz de adjacência/incidência:
A  = zeros(nv);
A(1,2) = 1; A(2,1) = 1;
A(1,3) = 1; A(3,1) = 1;
A(1,4) = 1; A(4,1) = 1;
A(2,3) = 1; A(3,2) = 1;
A(3,4) = 1; A(4,3) = 1;
A(3,5) = 1; A(5,3) = 1;
A(4,5) = 1; A(5,4) = 1;
A(6,7) = 1; A(7,6) = 1;

%% Plotando o grafo inicial com base na matriz acima.
%% Note que o grafo inicial possui duas comp. conexas.
figure(1)
G = graph(A,'upper');
h1 = plot(G);

%% Calculando todas as componentes conexas (duas, neste caso).
[nc nvc vc] = Split_Edges(nv,A);
vc = vc(1,:); % selecionando a maior componente conexa.
%% Construindo a matriz de adjacência da maior componente conexa.
A = A(vc,vc);
%% Escolhendo alguns vértices no grafo da maior componente conexa.
p = [1 5]; 
%% Plotando o grafo da maior componente conexa e destacando os vértices selecionados.
figure(2)
G = graph(A,'upper');
h2 = plot(G);
highlight(h2,p,'NodeColor','red')
