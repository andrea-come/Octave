clc
clear
% Esercizio 1.3: si dimostri che i^i è un numero reale
printf("Esercizio 1.3\n")
printf("i^i = %f, Im(i^i) = %f\n", i^i, imag(i^i))
display("")



% Esercizio 1.4: si costruiscano una matrice triangolare superiore
%  ed una triangolare inferiore di dimensione 10
%  con 2 sulla diagonale principale e -3 sulla seconda sopra
%  (rispettivamente, sotto) diagonale.
printf("Esercizio 1.4\n")
n = 10;
L = 2*eye(n) -3 * diag(ones(n-1,1),-1)
U = 2*eye(n) -3 * diag(ones(n-1,1),+1)
display("")



% Esercizio 1.5: si scrivano le istruzioni che consentono di scambiare fra
%  loro la terza e la settima riga delle matrici costruite nell’Esercizio 1.4,
%  indi quelle per scambiare l’ottava con la quarta colonna
printf("Esercizio 1.5\n")
Lr = L;
temp = Lr(7,:);
Lr(7,:) = Lr(3,:);
Lr(3,:) = temp
% Alternativa, lo stesso vettore può essere usato per più matrici
r = [1:n];
r(3) = 7;
r(7) = 3;
Ur = U(r,:)

c = [1:n];
c(8) = 4;
c(4) = 8;
Lc = L(:,c)
Uc = U(:,c)
display("")



% Esercizio 1.6: si stabilisca se i seguenti vettori di R4 sono fra loro
%  linearmente indipendenti
printf("Esercizio 1.6\n")
v1 = [0 1 0 1];
v2 = [1 2 3 4];
v3 = [1 0 1 0];
v4 = [0 0 1 1];

printf("I vettori %s, %s, %s e %s sono ", mat2str(v1),...
        mat2str(v2), mat2str(v3), mat2str(v4))
if det([v1; v2; v3; v4]) == 0
  printf("linearmente dipendenti.\n")
else
  printf("linearmente indipendenti.\n")
endif
display("")



% Esercizio 1.7: si scrivano le seguenti funzioni e si calcolino con
%  il toolbox simbolico derivata prima e seconda ed integrale indefinito
printf("Esercizio 1.7\n")
pkg load symbolic

syms x;
f = sqrt(x^2 + 1);
g = sin(x^3) + cosh(x);

display(f)
printf("f' = ")
display(simplify(diff(f)))
printf("f'' = ")
display(simplify(diff(f, 2)))
printf("∫f = ")
display(simplify(int(f)))

display(g)
printf("g' = ")
display(diff(g))
printf("g'' = ")
display(simplify(diff(g, 2)))
printf("∫g = ")
display(int(g))
printf("F indica la Funzione Ipergeometrica Generalizzata\n")
display("")
% La F che appare nell'ultimo integrale è la
%  Funzione Ipergeometrica Generalizzata con un parametro superiore (2/3)
%  e due parametri inferiori (3/2, 5/3).
% Dare un'occhiata a "Meijer G-function".



% Esercizio 1.8: dato un vettore v di dimensione n, scrivendo c=poly(v) è pos-
%  sibile costruire gli n+1 coefficienti del polinomio p(x) che ha come
%  radici le componenti del vettore v.
% In aritmetica esatta si ha v = roots(poly(v)), tuttavia ciò potrebbe non
%  verificarsi a causa degli errori di arrotondamento, come si può constatare
%  richiamando il comando roots(poly([1:n])), dove n varia da 2 fino a 25.
printf("Esercizio 1.8\n")
format long
for n = 12:3:21 % gli errori di arrotondamento iniziano per n = 13
  printf("n = %d\n", n)
  display(roots(poly([1:n])))
  display("")
endfor
format short
% Da n = 21 iniziano ad apparire radici complesse



% Esercizio 1.9: si scriva un programma per il calcolo della uccessione
%
%  I_0 = (e − 1)/e
%  I_(n+1) = 1 − (n + 1)I_n , per n = 0, 1, . . . , 21.
%
% Sapendo che In → 0 per n → ∞, si commentino i risultati ottenuti.
printf("Esercizio 1.9\n")
n_max = 21;
I = (e - 1)/e;
printf("n =  0, I = %.4e\n", I)
for n = 1:(n_max + 1)
  I = 1 - n*I;
  printf("n = %2d, I = %.4e\n", n, I)
endfor
printf("I = %.4e\n\n", I)

% Alternativa che salva i valori della successione
I = zeros(n_max + 2, 1);
I(1) = (exp(1) - 1)/exp(1);
for i = 0:n_max
  I(i+2) = 1 - (i+1)*I(i+1);
endfor



% Esercizio 1.10: si spieghi il comportamento della successione (1.4)
printf("Esercizio 1.10\n")
n_max = 30;
zs = zeros(n_max, 1);
errors = zeros(n_max, 1);
zs(1) = 2;
errors(1) = abs(2 - pi);

for n = 2:n_max
  zs(n) = 2^(n-0.5) * sqrt(1 - sqrt(1 - 4^(1-n) * (zs(n-1))^2 ) );
  errors(n) = abs(zs(n) - pi);
endfor

fprintf('n = %3d | z = %.6f | err = %e\n', [(1:n_max)', zs, errors]');

% Creazione della figura
figure(1); clf;
semilogy(1:n_max, errors, '-o', 'LineWidth', 2, 'MarkerSize', 4);

% Limiti e griglia asse y
ylim([1e-15 1]);
esponenti = -15:2:0;
tacche = 10.^esponenti;
set(gca, 'YTick', tacche);
grid on;
set(gca, 'YMinorGrid', 'off');

% Etichette e Titolo
xlabel('n', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('|z_n - \pi|', 'FontSize', 12, 'FontWeight', 'bold');
title('Instabilità numerica nel calcolo di \pi', 'FontSize', 14);

% Leggibilità degli assi
set(gca, 'FontSize', 12); % Aumenta la dimensione del testo degli assi
xlim([1 n_max]); % Fissa i limiti dell'asse X
%legend('Errore Assoluto', 'Location', 'northeast');

% Il problema della cancellazione numerica è presente nell'istruzione:
%
%  1 - sqrt(1 - 4^(1-n))
%
% essendo 4^(1-n) piccolo rispetto a 1.
% Conviene portare la radice al denominatore
printf("\nRiduzione dell'errore ottenuta con formula eqivalente\n")
for n = 2:n_max
  term = 4^(1-n) * (zs(n-1))^2;
  zs(n) = 2^(n-0.5) * sqrt( term / (1 + sqrt(1 - term)) );
  errors(n) = abs(zs(n) - pi);
endfor

fprintf('n = %3d | z = %.6f | err = %e\n', [(1:n_max)', zs, errors]');

hold on
semilogy(1:n_max, errors, '-o', 'LineWidth', 2, 'MarkerSize', 4);
display("")



% Esercizio 1.11: per il calcolo di π si può usare la seguente tecnica:
%  si generano n coppie {(x_k , y_k )} di numeri casuali compresi fra 0 e 1
%  e di questi si calcola il numero m di punti che cadono nel primo quarto
%  del cerchio di centro l’origine e raggio 1.
% Si ha che π è il limite per n che tende all’infinito dei rapporti
%  π_n = 4m/n. Si scriva un programma che esegua questo calcolo e si
%  verifichi la correttezza del risultato al crescere di n.
printf("Esercizio 1.11\n")
printf("*********** Singola esecuzione ***********\n")
n = 6000000;

tic;
xs = rand(n, 1); ys = rand(n, 1);
pi_approx = 4 * sum(xs.^2 + ys.^2 <= 1) / n;
tempo = toc;

fprintf('pi:       %.8f\n', pi_approx);
fprintf('Errore:   %.4e\n', abs(pi_approx - pi));
fprintf('Tempo:    %.2f secondi\n', tempo);
printf("__________________________________________\n\n")

% Usando il metodo a blocchi
printf("*********** Metodo a blocchi *************\n")
n_total = 100000000;
block_size = 5000000; % 5 milioni alla volta (leggero per la RAM)
hits = 0; % Contatore successi accumulati

% Calcoliamo quanti blocchi servono
num_blocks = ceil(n_total / block_size);

fprintf('Calcolo in corso su %d blocchi...\n', num_blocks);

tic; % Avvia cronometro
for i = 1:num_blocks
    % Se è l'ultimo blocco, potremmo dover generare meno punti
    current_n = min(block_size, n_total - (i-1)*block_size);

    xs = rand(current_n, 1);
    ys = rand(current_n, 1);
    hits = hits + sum(xs.^2 + ys.^2 <= 1);
end
pi_approx = 4 * hits / n_total;
tempo = toc;

fprintf('pi:       %.8f\n', pi_approx);
fprintf('Errore:   %.4e\n', abs(pi_approx - pi));
fprintf('Tempo:    %.2f secondi\n', tempo);
printf("__________________________________________\n\n")

% Facendo la media
printf("*********** Calcolando la media **********\n")
n_per_run = 2000000;
num_runs = 50;
pi_estimates = zeros(num_runs, 1); % Vettore per salvare i risultati

fprintf('Avvio di %d simulazioni Montecarlo...\n', num_runs);

tic;
for i = 1:num_runs
    xs = rand(n_per_run, 1);
    ys = rand(n_per_run, 1);
    % Stima locale
    pi_estimates(i) = 4 * sum(xs.^2 + ys.^2 <= 1) / n_per_run;
end
% Calcolo delle statistiche
pi_mean = mean(pi_estimates);  % Media
pi_std  = std(pi_estimates);   % Deviazione standard
tempo = toc;

fprintf('Valore medio di pi:    %.8f\n', pi_mean);
fprintf('Errore:                %.4e\n', abs(pi_mean - pi));
fprintf('Deviazione Standard:   %.6f\n', pi_std);
fprintf('Tempo:                 %.2f secondi\n\n', tempo);

% Grafico della distribuzione
figure(2); clf;
hist(pi_estimates, 10);
title('Distribuzione delle stime di \pi');
xlabel('Valore stimato');
ylabel('Frequenza');
grid on;
hold on;

yl = ylim;

h1 = plot([pi_mean pi_mean], yl, '-r', 'LineWidth', 2);   % rossa continua
h2 = plot([pi pi],           yl, '--g', 'LineWidth', 2);  % verde tratteggiata

legend('Distribuzione', 'Media Calcolata', 'Vero \pi');
hold off;



% Esercizio 1.12 Sempre per il calcolo di π si può utilizzare una troncata
%  della serie con termine generale:
%
%  16^(-n) * ( 4/(8n + 1) - 2/(8n + 4) - 1/(8n + 5) - 1/(8n + 6) )
%
% Si realizzi una function che ne calcoli la somma fino ad un certo n fissato.
% Quanto grande deve essere n per ottenere un valore di π confrontabile con
%  quello memorizzato nella variabile pi?
% La formula usata è detta "formula di Bailey-Borwein-Plouffe".
printf("Esercizio 1.12\n")
function pi_approx = es_1_12(n)
  i = 0:n;
  terms = 16.^(-i) .* (4./(8*i + 1) - 2./(8*i + 4) - 1./(8*i + 5) - 1./(8*i + 6));
  pi_approx = sum(terms);
end
n = 6;
pi_approx = es_1_12(n)



% Esercizio 1.13: si scriva un programma per il calcolo del coefficiente
%  binomiale n!/(k!(n − k)!), dove n e k sono numeri naturali con k ≤ n.
printf("Esercizio 1.13\n")
n = 50; k = 6;
% simmetria del binomiale per ridurre i calcoli
if k > n - k
  k = n - k;
endif
% Calcolo iterativo per evitare overflow intermedi
bin = 1;
for i = 1:k
    bin = bin * (n - i + 1) / i;
endfor
printf("Il binomiale 50 su 6 vale: %d\n\n", bin)



% Esercizio 1.14: si realizzi una function che calcoli l’elemento f_n della
%  successione di Fibonacci in forma ricorsiva. Osservando poi che
%
% [f_i; f_{i-1}] = [1 1; 1 0] [f_{i-1}; f_{i-2}]
%
%  si realizzi un’altra function che calcoli f_n sfruttando questa relazione.
% Si confrontino i relativi tempi di calcolo.
printf("Esercizio 1.14\n")
function f_n = fib_ric(n)
  persistent valori; % Il vettore per memorizzare valori già calcolati

  if isempty(valori)
    valori = [0, 1]; % casi base
  endif

  % Se n+1 è dentro i limiti del vettore allora è stato già calcolato
  if (n + 1) <= length(valori)
    f_n = valori(n + 1);
    return;
  endif

  val = fib_ric(n - 1) + fib_ric(n - 2);
  valori(n + 1) = val;
  f_n = val;
endfunction

function f_n = fib_matrice(n)
  A = [1 1; 1 0];
  f = A^n * [0; 1];
  f_n = f(1);
endfunction

n = 40;

t_0 = cputime;
f_n = fib_ric(n);
t = cputime - t_0;
printf("Usando la ricorsione, con dizionario:\n")
printf("f_%d = %d, tempo impiegato: %f\n\n", n, f_n, t)

t_0 = cputime;
f_n = fib_matrice(n);
t = cputime - t_0;
printf("Usando la matrice:\n")
printf("f_%d = %d, tempo impiegato: %f\n\n", n, f_n, t)



% Esercizio 1.15: sappiamo che lim_{n→∞} (1 + 1/n )^n = e.
% È sensato approssimare numericamente il numero di Nepero e con un
%  valore molto alto di n?
printf("Esercizio 1.15\n")
exp_max = 18;
n = 10.^[1:exp_max];
e_n = (1 + 1./n).^n;
err_n = abs(e_n - e);

for i = 1:length(e_n);
  printf("n = %.2e, e_n = %f, errore = %e\n", 10^i, e_n(i), err_n(i))
endfor

figure(3)
semilogx(n, e_n, 'bo');
hold on; grid on;
semilogx([1, 10^exp_max], [e, e], 'r--')

xlim([1 10^exp_max]);
ylim([1 3.2])
esponenti = 0:3:exp_max;
tacche = 10.^esponenti;
set(gca, 'XTick', tacche);
set(gca, 'YMinorGrid', 'off');
set(gca, 'XMinorGrid', 'off');
