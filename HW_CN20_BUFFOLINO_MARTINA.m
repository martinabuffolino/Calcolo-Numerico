% MARTINA BUFFOLINO -- n.matr. 0124001850 
% le function implementate sono:
% -- sz1 = Secanti(f, 105, -0.7, 1e-10, 10);
% -- z1 = fzero(f, 105);
% -- [L,U] = lu(A);
% -- mx=media_pesata_001850(x)
% -- s=saxpy_001850(x,y,a)
% 
% Punto 1.
% Inizialmente dichiaro 3 variabili: a b (estremi dell'intervallo del 
% vettore T) e seed.
% Successivamente assegno un valore a piacere al seed per poi utilizzarlo
% nella funzione rng.
% Genero un vettore T di numeri casuali utilizzando randi(generatore di
% numeri interi).
% Randi in input prende l'intervallo [a b] e il numero di componenti da
% generare.
a = 100; 
b = 999;
seed = 438519;
rng(seed);
T = randi ([a b],1,30);

% Punto 2.
% Inizialemnte dichiaro la variabile divisore.
% Ordino in modo crescente mio vettore T ultilizzando sort ed elimino 
% eventuali ripetizioni utilizzando unique.
% Per eliminare tutte le componenti multiple di 6 utilizzo mod: se il
% resto della divisione tra l'elemento del vettore e il divisore è
% uguale 0 elimino l'elemento.
% Infine vado a definire il minino(m) e il massimo(M) del vettore X.
% Per controllare se il vettore X ha un numero di componeti >20 e <30
% uso l'operatore if; nel caso in cui la condizione non si verifica
% modifico il valore del seed.
divisore = 6;
X = sort(T);
X = unique(X);
X(~mod(X,divisore))=[]
M = max(X);
m = min(X);
if (length(X)>20) && (length(X)<30)
    length(X)
else 
    seed = 527194;
end

% Punto 3.
% Definisco la griglia fitta di X utilizzando il min e max su 20 punti.
% Definisco un handle per la funzione f composta da una funzione 
% polinomiale p e una funzione trigonomatrica g.
% Per verificare se f cambia segno lancio la Y con la funzione sign.
% Sign mi ritorna un array della stessa lunghezza di X; stampa :
% 1 se l'elemento è maggiore di 0;
% 0 se l'elemento è uguale a 0;
% -1 se l'elemento è minore di 0;
X=linspace(m,M,20);
p=@(X) 2*X.^2 + 5*X + 4;
g=@(X) sin(X);
f=@(x) g(p(X));
Y = sign(X);

% Punto 4.
% Metodo delle Secanti
% Alla function Secanti passo in input la funzione, l'intervallo di
% approsimazione, il delta e il numero massimo di interazioni.
% L'accuratezza è uguale alla differenza tra la soluzione esatta(z1) e la
% soluzione approssimata(sz1).
% Con acc<delta si effettua un controllo logico. Se è uguale a 1 vuol dire
% che l'accuratezza è più piccola della tolleranza richiesta.
delta = 4e-10;
sz1 = Secanti(f, 105, -0.7, 1e-10, 10);
z1 = fzero(f, 105);
acc = abs(z1-sz1);
acc<delta;

% Punto 5.
% Per disegnare graficamente la funzione basta plottare x,y 
X=linspace(m,M,20);
p=@(X) 2*X.^2 + 5*X + 4;
g=@(X) sin(X);
f=@(x) g(p(X));
Y=f(X);
hold on
grid on
plot(X,Y, 'b','LineWidth',1.5);

text(sz1,0, 'szero');

% Punto 6.
% Inizialmente dichiaro due variabili L e R che utilizzo come intervallo
% per generare la matrice. Una volta generata la matrice quadrata A
% verifico che non sia singolare. Per verificare che la matrice non sia
% singolare calcolo il determinante; se il determinate è uguale a 0
% reiposto il seed e genero una nuova matrice, in caso contrario vuol dire
% che la matrice non è singolare.
L = 2.8;
R = 11.7;
A = L+(R-L)*rand(5,5);
if det(A)==0
    seed = 842637;
    rng(seed);
    A= L+(R-L)*rand(5,5);
end

% Punto 7.
% Inizialmente dichiaro due variabili [c d] che uso come estremi 
% dell'intervallo del vettore colonna di 5 componenti.
% Per determinare il vettore colonna b moltiplico la matrice A per il 
% vettore colonna x_true.
c = 10;
d = 50;
x_true = randi ([c d],1,5)';
b = A*x_true;

% Punto 8.
% Metodo di Fattorizzazione LU
% Possiamo considerare A come il prodotto di due matrici (L*U) dello
% stesso ordine dove L è la matrice triangolare inferiore invece U è la 
% matrice triangolare superiore.
% L'algoritmo si divide in 3 passi: 
% 1.bisogna fattorizzare A=L*U --> [L,U] = lu(A);
% 2.calcoliamo la p(con STriangInf) Lp=b sapendo che LUx=b U*x=p --> y=L\b;
% 3.calcoliamo la x(con STriangSup) Ux=p --> x_finale = U\y;
% Per calcolare l'errore confronto la x (calcolata col sistema lineare) 
% con la x_finale (calcolata con STriangSup) sono uguali;
b = A*x_true;
[L,U] = lu(A);
x = A\b;
y = L\b;
x_sol = U\y;

errore = abs(x_sol-x_true);

% Punto 9.
% Per invertire due righe della matrice A basta cambiare gli indici di
% riga. Una volta calcolata la media tra L e R scorro con due cicli
% annidati sulle righe e sulle colonne; se i valori dell'i-esima riga e
% dell'i-esima colonna sono minori della media devo sostutuili con essa.
A([1 2],:) = A([2 1],:);
intervallo = [L,R];
media = mean(intervallo);
for i = 1:5
    for j = 1:5
        if A(i,j) < media
            A(i,j)=media;
        end
    end
end

% Punto 10.
% Per creare il vettore x estraggo 10 componenti non consecutive di X al 
% al passo 2. Dopodichè creo la griglia fissa passando in input gli estremi
% e N punti.
x = X(1:2:20);
min = min(x);
max = max(x);
xx = linspace(min,max,150);

% Punto 11.
% Creo la griglia di valutazione passando gli estremi e i numero di punti; 
% in questo caso utilizziamo 10 punti perchè il grado del polinomio è 9 
% quindi n-1. 
% Con polyfit calcoliamo i coefficienti passando ascisse,ordinate e
% lunghezza di x. 
% Con polyval valutiamo il polinomio interpolante sui coefficienti
x=linspace(min,max,10);
f=@(x) g(p(x));
y=f(x);
xx = linspace(min,max,100);
coeff_int= polyfit(x,y,length(x)-1);
yy_int = polyval(coeff_int,xx);

% Creo la griglia di valutazione passando gli estremi e i numero di punti; 
% in questo caso utilizziamo 4 punti perchè il grado del polinomio è 3 
% quindi n-1. 
% Con polyfit calcoliamo i coefficienti passando ascisse,ordinate e
% lunghezza di x. 
% Con polyval valutiamo il polinomio interpolante sui coefficienti
x1 = linspace(min,max,4);
f = @(x1) g(p(x));
y = f(x);
coeff_app = polyfit(x1,y,length(x1)-1);
yy_app = polyval(coeff_app,xx);

% Creo la griglia di valutazione passando gli estremi e i numero di punti; 
% in questo caso utilizziamo 3 punti perchè il grado del polinomio è 2 
% quindi n-1. 
% Con polyfit calcoliamo i coefficienti passando ascisse,ordinate e
% lunghezza di x. 
% Con polyval valutiamo il polinomio interpolante sui coefficienti
x2 = linspace(min,max,3);
f = @(x2) g(p(x));
y = f(x);
coeff = polyfit(x2,y,length(x2)-1);
yy_trig = polyval(coeff,xx);

% Punto 12.
figure (1)
hold on
grid on
%f in x
plot(x,y,'o-','MarkerSize',6,'LineWidth',1.5) 

%pol_int in xx
plot(xx,yy_int,'red','LineWidth',1.5) 

%pol_app in xx
plot(xx,yy_app,'g-','LineWidth',1.5) 

%pol_int in xx
plot(xx,yy_trig,'m-','LineWidth',1.5) 

% Punto 13.
mx=media_pesata_001850(x)

% Punto 14.
s=saxpy_001850(x,y,a)

% Punto 15.
% La function prende in input la il vettore x di cui calcolo la lunghezza.
% Attraverso il for scorro su tutto il vettore impostando all'interno tre
% condizioni.
function mx=media_pesata_001850(x)
    n = length(x);
    for i = 1: n
        if(i==1) 
            mx(i) = (x(1)+x(2))/2;
        end
        if(i>=2 && i <= n-1)
            mx(i) = (x(i-1)+2*x(i)+x(i+1))/4;
        end
        if(i==n)
            mx(i) = (x(n-1)+x(n))/2;
        end         
    end
end

% Punto 16.
% La function prende in input il vettore x di cui calcolo la lunghezza,
% il vettore y e un numero reale a.
function s=saxpy_001850(x,y,a)
    n = length(x);
    for i= 1: n
        s(i) = (a*x(i)+ y(i));
    end
end

