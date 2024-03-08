//Definition du systeme d'equations differentielles, au debut fichier
function [v] = auto_repressor_rhs(t,x,pr)
    //t: instant de temps
    //x: variables d'etat
    //pr: parametres du systeme; reprendre depuis la definition
    //Renommer les parametres et les variables
    alpha = pr(1);
    beta = pr(2);
    kP = pr(3);
    gM = pr(4);
    gP = pr(5);
    n = pr(6);
    M=x(1);
    P=x(2);
    //Initialisation du vecteur des derivees
    v = [];
    v(1) = alpha*(kP^n) / (kP^n + P^n) - gM * M;
    v(2) = beta*M - gP * P;
endfunction;
//Partie principale
//Definition de parametres; vecteur de parametres;
//BIEN RESPECTER L’ORDRE DE DEFINITION DE pr
alpha = 3.1; beta = 2.3;
kP = 10;
gM = 0.1; gP = 0.01;
n = 2;
par = [alpha,beta,kP,gM,gP,n];
//Definition de condition(s) initiale(s),
//instant initial et vecteur de temps de simulation
x0 = [10;1];
t0 = 0;
dt = 0.05;
tvec = t0:dt:30;
//Commande de resolution du systeme d'equations differentielles
//methodes numeriques disponibles: “rk”, “stiff”,...
sol=ode("stiff",x0,t0,tvec, list(auto_repressor_rhs,par) );
//Faire le graphe de la solution au cour du temps, dans la figure 1
//Options: couleur (r=red, b=blue,k=black,...), forme de ligne (-,--,:,-.)
figure(1);
plot(tvec,sol(1,:),'b-',tvec,sol(2,:),'r-');