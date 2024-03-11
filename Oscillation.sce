function hp = hillp(x,k,nn)
    hp = x^nn /(x^nn + k^nn);
endfunction

function [v] = equations(t,x,pr)
    //t: instant de temps
    //x: variables d'etat
    //pr: parametres du systeme; reprendre depuis la definition
    //Renommer les parametres et les variables

    W = x(1);
    Frq = x(2); // ARNm
    Pcy = x(3);
    Pnu = x(4);

    a1 = par(1);
    a2 = par(2);
    a3 = par(3);
    a4 = par(4);
    gamma = par(5);
    gM = par(6); // Correspond à gamma M
    gP = par(7); // doit être inférieur à  1/2 * a4/K4

    K1 = par(8);
    K4 = par(9);

    n = par(10);
    deg = gamma * Pnu * W
    //Initialisation du vecteur des derivees
        v = [];
    v(1) = a4 * hillp(W, K4, n) - deg - gP * W;
    v(2) = a1* hillp(W, K1, n) -gM * Frq ;
    v(3)= a2 * Frq -a3 * Pcy - gP * Pcy;
    v(4) = a3 * Pcy -gP * Pnu -deg ;
endfunction

a1 = 8;
a2 = 7;
a3 = 0.7;
a4 = 40;
gamma = 0.15;
gP = 0.01; // Correspond à gamma P
gM = 0.4; // doit être inférieur à  1/2 * a4/K4
K1 = 50;
K4 = 0.7;
n = 2;

par = [a1, a2, a3, a4, gamma, gM, gP, K1, K4, n];

x0 = [0.1; 1; 1; 5]; // Conditions initiales
t0 = 0;
dt = 0.01;
tvec = t0:dt:150;


sol=ode("stiff",x0,t0,tvec, list(equations,par) );
//Faire le graphe de la solution au cour du temps, dans la figure 1
//Options: couleur (r=red, b=blue,k=black,...), forme de ligne (-,--,:,-.)
figure(1);
plot(tvec,sol(1,:),'b-',tvec,sol(2,:),'r-',tvec,sol(3,:),'g-',tvec,sol(4,:),'k-');
legend('W', 'Frq', 'Pcy', 'Pnu');
