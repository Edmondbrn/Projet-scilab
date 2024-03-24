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

function Periode= calculperiode(X,tvec)
    //Initialisation
    Periode = 0 ;
    //X est un vecteur X=[X1,...,Xn]
    L=length(X);
    S=floor(0.3*L);
    MaxX=max(X(S:L));
    MinX=min(X(S:L));
    theta=0.5*(MaxX+MinX);
    j = 0;
    for i=S : L-1 
        auxg=X(i)-theta;
        auxd=X(i+1)-theta; 
        if  (auxg<0) &  (auxd>0)
            j=j+1;
            Pvec(j)=tvec(i);
        end 
    end
    
    //Continuation 
    //Pvec contient les instants de tous les croisements(t1,t2,t3,etc)
    Periode=mean(diff(Pvec));
    Pstd = stdev(diff(Pvec)); // vérifier ecart type entre les temps
    disp("L ecart type entre les différents temps est : ", Pstd); 
    
    //diff(Pvec)= [Pvec(2)-Pvec(1),Pvec(3)-Pvec(2),....,Pvec(n)-Pvec(n-1)] Explication de diff, pas besoin de le mettre dans le code. 
endfunction

function S = sensibilite(par, tvec)
    //Calcule de la période d'oscillation d'origine
    t0 = 0;
    X0 = [1;5;5;1];
    xsol=ode("stiff",X0,t0,tvec, list(equations,par) );
    periode = calculperiode(xsol(1,:), tvec);
    
    // Calcul de la sensibilité
    for j = 1:10 // on parcourt tous les paramètres
        parS = par; // on copie les paramètres pour les comparer plus tard
        parS(j) = par(j) * 4; // on augmente le paramètre j de 5%
        xsolS = ode("stiff", X0, t0, tvec, list(equations, parS)); // On résout le système avec les nouveaux paramètres
        Periode_S = calculperiode(xsolS(1,:), tvec); // ON calcule la période d'oscillation du nouveau système
        S(j) = (Periode_S - periode) / (parS(j) - par(j)); // On calcule la sensibilité
        disp(Periode_S, periode);
    end
    moyenne_S = mean(S.*S); // On calcule la moyenne des sensibilités
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

x0 = [1; 5; 5; 1]; // Conditions initiales
t0 = 0;
dt = 0.01;
tvec = t0:dt:150;


sol=ode("stiff",x0,t0,tvec, list(equations,par) );
//Faire le graphe de la solution au cour du temps, dans la figure 1
//Options: couleur (r=red, b=blue,k=black,...), forme de ligne (-,--,:,-.)
figure(1);
plot(tvec,sol(1,:),'b-',tvec,sol(2,:),'r-',tvec,sol(3,:),'g-',tvec,sol(4,:),'k-');
legend('W', 'Frq', 'Pcy', 'Pnu');

// Calcule de la période d'oscillation

 // Trace le diagramme de phase W vs Pnu
figure(2); // Crée une nouvelle figure
plot(sol(1,:), sol(4,:), 'b-'); // Trace W en fonction de Pnu
xlabel('W'); // Ajoute un label à l'axe des x
ylabel('Pnu'); // Ajoute un label à l'axe des y
title('Diagramme de phase W vs Pnu'); // Ajoute un titre au graphique

// 4.2886148     42.282118      7.5341680
// m_st,            Pc_st,      Pn_st
// W_st vaut environ 26.13 et feq 0 donc le modèle est valide (voir Projet.sce)
Pn_st = 7.5341680;


// Récupération des fonctions W et Pnu
W = sol(1,:);
Pnu = sol(4,:);

// Calcule la différence entre Pnu(t) et Pn_st. On récupère l'indice du temps si on vient de passer au dessus de Pn_st (si la soustraction change de signe)
// 1 si on passe de - à +, 0 sinon
crossings = find(diff((Pnu - Pn_st) > 0) == 1);

crossing_times = tvec(crossings); // récupère les temps correspondant aux indices

// Calcule la différence entre les temps consécutifs pour obtenir les périodes
periods = diff(crossing_times);

// On fait la moyenne des ces périodes
oscillation_period = mean(periods);

// Affiche la période d'oscillation
disp("La période d oscillation du système est: " , oscillation_period);

// Méthode de la période de la prof
Periode=calculperiode(sol(1,:),tvec);
disp(Periode);

// Méthode de la sensibilité
moyenne_S = sensibilite(par, tvec);
disp("La sensibilité est : ", moyenne_S);