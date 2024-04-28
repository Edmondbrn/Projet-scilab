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

    a1 = pr(1);
    a2 = pr(2);
    a3 = pr(3);
    a4 = pr(4);
    gamma = pr(5);
    gM = pr(6); // Correspond à gamma M
    gP = pr(7); // doit être inférieur à  1/2 * a4/K4

    K1 = pr(8);
    K4 = pr(9);

    n = pr(10);
    deg = gamma * Pnu * W;
    
    // Afficher les valeurs des variables
    // disp('W = ', W);
    // disp('Frq = ', Frq);
    // disp('Pcy = ', Pcy);
    // disp('Pnu = ', Pnu);
    // disp('a1 = ', a1);
    // disp('a2 = ', a2);
    // disp('a3 = ', a3);
    // disp('a4 = ', a4);
    // disp('gamma = ', gamma);
    // disp('gM = ', gM);
    // disp('gP = ', gP);
    // disp('K1 = ', K1);
    // disp('K4 = ', K4);
    // disp('n = ', n);
    // disp('deg = ', deg);

    //Initialisation du vecteur des derivees
    v = [];
    v(1) = a4 * hillp(W, K4, n) - deg - gP * W;
    v(2) = a1* hillp(W, K1, n) -gM * Frq ;
    v(3)= a2 * Frq -a3 * Pcy - gP * Pcy;
    v(4) = a3 * Pcy -gP * Pnu -deg ;
    
    // Afficher la valeur de v
    // disp('v = ', v);
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
        parS(j) = par(j) * 1.05; // on augmente le paramètre j de 5%
        xsolS = ode("stiff", X0, t0, tvec, list(equations, parS)); // On résout le système avec les nouveaux paramètres
        Periode_S = calculperiode(xsolS(1,:), tvec); // ON calcule la période d'oscillation du nouveau système
        S(j) = (Periode_S - periode) / (parS(j) - par(j)); // On calcule la sensibilité
        // disp("La sensibilité du paramètre", Periode_S, periode);
    end
    moyenne_S = mean(S.*S); // On calcule la moyenne des sensibilités
    // S = sqrt(moyenne_S); // On retourne l'écart type
    for j = 1:10
        Sj_S = S(j) / moyenne_S;
        disp("La sensibilité du paramètre ", j, " est : ", Sj_S);
    end
endfunction


function amp = amplitude(par)
    X0 = [1;0.1;0.1;0.1];
    t0 = 0;
    Nc = 30; // Nombre de cycles
    tvec = [t0:0.05 : Nc * 24];
    L = length(tvec);
    L3 = floor(L * 0.3);
    xsol = ode("stiff", X0, t0, tvec, list(equations, par)); // resolution du système
    AmaxW = max(xsol(1,L3:L)) / max(xsol(1,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    AminW = min(xsol(1,L3:L)) / max(xsol(1,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    ampW = AmaxW - AminW;

    AmaxF = max(xsol(2,L3:L)) / max(xsol(2,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    AminF = min(xsol(2,L3:L))/ max(xsol(2,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    ampF = AmaxF - AminF;

    AmaxPc = max(xsol(3,L3:L))/ max(xsol(3,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    AminPc = min(xsol(3,L3:L))/ max(xsol(3,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    ampPc = AmaxPc - AminPc;

    AmaxPn = max(xsol(4,L3:L))/ max(xsol(4,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    AminPn = min(xsol(4,L3:L))/ max(xsol(4,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    ampPn = AmaxPn - AminPn;
    // vecteur pour stocker les amplitudes
    amp = [ampW, ampF, ampPc, ampPn];
endfunction

function Periode= calculperiode(X,tvec)
    //Initialisation
    Periode = 0 ;
    Pvec = []; // Ajoutez cette ligne pour initialiser Pvec
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
    if length(Pvec) > 1 then
        Periode=mean(diff(Pvec));
        Pstd = stdev(diff(Pvec)); // vérifier ecart type entre les temps
        disp("L ecart type entre les différents temps est : ", Pstd); 
    else
        disp("Pvec ne contient pas assez d éléments pour calculer la période et l écart-type.");
        Periode = 0/0; // Retourne une valeur par défaut
        Pstd = 0/0; // Retourne une valeur par défaut
    end
    endfunction

function oscillation_period = calculate_oscillation_period(sol, Pn_st, tvec)
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
endfunction


function c = fonction_cout(par)
    X0 = [1;0.1;0.1;0.1];
    t0 = 0;
    Nc = 30; // Nombre de cycles
    tvec = [t0:0.05 : Nc * 24];
    L = length(tvec);
    L3 = floor(L * 0.3);
    xsol = ode("stiff", X0, t0, tvec, list(equations, par)); // resolution du système

    periodeS = calculperiode(xsol(1,:), tvec); // calcul de la période d'oscillation,:
    periodeS_obs = 24; // période observée
    // c = (periodeS - periodeS_obs)^2;

    AmaxS = max(xsol(1,L3:L)) / max(xsol(1,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    AminS = min(xsol(1,L3:L)) / max(xsol(1,L3:L)); // calcul de l'amplitude on normalise par rapport à la valeur max
    Amax_obs = 14/16; // amplitude observée maximale
    Amin_obs = 1/16; // amplitude observée minimale

    //FRQ total
    Ftotal = xsol(2,L3:L);
    Ampmin_FRQ = min(Ftotal(L3:L)) / max(Ftotal(L3:L));
    Amin_FRQ_obs = 1/5;
    c = (periodeS - periodeS_obs)^2 + (AmaxS - Amax_obs)^2 + (AminS - Amin_obs)^2 + (Ampmin_FRQ - Amin_FRQ_obs)^2;
    endfunction


function [liste_periode, liste_a3] = modif_a3(par)
    liste_periode = [];
    liste_a3 = [];
    while par(3) < 10
        X0 = [1; 5; 5; 1]; // Conditions initiales
        t0 = 0;
        dt = 0.01;
        tvec = t0:dt:150;
        par(3) = par(3) * 1.2;
        xsol = ode("stiff", X0, t0, tvec, list(equations, par)); // On résout le système avec les nouveaux paramètres
        liste_periode = [liste_periode, calculperiode(xsol(1,:), tvec)]; // On calcule la période d'oscillation du nouveau système
        liste_a3 = [liste_a3, par(3)];
        disp(par(3));
    end
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

x0 = [1;5;5;1]; // Conditions initiales
t0 = 0;
dt = 0.01;
tvec = t0:dt:150;


sol=ode("stiff",x0,t0,tvec, list(equations,par) );
// //Faire le graphe de la solution au cour du temps, dans la figure 1
// //Options: couleur (r=red, b=blue,k=black,...), forme de ligne (-,--,:,-.)
// f1 = (figure(1));
// plot(tvec,sol(1,:),'b-',tvec,sol(2,:),'r-',tvec,sol(3,:),'g-',tvec,sol(4,:),'k-');
// legend('W', 'Frq', 'Pcy', 'Pnu');
// f1.background = color("white");
// xlabel('Temps');
// ylabel("Quantité");
// set(gca(), 'font_size', 4); // Changer la taille de la police à 4

// // Calcul de la période d'oscillation

//  // Trace le diagramme de phase W vs Pnu
// f2 = (figure(2)); // Crée une nouvelle figure
// plot(sol(1,:), sol(4,:), 'b-'); // Trace W en fonction de Pnu
// xlabel('W'); // Ajoute un label à l'axe des x
// ylabel('Pnu'); // Ajoute un label à l'axe des y
// title('Diagramme de phase W vs Pnu'); // Ajoute un titre au graphique
// f2.background = color("white"); // Change la couleur de fond de la figure
// xgrid(); // Ajoute une grille sur l'axe des x
// set(gca(), 'font_size', 4); // Changer la taille de la police à 4


// // 4.2886148     42.282118      7.5341680
// // m_st,            Pc_st,      Pn_st
// // W_st vaut environ 26.13 et feq 0 donc le modèle est valide (voir Projet.sce)
Pn_st = 7.5341680;



// Affiche la période d'oscillation
// disp("La période d oscillation du système est: " , oscillation_period);

// // Méthode de la période de la prof
Periode=calculperiode(sol(1,:),tvec);
disp(Periode);

// // Méthode de la sensibilité
// moyenne_S = sensibilite(par, tvec);
// disp("La sensibilité est : ", moyenne_S);

// disp("les amplitudes sont:" , amplitude(par));


// Paramètres initiaux
par_init = [8, 7, 0.5, 40, 0.15, 0.4, 0.01, 50, 0.7, 2]; // Remplacez par vos valeurs initiales
// Estimation des paramètres
par_est = fminsearch(fonction_cout, par_init);
disp("Les paramètres estimés sont : ", par_est);

// [periodes, a3s] = modif_a3(par_init);
// f3 = figure(3);
// plot(a3s, periodes, "ro-", 'LineWidth', 2); // Tracer la période en fonction de a3 avec une ligne rouge et une épaisseur de ligne de 2
// xtitle('Evolution de la période en fonction de a3', 'a3', 'Période');
// set(gca(), 'font_size', 4); // Changer la taille de la police à 4
// xgrid(); // Ajouter une grille sur l'axe des x
// f3.background = color("white"); // Change la couleur de fond de la figure


// dataFRQ = [];
// dataFrq = []; // données à prendre sur les données de départ
// dataTemps = []
// dt = 0.01;
// tvec = [0:dt:200];
// xsol = ode("stiff", x0, t0, tvec, list(equations, par));
// // Données et modèles, amplitudes normalisées

// Ptotal_modele = xsol(3, :) + xsol(4, :);
// PnTot = Ptotal_modele / max(Ptotal_modele);
// dataFRQn = dataFRQ / max(dataFRQ);
// data_frq = dataFrq / max(dataFrq);
// frq_modele = xsol(2, :) / max(xsol(2, :));

// //Identifier 2 cycles du modèles a partir d'un maximum de frq_modele
// L = length(frq_modele);
// L3 = floor(0.25 * L);
// [maxfrq, indicemaxfrq] = max(frq_modele(L3:L));

// i2 = indicemaxfrq + 48/dt;
// indices_T = [indicemaxfrq : 1 : i2];

// offset = tvec(indicemaxfrq) - dataTemps(1);



// figure(4);
// plot(tvec(indices_T) + offset, frq_modele(indices_T), 'r-');
// plot(dataTemps, data_frq, 'b.'); //  on veut comparerr le modèle aux données théoriques
