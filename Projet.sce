function W_eq = W_intersection(W0, par) 
    // Fonction qui trouve les zéros de f = g1(W0) -g2(W0)
    // Définition des tous les paramètres du système
    a1 = par(1);
    a2 = par(2);
    a3 = par(3);
    a4 = par(4);
    gamma = par(5);
    gM = par(6); // Correspond à gamma M
    gP = par(7); // doit être inférieur à  1/2 * a4/K4
    n = par(8);
    K1 = par(9);
    K4 = par(10);

    // Définition des différentes fonctions pour calculer Pn*
    facteur_1bis = a3 / (gP + (gamma * W0));
    facteur_2bis = a2 / (a3 + gP);
    facteur_3bis = (a1 / gM) * ((W0^n) / ((K1^n) +(W0^n)));
    g1 = facteur_1bis * facteur_2bis *facteur_3bis; // Définition de g1
    g2 = (1/gamma) * (a4 * ((W0^(n-1) / (K4^n + W0^n)) - gP)); // Définition de g2

    W_eq = g2 - g1; // Prend le milieu pour avoir la valeur à l'équilbre de la variation de W
endfunction

function hp = hillp(x,k,nn)
    hp = x^nn /(x^nn + k^nn);
endfunction

// Main
// Définition des paramaètres:
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

para = [a1, a2, a3, a4, gamma, gM, gP, n, K1, K4];


W = [0.0:0.01:100]; // liste de W qui part de 0 par pas de 0.01 avec 100 éléments
g1 = zeros(1, length(W)); // Initialisation d'un vecteur pour stocker les valeurs de g1
g2 = zeros(1, length(W)); // Initialisation d'un vecteur pour stocker les valeurs de g2
// W_intersection_test = zeros(1, length(W)); // Initialisation d'un vecteur pour stocker les valeurs de g2

for i = 1 : length(W) // boucle pour parcourir les listes
    facteur_1bis = a3 / (gP + (gamma * W(i)));
    facteur_2bis = a2 / (a3 + gP);
    facteur_3bis = (a1 / gM) * ((W(i)^n) / ((K1^n) +(W(i)^n)));

    g2(i) = facteur_1bis * facteur_2bis *facteur_3bis; // Stocke les valeurs de g1 et de g2
    g1(i) = (1/gamma) * (a4 * ((W(i)^(n-1) / (K4^n + W(i)^n)) - gP));
    W_intersection_test(i) = W_intersection(W(i), para);
end

// Affiche le graphique de g1 et de g2
f3 = (figure(3));
plot(W, g1, "r-", W, g2, "b-", W, W_intersection_test, "g-", 'linewidth', 3);
f3.background = color("white");
xlabel("W");
legend('dW/dt = 0', 'dPnu/dt = 0', "dW/dt - dPnu/dt");

[W_st, feq] = fsolve(K4 , list(W_intersection, para)); // K4 à la place de W0 ? K4 est proche de W0
disp(W_st, feq) // W_st vaut environ 26.13 et feq 0 donc le modèle est valide

scf(3);
plot(W_st, (1/gamma) * (a4 * ((W_st^(n-1) / (K4^n + W_st^n)) - gP)), "k*" ,"linewidth", 3);
xlabel("W");

// Définition de W0 ou "initial guess" Et définition de tous les points d'équilibre 
m_st = (a1/gM) * ( W_st^n /(W_st^n + K1^n));
Pc_st = (a2 / (a3 + gP)) * m_st;
Pn_st = (a3 / (gP + gamma * W_st)) * Pc_st;

// Matrice jacobienne

A = [];
// dW/dt
A(1,1) = (a4 * n * W_st^(n-1) / (K4^n + W_st^n)^2) -gP -gamma * Pn_st; // Dérivée partielle des équations de base par rapport à W
A(1,2) = 0; // Par rapport à m
A(1,3) = 0; // Par rapport à Pc
A(1,4) = -gamma * W_st; // Par rapport à Pn
// 2e ligne matrice / 2e équation dm/dt
A(2,1) = (a1 * n * K1^n * W_st^(n-1)) / (K1^n + W_st^n)^2;
A(2,2) = -gM;
A(2,3) = 0;
A(2,4) = 0;
// dPc/dt
A(3,1) = 0;
A(3,2) = a2;
A(3,3) = - (a3 + gP);
A(3,4) = 0;
// dPn/dt
A(4,1) = -gamma * Pn_st;
A(4,2) = 0;
A(4,3) = a3;
A(4,4) = -gP -gamma * W_st;

// Calcul des valeurs propres partie réelle négative --> poit équilibre stables
vp = spec(A);

// Calcul matrice jacobienne à W_st = Pc_st = Pn_st = m_st = 0
A0 = [];

A0(1,1) = -gP;
A0(1,2) = 0; 
A0(1,3) = 0; 
A0(1,4) = 0;
// 2e ligne matrice / 2e équation dm/dt
A0(2,1) = 0;
A0(2,2) = -gM;
A0(2,3) = 0;
A0(2,4) = 0;
// dPc/dt
A0(3,1) = 0;
A0(3,2) = a2;
A0(3,3) = - (a3 + gP);
A0(3,4) = 0;
// dPn/dt
A0(4,1) = 0;
A0(4,2) = 0;
A0(4,3) = a3;
A0(4,4) = -gP;

vp0 = spec(A0);
disp(vp0);  // Etat initial stable car toutes les parties réelles sont négatives
disp(vp);
disp(m_st, Pc_st, Pn_st); //   4.2886148     42.282118      7.5341680