function W_eq = W_intersection(W0, par)
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
    facteur_1 = a3 / (gP + (gamma * W0));
    facteur_2 = a2 / (a3 + gP);
    facteur_3 = (a1 / gM) * ((W^n) / ((K1^n) +(W^n)));
    g1 = facteur_1 * facteur_2 *facteur_3;

    g2 = (1/gamma) * (a4 * ((W^(n-1) / (K4^n + W0^n)) - gP)); // Définition de g2
    W_eq = g1 - g2; // Prend le milieu pour avoir la valeur à l'équilbre de la variation de W
endfunction

// Main
// Définition des paramaètres:
// a1 =
// a2 =
// a3 =
// a4 =
// gamma = 
// gM = 
// gP =
// n = 
// K1 = 
// K4 =

// para = [a1, a2, a3, a4, gamma, gM, gP, n, K1, K4];


// Définition de W0 ou "initial guess" Et définition de tous les points d'équilibre 
W0 = K4 * 4;
[W_st, feq] = fsolve(W0, list(W_intersection, para));
m_st = (a1/gM) / (W^n + K1^n);
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
A(4,1) = -gamma*Pn_st;
A(4,2) = 0;
A(4,3) = a3;
A(4,4) = -gP -gamma * W_st;

// Calcul des valeurs propres
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






