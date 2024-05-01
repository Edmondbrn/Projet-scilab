# import seaborn as sns
# import numpy as np
# import matplotlib.pyplot as plt

# par = ["a1", "a2", "a3", "a4", "gamma", "gM", "gP", "K1", "K4", "n"]


# data = np.array([0.2107143, 0.2367347, -8.2193878, -0.0271429, -3.4285714, -18.071429, -37.142857, -0.0217143, 0.2857143, 0.0285714])
# S = np.sqrt(np.mean(data**2))  # Calculer S

# normalized_data = data / S

# sns.set(style="darkgrid")
# fig, ax = plt.subplots(figsize=(10, 5))
# sns.lineplot(x=par, y=normalized_data, marker="o", color="b")
# plt.xlabel("Paramètres")
# plt.ylabel("Sensibilités normalisées")
# ax.set_title("Sensibilité de la période du système en fonction des paramètres", fontweight="bold", fontsize=15, pad=20)  # Décaler le titre vers le haut
# # Ajouter les valeurs des points sur le graphique
# for i, value in enumerate(normalized_data):
#     plt.text(i + 0.08, value + 0.01, f"{value:.2f}")

# plt.show()

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import differential_evolution
import seaborn as sns
import matplotlib.pyplot as plt

def hillp(x, K, n):
    return x**n / (K**n + x**n)

# Votre fonction d'équations
def equations(x, t, pr):
    W, Frq, Pcy, Pnu = x
    a1, a2, a3, a4, gamma, gM, gP, K1, K4, n = pr
    deg = gamma * Pnu * W
    v = np.zeros(4)
    v[0] = a4 * hillp(W, K4, n) - deg - gP * W
    v[1] = a1 * hillp(W, K1, n) - gM * Frq
    v[2] = a2 * Frq - a3 * Pcy - gP * Pcy
    v[3] = a3 * Pcy - gP * Pnu - deg
    return v

def calculperiode(X, tvec):
    L = len(X)
    S = int(0.3 * L)
    MaxX = max(X[S:L])
    MinX = min(X[S:L])
    theta = 0.5 * (MaxX + MinX)
    Pvec = []
    for i in range(S, L - 1):
        auxg = X[i] - theta
        auxd = X[i + 1] - theta
        if auxg < 0 and auxd > 0:
            Pvec.append(tvec[i])
    if len(Pvec) > 1:
        Periode = np.mean(np.diff(Pvec))
        Pstd = np.std(np.diff(Pvec))  # vérifier écart type entre les temps
    else:
        print("Pvec ne contient pas assez d'éléments pour calculer la période et l'écart-type.")
        Periode = np.nan  # Retourne une valeur par défaut
        Pstd = np.nan  # Retourne une valeur par défaut
    return Periode


# Votre fonction de coût
def fonction_cout(pr):
    X0 = [5, 1, 1, 5]
    t0 = 0
    Nc = 30  # Nombre de cycles
    tvec = np.arange(t0, Nc * 24, 0.05)
    L = len(tvec)
    L3 = int(L * 0.3)
    xsol = odeint(equations, X0, tvec, args=(pr,), rtol=1e-5, atol=1e-8)
    periodeS = calculperiode(xsol[:, 0], tvec)  # calcul de la période d'oscillation
    periodeS_obs = 22  # période observée

    AmaxS = max(xsol[L3:L, 0]) / max(xsol[L3:L, 0])  # calcul de l'amplitude on normalise par rapport à la valeur max (W)
    AminS = min(xsol[L3:L, 0]) / max(xsol[L3:L, 0])  # calcul de l'amplitude on normalise par rapport à la valeur max
    Amax_obs = 14 / 16  # amplitude observée maximale
    Amin_obs = 1 / 16  # amplitude observée minimale

    # FRQ total
    Ftotal = xsol[L3:L, 1]  #+ xsol[L3:L, 2]
    Ampmin_FRQ = min(Ftotal) / max(Ftotal)
    Amin_FRQ_obs = 1 / 5
    Amax_FRQ_obs = 5/5

    c = (periodeS - periodeS_obs)**2 + (AmaxS - Amax_obs)**2 + (AminS - Amin_obs)**2 + (Ampmin_FRQ - Amin_FRQ_obs)**2 
    if c < 10:  # Remplacez 10 par votre seuil
        print("Seuil atteint")
        print(pr)
        print(periodeS)
        print(c)
        sns.lineplot(xsol)
        plt.show()
        quit()
    # sns.lineplot(xsol)
    # plt.show()
    print(c)
    return c

def graphique(pr):
    X0 = [5, 1, 1, 5]
    t0 = 0
    Nc = 30  # Nombre de cycles
    tvec = np.arange(t0, Nc * 24, 0.05)
    L = len(tvec)
    L3 = int(L * 0.3)
    xsol = odeint(equations, X0, tvec, args=(pr,), rtol=1e-5, atol=1e-8)
    sns.lineplot(xsol)
    return plt.show()


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def graph(par):
    sns.set_style("darkgrid")  # Définit le style darkgrid 

    x0 = [5, 1, 1, 5]  # Conditions initiales
    t0 = 9 # petit décalage pour bien avoir les courbes superposées
    t0W = 0
    dt = 0.01
    tvec = np.arange(t0, 150, dt)
    tvecW = np.arange(t0W, 150, dt)
    # Supposons que vous ayez déjà défini les équations et les paramètres
    Xnouveau = odeint(equations, x0, tvec, args=(par,))
    XnouveauW = odeint(equations, x0, tvecW, args=(par,))
    dataFRQ = [1.5, 5.5, 9, 5, 4, 0.5, 3, 7, 5, 4, 2.5, 1]
    data_feq = [9, 4, 2, 0.5, 0.5, 9, 3, 3, 2, 1, 1, 0.5]
    data_W = [10,15,10,2,2,6,6,4,2,2,1,12]
    temps = [11,15,18,22,25,29,33,36,40,43,47,51]
    temps_W = [0,4,8,12,16,20,24,28,32,36,40,44]

    Ptotal = Xnouveau[:, 2] + Xnouveau[:, 3]
    P_tot_norma = Ptotal / np.max(Ptotal)
    dataFRQ_norma = np.array(dataFRQ) / np.max(dataFRQ)
    datafeq_norma = np.array(data_feq) / np.max(data_feq)
    data_W_norma = np.array(data_W) / np.max(data_W)
    frq_mod = Xnouveau[:, 1] / np.max(Xnouveau[:, 1])
    W_norma = XnouveauW[:,0] / np.max(XnouveauW[:,0])

    _, axs = plt.subplots(3, 1, figsize=(10, 18))  # Crée une grille de 3 lignes et 1 colonne

    offset = 200
    fin = 5050
    axs[0].plot(tvec[offset:fin], P_tot_norma[offset:fin], 'blue', label='FRQ total')
    axs[0].plot(temps, dataFRQ_norma, "go", label='Data FRQ')
    axs[0].set_xlim(11, 55)
    axs[0].set_xlabel('Temps', fontsize = 12)
    axs[0].set_ylabel("Quantité normalisée", fontsize = 12)
    axs[0].set_title("Comparaison entre le modèle et les données théoriques", fontweight = "bold", fontsize = 15)
    axs[0].legend(loc = "upper right")
    for i, txt in enumerate(dataFRQ):
        axs[0].text(temps[i] + 0.5, dataFRQ_norma[i], str(round(dataFRQ_norma[i], 2)), fontsize=8)

    axs[1].plot(tvec[offset:fin], frq_mod[offset:fin], 'red', label='ARNm Frq')
    axs[1].plot(temps, datafeq_norma, "bo", label='Data ARNm')
    axs[1].set_xlim(11, 55)
    axs[1].set_xlabel('Temps', fontsize = 12)
    axs[1].set_ylabel("Quantité normalisée", fontsize = 12)
    axs[1].legend(loc = "upper right")
    for i, txt in enumerate(datafeq_norma):
        axs[1].text(temps[i] +0.5, datafeq_norma[i], str(round(datafeq_norma[i],2)), fontsize=8)

    axs[2].plot(tvecW[0:fin], W_norma[0:fin], 'black', label='Protéine WCC')
    axs[2].plot(temps_W, data_W_norma, "ro", label='Data WCC')
    axs[2].set_xlim(0, 55)
    axs[2].set_xlabel('Temps', fontsize = 12)
    axs[2].set_ylabel("Quantité normalisée", fontsize = 12)
    axs[2].legend(loc = "upper right")
    for i, txt in enumerate(data_W_norma):
        axs[2].text(temps_W[i] +0.5, data_W_norma[i], str(round(data_W_norma[i],2)), fontsize=8)

    plt.show()




# Limites pour les paramètres
bounds = [(5, 10), (5, 10), (0.5, 1), (20, 40), (0.1, 0.2), (0.36, 0.44), (0.01, 0.01), (20, 50), (0.5, 1), (2,2)]  # Remplacez par vos limites
# para = [7.98997769e+00, 8.29173119e+00, 5.38781314e-01, 2.20830606e+01,1.74807138e-01, 3.81164483e-01, 1.06649477e-02, 3.22598465e+01,5.31703912e-01, 1.92813431e+00]
para = [6.85776404e+00 ,7.96169636e+00, 5.16115481e-01, 2.39993029e+01,
 1.53461874e-01, 3.81329247e-01, 1.00000000e-02, 2.73160696e+01,
 6.27279798e-01, 2.00000000e+00]

graph(para)
# graph([8.84551167e+00, 9.11727604e+00, 5.13363343e-01, 2.17881216e+01,1.52438272e-01, 3.63689338e-01, 1.03896957e-02, 4.79744318e+01,6.06026462e-01, 2.08140361e+00])

# Optimisation globale
# result = differential_evolution(fonction_cout, bounds)
