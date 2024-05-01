
# Modélisation du cycle circadien de *Neurospora Crassa*

## Présentation
   Ce projet vise à modéliser l'évolution des concetrations des deux grandes protéines régulant le cycle circadien du champignon ***Neurospora Crassa***. Nous avons sélectionné la protéine WCC et FRQ, pour FRQ nous avons différencié la protéine cytosolique et nucléaire. L'ARNm de FRQ a aussi été modélisé.

## Outils
- Scilab deskstop pour l'analyse du système différentiel et pour certains graphiques simples
- Python 3.12 pour l'optimisation des paramètres du système et pour les graphiques complexes

## Modèle
 $$\frac{dW}{dt} = \frac{a4 \cdot Wn}{Kn4 + Wn} - \gamma_P \cdot W - \gamma_W \cdot P_{nu}$$

 $$\frac{dm}{dt} = \frac{a1 \cdot Wn}{Kn1 + Wn} - \gamma_m \cdot m$$

 $$\frac{dP_{cy}}{dt} = a2 \cdot m - a3 \cdot P_{cy} - \gamma_P \cdot P_{cy}$$

 $$\frac{dP_{nu}}{dt} = a3 \cdot P_{cy} - \gamma_P \cdot P_{nu} - \gamma_W \cdot P_{nu}$$

### Intervalle conseillé pour les paramètres:
- a1,a2 ∈[5,10]
- a3 ∈[0.5,1]
- a4 ∈[20,40]
-  K1 ∈[20,50]
- K4 ∈[0.5,1]
- γ ∈[0.1,0.2]
-  γm = 0.4
-  γP = 0.01
-  n = 2

## Analyse
### Etude des équilibres par la matrice jacobienne linéarisée
- Equilibre nul
Valeurs propres

| Partie réelle | Partie imaginaire |
| --------- | --------- |
|-0.01   | 0 |   
| -0.01   | 0  |
| -0.71   | 0 | 
| -0.4   | 0  |

Il s'agit d'un équilibre stable.

- Equilibre d'intérêt

![Logo]("Intersection_W_st.png")

1. W_eq = 26.13
2. m_eq = 4.29
3. Pc_eq = 42.28
4. Pn_eq = 7.53

| Partie réelle | Partie imaginaire |
| --------- | --------- |
|-5.04   | 0 |   
| 0.04  | 0.57i  |
| 0.04   | -0.57i | 
| -1.22  | 0  |

Il s'agit d'un équilibre instable.

## Comportement du système
### Capacité oscillatoire
Notre système doit pouvoir osciller puisque nous voulons modéliser un cycle jour/nuit.

![Logo]("graph/oscillation_pre_optimisation.png")
## Période d'oscillation

Nous avons calculé la période d'oscillation du système de deux manières différentes:
- En se basant sur le portrait de phase de W et de P_nu (FRQ nucléaire), nous avons fixé le point de coordonnées (W_eq, Pn_eq). Nous avons ensuite extrait les intervalles de temps entre chaque passage de ce point et calculer la moyenne de ce dernier. Nous avons exclus les premiers temps (30%), afin d'éviter les biais par la mise en place de l'oscillation.
![App Screenshot]("graph/Portrait_de_phase_W_Pnu.png")

- La deuxième méthode se base sur une chronique de W en fonction du temps. Nous avons fixé un point (0.5*amplitude maxiamle) et nous avons également calculer les différents intervalles de temps pour l'atteindre. Nous avons aussi exclus les 30% premiers temps pour éviter les biais.

Ces deux méthodes nous donnent des résultats similaires avec un période d'environ 12.2 avant optimisation des paramètres (période théorique ~22)

## Sensibilité des paramètres du système

Nous avons voulu ensuite déterminer les paramètres les plus sensibles de notre système. Pour ce faire, nous les avons faits varier de +/- 5% et nous avons ensuite calculer la différence entre la période initiale et la nouvelle période.

![App Screenshot]("graph/Sensibilite_parametre.png")

Les paramètres les plus intéressants sont donc :
-γm et γp qui représentent les taux de dégradation protéiques et de l'ARNm de FRQ
- a3 qui est le taux de transport cyrosolique-nucléaire de FRQ

## Analyse de l'impact de l'évolution de a3

![App Screenshot]("graph/Effet_variation_a3_periode.png")

Par la suite, nous avons voulu quantifier l'impact de l'évolution de a3 sur la période du système. Nous pouvons voir que a3 influe énormémement sur la période lorsqu'il est petit (entre 0 et 1). Cela confirme l'intervalle de départ et sera utile lors de l'optimisation des paramètres du système.

## Optimisation des paramètres
Pour ce faire, nous avons utilisé une fonction coût qui prend en compte :
- La période
- L'amplitude maximale normalisée
- L'amplitude minimale normalisée
- L'amplitude minimale de FRQ normalisée

Après optimisation (valeur fonction coût ~8), nous avons obtenu ce graphique:

![App Screenshot]("graph/Oscillation_post_optimisation.png")

## Comparaison du modèle et des données observées

Nous avons ensuite comparé notre modèle aux données expérimentales:
![App Screenshot]("graph/Graphique_final.png")

Les points de départ de chaque données ont été décalés afin qu'ils partent du même endroit. 
## Auteurs
Edmond BERNE
- [@Edmondbrn](https://github.com/Edmondbrn)

Julie FERIAU

