
# Modélisation du cycle circadien de *Neurospora Crassa*

## Présentation
   Ce projet vise à modéliser l'évolution des concetrations des deux grandes protéines régulant le cycle circadien du champignon ***Neurospora Crassa***. Nous avons sélectionné la protéine WCC et FRQ, pour FRQ nous avons différencié la protéine cytosolique et nucléaire. L'ARNm de FRQ a aussi été modélisé.

## Outils
- Scilab deskstop pour l'analyse du système différentiel et pour certains graphiques simples
- Python 3.12 pour l'optimisation des paramètres du système et pour les graphiques complexes

## Modèle

$ \frac{dW}{dt} = a_4 W_n K_{n4} + W_n - \gamma_P W - \gamma_{WP} W_{nu} $

$ \frac{dm}{dt} = a_1 W_n K_{n1} + W_n - \gamma_{mm} $

$ \frac{dP_{cy}}{dt} = a_2 m - a_3 P_{cy} - \gamma_P P_{cy} $

$ \frac{dP_{nu}}{dt} = a_3 P_{cy} - \gamma_P P_{nu} - \gamma_{WP} P_{nu} $

## Auteurs
Edmond BERNE
- [@Edmondbrn](https://github.com/Edmondbrn)
Julie FERIAU


## Screenshots

![App Screenshot](https://via.placeholder.com/468x300?text=App+Screenshot+Here)

