
# Modélisation du cycle circadien de *Neurospora Crassa*

## Présentation
   Ce projet vise à modéliser l'évolution des concetrations des deux grandes protéines régulant le cycle circadien du champignon ***Neurospora Crassa***. Nous avons sélectionné la protéine WCC et FRQ, pour FRQ nous avons différencié la protéine cytosolique et nucléaire. L'ARNm de FRQ a aussi été modélisé.

## Outils
- Scilab deskstop pour l'analyse du système différentiel et pour certains graphiques simples
- Python 3.12 pour l'optimisation des paramètres du système et pour les graphiques complexes

## Modèle

Les équations sont les suivantes :

1. $$\frac{dW}{dt} = \frac{a4 \cdot Wn}{Kn4 + Wn} - \gamma_P \cdot W - \gamma_W \cdot P_{nu}$$

2. $$\frac{dm}{dt} = \frac{a1 \cdot Wn}{Kn1 + Wn} - \gamma_m \cdot m$$

3. $$\frac{dP_{cy}}{dt} = a2 \cdot m - a3 \cdot P_{cy} - \gamma_P \cdot P_{cy}$$

4. $$\frac{dP_{nu}}{dt} = a3 \cdot P_{cy} - \gamma_P \cdot P_{nu} - \gamma_W \cdot P_{nu}$$


## Auteurs
Edmond BERNE
- [@Edmondbrn](https://github.com/Edmondbrn)
Julie FERIAU


## Screenshots

![App Screenshot](https://via.placeholder.com/468x300?text=App+Screenshot+Here)

