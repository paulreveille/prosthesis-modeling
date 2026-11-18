# -*- coding: utf-8 -*-
"""
TIPE - Prothèse tibiale en fibre de carbone
Analyse de l'essai en compression n°3 : ligne de charge à 6 mm du bord
de la lame (correspond à un angle de 14°).

Ce script :
  1. lit les mesures du banc d'essai (déplacement en mm, effort en N),
  2. ajuste un polynôme de degré 2 par moindres carrés sur les 20 derniers
     points du fichier brut (correction du décalage de zéro),
  3. trace la ligne de charge (effort en fonction du déplacement),
  4. calcule l'énergie stockée, l'énergie restituée et la raideur.

Auteur : Paul Réveille
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression

# Dossiers du dépôt (le script peut être lancé depuis n'importe où)
DOSSIER = os.path.dirname(os.path.abspath(__file__))
DONNEES = os.path.join(DOSSIER, "..", "data", "test_3_14deg_6mm")
RESULTATS = os.path.join(DOSSIER, "..", "results")
os.makedirs(RESULTATS, exist_ok=True)


with open(os.path.join(DONNEES, "raw.txt")) as fichier:
    lignes = fichier.readlines()[2:]  # on saute les 2 lignes d'en-tête

AB = []  # déplacement (mm), fichier brut
EB = []  # effort (N), fichier brut
for ligne in lignes:
    x1, x2, x3, x4 = ligne.split()
    AB.append(-float(x3))
    EB.append(-float(x4))


X = AB[-20:]
Y = EB[-20:]
X2 = [X[k] ** 2 for k in range(len(X))]

# Matrice A du modèle : colonnes x^2, x, 1
A = np.vstack([X2, X, [1 for _ in range(20)]]).T

modele = LinearRegression(fit_intercept=False)
modele.fit(A, Y)
coeffs = modele.coef_
r2 = modele.score(A, Y)
print("Régression : a = %.4f, b = %.4f, c = %.4f, R2 = %.5f"
      % (coeffs[0], coeffs[1], coeffs[2], r2))


with open(os.path.join(DONNEES, "corrected.txt")) as fichier:
    lignes = fichier.readlines()[2:]

A3 = []  # déplacement (mm)
E3 = []  # effort (N)
for ligne in lignes:
    x1, x2 = ligne.split()
    A3.append(-float(x1))
    E3.append(-float(x2))


def maximum(L):
    """Renvoie [valeur maximale de L, indice de ce maximum]."""
    M = L[0]
    i = 0
    n = len(L)
    for k in range(n):
        if L[k] >= M:
            M = L[k]
            i = k
    return [M, i]


Emax = maximum(E3)[0]
imax = maximum(E3)[1]
Amax = maximum(A3)[0]

# Tracé de la ligne de charge
plt.figure()
plt.title("Essai en compression")
plt.grid()
plt.plot(A3, E3, label="LDC 6mm", color="darkorange")
plt.xlabel("Déplacement en mm")
plt.ylabel("Effort en N")
plt.legend(loc="best")
plt.savefig(os.path.join(RESULTATS, "load_line_6mm.png"), dpi=150)
plt.show()


N = len(A3)
AireH = 0
AireB = 0
for k in range(imax + 1):
    AireH += ((A3[k + 1] - A3[k]) * (10 ** -3)) * (E3[k + 1] + E3[k]) / 2
for k in range(imax + 1, N - 1):
    AireB += ((A3[k + 1] - A3[k]) * (10 ** -3)) * (E3[k + 1] + E3[k]) / 2
Ep = AireH + AireB  # énergie perdue (hystérésis)


def raideur(E, A, Emax):
    """Moyenne des rapports effort / déplacement jusqu'à l'effort maximal."""
    rm = 0
    i = 0
    while E[i] != Emax:
        rm += abs(E[i]) / abs(A[i])
        i += 1
    return rm / i


k = raideur(E3, A3, Emax) * 10 ** 3  # raideur en N/m

print("Effort max : %.1f N pour %.2f mm" % (Emax, Amax))
print("Energie stockée   : %.2f J" % AireH)
print("Energie restituée : %.2f J" % (-AireB))
print("Energie perdue    : %.2f J" % Ep)
print("Raideur           : %.0f N/m" % k)
