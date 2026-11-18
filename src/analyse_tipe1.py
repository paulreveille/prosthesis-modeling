# -*- coding: utf-8 -*-
"""
TIPE - Prothèse tibiale en fibre de carbone
Analyse de l'essai en compression n°1 : ligne de charge à 10 mm du bord
de la lame (correspond à un angle de 20°).

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
DONNEES = os.path.join(DOSSIER, "..", "data", "test_1_20deg_10mm")
RESULTATS = os.path.join(DOSSIER, "..", "results")
os.makedirs(RESULTATS, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Fichier brut (non corrigé) : 4 colonnes
#    temps, déplacement vérin, déplacement remis à zéro, effort vérin
# ---------------------------------------------------------------------------
with open(os.path.join(DONNEES, "raw.txt")) as fichier:
    lignes = fichier.readlines()[2:]  # on saute les 2 lignes d'en-tête

AB = []  # déplacement (mm), fichier brut
EB = []  # effort (N), fichier brut
for ligne in lignes:
    x1, x2, x3, x4 = ligne.split()
    AB.append(-float(x3))
    EB.append(-float(x4))


# ---------------------------------------------------------------------------
# 2. Régression polynomiale par moindres carrés : F = a x^2 + b x + c
#    sur les 20 derniers points. Résolution de Y = A (a, b, c)^T
#    avec (a, b, c) = (A^T A)^-1 A^T Y
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 3. Fichier corrigé : 2 colonnes (déplacement remis à zéro, effort)
# ---------------------------------------------------------------------------
with open(os.path.join(DONNEES, "corrected.txt")) as fichier:
    lignes = fichier.readlines()[2:]

A1 = []  # déplacement (mm)
E1 = []  # effort (N)
for ligne in lignes:
    x1, x2 = ligne.split()
    A1.append(-float(x1))
    E1.append(-float(x2))


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


Emax = maximum(E1)[0]
imax = maximum(E1)[1]
Amax = maximum(A1)[0]

# Tracé de la ligne de charge
plt.figure()
plt.title("Essai en compression")
plt.grid()
plt.plot(A1, E1, label="LDC 10mm", color="deepskyblue")
plt.xlabel("Déplacement en mm")
plt.ylabel("Effort en N")
plt.legend(loc="best")
plt.savefig(os.path.join(RESULTATS, "load_line_10mm.png"), dpi=150)
plt.show()


# ---------------------------------------------------------------------------
# 4. Intégration (méthode des trapèzes) pour déterminer l'énergie perdue
#    AireH : aire sous la courbe à la montée (énergie stockée, en J)
#    AireB : aire sous la courbe à la descente (négative : énergie restituée)
# ---------------------------------------------------------------------------
N = len(A1)
AireH = 0
AireB = 0
for k in range(imax + 1):
    AireH += ((A1[k + 1] - A1[k]) * (10 ** -3)) * (E1[k + 1] + E1[k]) / 2
for k in range(imax + 1, N - 1):
    AireB += ((A1[k + 1] - A1[k]) * (10 ** -3)) * (E1[k + 1] + E1[k]) / 2
Ep = AireH + AireB  # énergie perdue (hystérésis)


def raideur(E, A, Emax):
    """Moyenne des rapports effort / déplacement jusqu'à l'effort maximal."""
    rm = 0
    i = 0
    while E[i] != Emax:
        rm += abs(E[i]) / abs(A[i])
        i += 1
    return rm / i


k = raideur(E1, A1, Emax) * 10 ** 3  # raideur en N/m

print("Effort max : %.1f N pour %.2f mm" % (Emax, Amax))
print("Energie stockée   : %.2f J" % AireH)
print("Energie restituée : %.2f J" % (-AireB))
print("Energie perdue    : %.2f J" % Ep)
print("Raideur           : %.0f N/m" % k)
