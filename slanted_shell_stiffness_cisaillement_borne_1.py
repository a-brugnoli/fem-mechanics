import numpy as np

def slanted_shell_stiffness_cisaillement_borne_1(n, E1, E3, v13, v31, rm, t, L, phi, r2):
    coef = 5 / 6  # coeff de correction cisaillement
    r = rm + r2 * n

    # Définition des composantes de la matrice B
    B = np.zeros((5, 6))

    B[0, 0] = -1 / L
    B[0, 3] = 1 / L

    B[1, 0] = ((1 - n) * np.sin(phi)) / (2 * r)
    B[1, 1] = ((1 - n) * np.cos(phi)) / (2 * r)
    B[1, 3] = ((1 + n) * np.sin(phi)) / (2 * r)
    B[1, 4] = ((1 + n) * np.cos(phi)) / (2 * r)

    B[2, 2] = 1 / L
    B[2, 5] = -1 / L

    B[3, 2] = -((1 - n) * np.sin(phi)) / (2 * r)
    B[3, 5] = -((1 + n) * np.sin(phi)) / (2 * r)

    B[4, 1] = -1 / L
    B[4, 2] = -(1 - n) / 2
    B[4, 4] = 1 / L
    B[4, 5] = -(1 + n) / 2

    # Passage au repère global
    P1 = np.array([[np.cos(phi), -np.sin(phi), 0],
                   [np.sin(phi), np.cos(phi), 0],
                   [0, 0, -1]])
    P2 = np.zeros((3, 3))
    Passage = np.block([[P1, P2], [P2, P1]])
    P = np.linalg.inv(Passage)

    # Définition des composantes de la matrice D
    D = np.zeros((5, 5))

    D[0, 0] = E1 * t / (1 - v13 * v31)
    D[0, 1] = v13 * E3 * t / (1 - v13 * v31)

    D[1, 0] = v31 * E1 * t / (1 - v13 * v31)
    D[1, 1] = E3 * t / (1 - v13 * v31)

    D[2, 2] = E1 * t**3 / (12 * (1 - v13 * v31))
    D[2, 3] = v13 * E3 * t**3 / (12 * (1 - v13 * v31))

    D[3, 2] = v31 * E1 * t**3 / (12 * (1 - v13 * v31))
    D[3, 3] = E3 * t**3 / (12 * (1 - v13 * v31))

    D[4, 4] = (E1 * coef * t) / (2 * (1 + v13))  # E1 et v13 a corig si composite

    # Matrice de raideur d'un élement coque
    K = np.dot(np.dot(B.T, D), B) * r * L * np.pi
    Kp = np.dot(np.dot(Passage, K), P)

    return Kp
