import numpy as np
from numerical_quadrature import integral_gauss, integral_gauss_v
# Construction de la géométrie et maillage
from scipy.sparse import lil_matrix


# Données
Fmax = 480664
E1 = 210000
E2 = 210000
v12 = 0.33
v21 = v12 * E2 / E1
r_int = 1020
r_ext = 1064
rr = 4  # rayon de raccordement
L = 40
h = 296  # 300
t = 5
N = 20
amp = 100
angle_virole = 90 * np.pi / 180

taille_elt_pied_bride = 2
taille_elt_virole = 2
nbr_element_raccr = 6
angle_sect_maillage_raccr = (90 / nbr_element_raccr) * np.pi / 180
nbr_elem_pied_bride = round(L / taille_elt_pied_bride)
nbr_elem_virole = round(h / taille_elt_virole)

import numpy as np
import matplotlib.pyplot as plt

# Maillage
noeud = []
elem = []
w = 0
noeud.append({'r': r_int, 'z': 300})

# Maillage virole
for k in range(2, nbr_elem_virole + 2):
    r = r_int
    z = noeud[0]['z'] - np.sin(angle_virole) * (h / nbr_elem_virole * (k - 1))
    noeud.append({'r': r, 'z': z})
    elem.append({
        'noeud1': k - 1,
        'noeud2': k,
        'E1': E1,
        'E2': E2,
        'nu12': v12,
        'type': 1,  # coque cylindrique
        'epaiss': t,
        'longueur': h / nbr_elem_virole,
        'rd': noeud[k-1]['r'],
        'rg': noeud[k-2]['r'],
        'angle': 0
    })

angle = [168.749161, 146.250778, 123.750511, 101.249802]

# Maillage raccord
for j in range(k + 1, nbr_element_raccr + nbr_elem_virole + 2):
    r = noeud[k-1]['r'] + (rr - np.cos(angle_sect_maillage_raccr * (j - k)) * rr)
    z = noeud[k-1]['z'] - np.sin(angle_sect_maillage_raccr * (j - k)) * rr
    noeud.append({'r': r, 'z': z})
    elem.append({
        'noeud1': j - 1,
        'noeud2': j,
        'E1': E1,
        'E2': E2,
        'nu12': v12,
        'type': 2,
        'epaiss': t,
        'angle': np.arctan((noeud[j-1]['r'] - noeud[j - 2]['r']) / np.abs(noeud[j-1]['z'] - noeud[j - 2]['z'])),
        'longueur': np.sqrt((noeud[j-1]['r'] - noeud[j - 2]['r'])**2 + (noeud[j-1]['z'] - noeud[j - 2]['z'])**2),
        'rd': noeud[j-1]['r'],
        'rg': noeud[j - 2]['r']
    })

# Maillage pied de bride
for i in range(j + 1, nbr_elem_pied_bride + nbr_element_raccr + nbr_elem_virole + 2):
    r = noeud[j-1]['r'] + (i - j) * (L / nbr_elem_pied_bride)
    z = 0
    noeud.append({'r': r, 'z': z})
    elem.append({
        'noeud1': i - 1,
        'noeud2': i,
        'E1': E1,
        'E2': E2,
        'nu12': v12,
        'type': 3,
        'epaiss': t,
        'longueur': L / nbr_elem_pied_bride,
        'angle': 90,
        'rg': noeud[i-2]['r'],
        'rd': noeud[i-1]['r']
    })

xsm1 = [n['r'] for n in noeud]
ysm1 = [n['z'] for n in noeud]

# # Visualization
# plt.figure()
# plt.scatter(xsm1, ysm1, color='r')
# plt.axis([r_int - 100, r_ext + 100, -100, h + 100])
# plt.axis('equal')

# plt.figure()
# plt.plot(xsm1, ysm1, 'r-o')
# plt.axis([r_int - 100, r_ext + 100, -100, h + 100])
# plt.axis('equal')
# plt.show()


Nombre_noeuds = len(noeud)
ddl = 3
taille = Nombre_noeuds * ddl
Global_stiffness = lil_matrix((taille, taille))

for i in range(len(elem)):
    if elem[i]['type'] == 3:  # pied de bride
        E1 = elem[i]['E1']
        E3 = elem[i]['E2']
        v13 = elem[i]['nu12']
        v31 = elem[i]['nu12'] * E3 / E1
        t = elem[i]['epaiss']
        L = taille_elt_pied_bride
        phi = elem[i]['angle'] * np.pi / 180
        rm = elem[i]['rg'] + (elem[i]['rd'] - elem[i]['rg']) / 2
        r2 = (elem[i]['rd'] - elem[i]['rg']) / 2

        k_elem = integral_gauss(E1, E3, v13, v31, rm, t, L, phi, r2)

    elif elem[i]['type'] == 2:  # raccord
        E1 = elem[i]['E1']
        E3 = elem[i]['E2']
        v13 = elem[i]['nu12']
        v31 = elem[i]['nu12'] * E3 / E1
        rm = elem[i]['rg'] + (elem[i]['rd'] - elem[i]['rg']) / 2
        r2 = (elem[i]['rd'] - elem[i]['rg']) / 2
        t = elem[i]['epaiss']
        L = elem[i]['longueur']
        phi = elem[i]['angle']

        k_elem = integral_gauss(E1, E3, v13, v31, rm, t, L, phi, r2)

    elif elem[i]['type'] == 1:  # virole
        E1 = elem[i]['E1']
        E3 = elem[i]['E2']
        v13 = elem[i]['nu12']
        v31 = elem[i]['nu12'] * E3 / E1
        t = elem[i]['epaiss']
        L = elem[i]['longueur']
        phi = elem[i]['angle'] * np.pi / 180
        rm = elem[i]['rg']

        k_elem = integral_gauss_v(E1, E3, v13, v31, rm, t, L, phi)

    for l in range(1, 7):
        if l <= 3:
            a = 3 * (elem[i]['noeud1'] - 1) + l
        else:
            a = 3 * (elem[i]['noeud2'] - 1) + (l - 3)
        for m in range(1, 7):
            if m <= 3:
                b = 3 * (elem[i]['noeud1'] - 1) + m
            else:
                b = 3 * (elem[i]['noeud2'] - 1) + (m - 3)
            Global_stiffness[a-1, b-1] += k_elem[l-1, m-1]

# Removing the first node from assembly (CL: blocking node 1)
K = np.delete(Global_stiffness.toarray(), [0, 1, 2], axis=0)
K = np.delete(K, [0, 1, 2], axis=1)

K = lil_matrix(K)

# Creation du vecteur effort
Effort = np.zeros((Global_stiffness.shape[0] - 3, 1))

print("Stiffness construction part to be revised")
