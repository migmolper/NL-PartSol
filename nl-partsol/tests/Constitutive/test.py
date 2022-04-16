import numpy as np

Ndim = 2

d_phi_1D = np.array([10, 20, 60, 40])
d_phi_2D = np.reshape(d_phi_1D, (Ndim, Ndim))

b_e_1D = np.array([157, 671, 671, 1561])
b_e_2D = np.reshape(b_e_1D, (Ndim, Ndim))


b_e_new_1D = np.zeros((Ndim*Ndim,1))


for i in range(0, Ndim):
    for j in range(0, Ndim):
        for k in range(0, Ndim):
            for l in range(0, Ndim):
                b_e_new_1D[i * Ndim + j] = (
                    b_e_new_1D[i * Ndim + j]
                    + d_phi_1D[i * Ndim + k]
                    * b_e_1D[k * Ndim + l]
                    * d_phi_1D[j * Ndim + l]
                )


b_e_new_2D = np.einsum("ik,kl,lj", d_phi_2D, b_e_2D, d_phi_2D.transpose())

print(b_e_new_2D - np.reshape(b_e_new_1D, (Ndim, Ndim)))

