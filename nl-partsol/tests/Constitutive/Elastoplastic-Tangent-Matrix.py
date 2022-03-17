import numpy as np


dN_alpha = np.array([0.35, 0.42])
dN_beta = np.array([0.1, 0.8])
tau = np.array([[200, 345], [345, 652]])
D_phi = np.array([[1.2, 0.5], [0.7, 0.8]])
D_phi_mT = np.linalg.inv(D_phi).transpose()
b_e = np.array([[0.8, 0.6], [0.6, 0.1]])
a_ep = np.array([[10.11, 35.3], [35.3, 47.05]])

FPK = np.einsum("ij,jk", tau, D_phi_mT)

u = dN_beta
v = np.einsum("ij,j", D_phi_mT, dN_alpha)

eigval_b_e, eigvec_b_e = np.linalg.eigh(b_e)

eigval_T = np.linalg.eigvalsh(tau)

n_I = eigvec_b_e[:, 0]
n_II = eigvec_b_e[:, 1]

m_I = np.einsum("i,j", n_I, n_I)
m_I_II = np.einsum("i,j", n_I, n_II)
m_II_I = np.einsum("i,j", n_II, n_I)
m_II = np.einsum("i,j", n_II, n_II)

mIu = np.einsum("ij,j", m_I, u)
mI_IIu = np.einsum("ij,j", m_I_II, u)
mII_Iu = np.einsum("ij,j", m_II_I, u)
mIIu = np.einsum("ij,j", m_II, u)

mIv = np.einsum("ij,j", m_I, v)
mI_IIv = np.einsum("ij,j", m_I_II, v)
mII_Iv = np.einsum("ij,j", m_II_I, v)
mIIv = np.einsum("ij,j", m_II, v)

print("A_ep_AB:")
print(
    0.5
    * ((eigval_T[1] - eigval_T[0]) / (eigval_b_e[1] - eigval_b_e[0]))
    * (
        eigval_b_e[1] * np.einsum("i,j", mI_IIv, mI_IIu)
        + eigval_b_e[0] * np.einsum("i,j", mI_IIv, mII_Iu)
    )
)

print("A_ep_AB:")
print(
    0.5
    * ((eigval_T[0] - eigval_T[1]) / (eigval_b_e[0] - eigval_b_e[1]))
    * (
        eigval_b_e[0] * np.einsum("i,j", mII_Iv, mII_Iu)
        + eigval_b_e[1] * np.einsum("i,j", mII_Iv, mI_IIu)
    )
)


A_ep = (
    a_ep[0][0] * np.einsum("i,j", mIv, mIu)
    + a_ep[0][1] * np.einsum("i,j", mIv, mIIu)
    + a_ep[1][0] * np.einsum("i,j", mIIv, mIu)
    + a_ep[1][1] * np.einsum("i,j", mIIv, mIIu)
)

A_ep = (
    A_ep
    + 0.5
    * ((eigval_T[1] - eigval_T[0]) / (eigval_b_e[1] - eigval_b_e[0]))
    * (
        eigval_b_e[1] * np.einsum("i,j", mI_IIv, mI_IIu)
        + eigval_b_e[0] * np.einsum("i,j", mI_IIv, mII_Iu)
    )
    + 0.5
    * ((eigval_T[0] - eigval_T[1]) / (eigval_b_e[0] - eigval_b_e[1]))
    * (
        eigval_b_e[0] * np.einsum("i,j", mII_Iv, mII_Iu)
        + eigval_b_e[1] * np.einsum("i,j", mII_Iv, mI_IIu)
    )
)

print("A_ep:")
print(A_ep)
