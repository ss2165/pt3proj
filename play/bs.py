from scipy.integrate import odeint
import numpy as np

import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

bloch = np.array([0, 0, 1])
TF = 5*2*np.pi
#TF = 200
NT = 2000
d = 1
Gam = 0.1

rabi0 = np.array([1, 0, d])

damped = True


def rabi(t):
    return [1, 0, -10+t*(20/TF)]

def f(a, t):
    der = np.cross(a, rabi(t))

    if damped:
        der += [-Gam*a[0]/2., -Gam*a[1]/2., -Gam*(a[2]-1)]
    return der


t = np.linspace(0, TF, NT)
y = odeint(f, bloch, t)

L1 = y.reshape(-1)
L2 = L1*L1
L3 = L2.reshape((len(L2)/3, 3))
L = np.sqrt(np.sum(L3, axis=1))


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')

X = y[:, 0]
Y = y[:, 1]
Z = y[:, 2]

# Create cubic bounding box to simulate equal aspect ratio
max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
# Comment or uncomment following both lines to test the fake bounding box:
for xb, yb, zb in zip(Xb, Yb, Zb):
   ax.plot([xb], [yb], [zb], 'w')

plt.grid()

ax.set_xlim(-1, 1)
ax.set_ylim(bottom=-1, top=1)
ax.set_zlim(bottom=-1, top=1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title("Bloch vector oscillations")
ax.plot(X, Y, zs=Z)



# plt.plot(t, L)
plt.show()
