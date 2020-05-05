import VBBinaryLensing
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# 2D magnification

VBBL = VBBinaryLensing.VBBinaryLensing()

scale = 0.8
index = 512

s = 1.2
q = 0.005

rho = 0.01

x = np.linspace(-scale,scale,index)
y = np.linspace(-scale,scale,index)

X,Y = np.meshgrid(x,y)

def mag(x,y):
    return VBBL.BinaryMag2(s,q,x,y,rho)

result = []

for j in range(index):
    temp = []
    for i in range(index):
        temp.append(mag(X[j][i],Y[j][i]))
    result.append(temp)

plt.figure()
plt.pcolormesh(X,Y,result)
plt.axis('scaled')
plt.colorbar()
plt.show()
