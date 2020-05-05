import VBBinaryLensing
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

VBBL = VBBinaryLensing.VBBinaryLensing()
VBBL.RelTol = 1e-03

t = np.linspace(7120,7220,2000)

s = 1.2
q = 0.5e-2
rho = 0.00001
alpha = 0.93
tE = 100.3
t0 = 7154
u0 = -0.3

tau = (t-t0)/tE
y1 = -u0*np.sin(alpha)+tau*np.cos(alpha)
y2 = u0*np.cos(alpha)+tau*np.sin(alpha)

mag = np.zeros(len(tau))

params = [np.log(s), np.log(q), u0, alpha, np.log(rho), np.log(tE), t0]
mag = VBBL.BinaryLightCurve(params, t, y1, y2)

solutions = VBBL.PlotCrit(s, q) # Returns _sols object containing n crit. curves followed by n caustic curves
print(type(solutions))


def iterate_from(item):
    while item is not None:
        yield item
        item = item.next


curves = []
for curve in iterate_from(solutions.first):
    for point in iterate_from(curve.first):
        curves.append((point.x1, point.x2))

critical_curves = np.array(curves[:int(len(curves) / 2)])
caustic_curves = np.array(curves[int(len(curves) / 2):])


fig, ax = plt.subplots(figsize=(12,7))
ax.plot(t, mag, 'k-')
ax.grid(True)

ax2 = fig.add_axes([.54, .44, .4, .4], aspect=1)
ax2.plot(caustic_curves[:, 0], caustic_curves[:, 1], 'k-')
ax2.plot(y1, y2, 'k--')
ax2.grid(True)
ax2.set_xlim(-1,1)
ax2.set_ylim(-1,1)
plt.title('q = '+str(q)+',s='+str(s)+"(VBBL)")
plt.show()