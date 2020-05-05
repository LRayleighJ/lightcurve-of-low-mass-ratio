import MulensModel as mm
import numpy as np
import matplotlib.pyplot as plt

t_left = 7120
t_right = 7220

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

my_1s2l_model = mm.Model({'t_0': t0,
                          'u_0': u0,
                          't_E': tE,
                          'rho': rho,
                          'q': q,
                          's': s,
                          'alpha': np.rad2deg(alpha)})

my_1s2l_model.set_magnification_methods([t_left, 'VBBL', t_right])

plt.figure()

my_1s2l_model.plot_magnification(t_range=[t_left, t_right], color='black')

plt.grid()
plt.title('q = '+str(q)+',s='+str(s)+"(mulensmodel)")
plt.show()