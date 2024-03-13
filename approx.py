import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib import rcParams 
rcParams['font.family']='serif'
rcParams['xtick.top']=True
rcParams['xtick.direction']='in'
rcParams['ytick.right']=True
rcParams['ytick.direction']='in'
rcParams['ytick.major.size']=6
rcParams['ytick.minor.size']=3
rcParams['ytick.major.width']=1.
rcParams['ytick.minor.width']=1.
rcParams['ytick.minor.visible']=True
rcParams['xtick.major.size']=6
rcParams['xtick.minor.size']=3
rcParams['xtick.major.width']=1.
rcParams['xtick.minor.width']=1.
rcParams['xtick.minor.visible']=True
rcParams['lines.linewidth']=1.5


font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 19}
matplotlib.rc('font', **font)
matplotlib.rc('mathtext',**{'fontset' : 'cm'})

## IMPUT PARAMETERS ##
j = 0.1         # dimensionless angular momentum j = J / M^2
q_in = 3        # dimensionless quadrupole moment q = Q / M^3 * j^2
M = 1.4         # mass in solar mass units
b = 1.0         # beta to beta_cusp ratio - from 0 to 1
## END OF PARAMETERS ##


## CONSTANTS ##
c = 2.99792458e10  # speed of light in [cm/s]
G = 6.67428e-8  # gravitational constant in [cm^3/(g*s^2)]
Msun = 1.98892e33  # solar mass in [g]


R = c ** 3 / (G * M * Msun)  # unit conversion

nu0 = 1 / (2 * np.pi) * R / 6 ** (3 / 2)
q = q_in * np.square(j)

RADAXI = []
RADPR = []
VERTAXI = []
VERTPR = []
NK = []

nuk = 100  # Keplerian frequency in Hertz = variable

while nuk < 3000:

    r0 = np.power((np.power(6,1.5) - j * nuk / nu0 + 0.5 * (q - np.square(j)) * np.square(nuk / nu0)) / (nuk / nu0),(2 / 3))

    RR = 1 - np.square(b) * (1 / 33 + 10 / nu0 * np.square(np.power(r0, 1.5) - 18.4 + 15 * j - 4 * q) / (1 - 1.5 * j + 0.6 * q))
    BB = 1 - 0.2 * np.square(b) * (1 + j)
    KK = 1 - np.square(b) * (14 / nu0 * np.square(np.power(r0,1.5) - 40 / 3 + 10 * j - 2.5 * q) / (1 - 1.7 * j + 0.7 * q))
    TT = 1 - np.square(b) * (5 * j / nu0 * (1 - 2 / (1 - 2 * j - 2 * np.square(j) + 3 * q) * np.square(r0 - 6.1 + 2.5 * j)) + 9 * (
                q - np.square(j)) / nu0 * (r0 - 7.2))

    try:
        radaxisym = nuk * RR * np.sqrt(
            1 - 6 / r0 + 8 * j / np.power(r0,(3 / 2)) + 57 * np.square(j) / np.power(r0,(7 / 2)) - 3 * q / np.square(r0) * (
                        1 + 19 / np.power(r0,(3 / 2))))
    except:
        break
    radprec = nuk * (1 - BB * np.sqrt(
        1 - 6 / r0 + 8 * j / np.power(r0,(3 / 2)) + 57 * np.square(j) / np.power(r0, (7 / 2)) - 3 * q / np.square(r0) * (1 + 19 / np.power(r0, (3 / 2)))))
    vertaxisym = nuk * KK * np.sqrt(
        1 - 4 * j / np.power(r0, (3 / 2)) - 24 * np.square(j) / np.power(r0, (7 / 2)) + 3 * q / np.square(r0) * (1 + 8 / np.power(r0, (3 / 2))))
    vertprec = nuk * (1 - TT * np.sqrt(
        1 - 4 * j / np.power(r0,(3 / 2)) - 24 * np.square(j) / np.power(r0, (7 / 2)) + 3 * q / np.square(r0) * (1 + 8 / np.power(r0, (3 / 2)))))

    RADAXI.append(radaxisym)
    RADPR.append(radprec)
    VERTAXI.append(vertaxisym)
    VERTPR.append(vertprec)
    NK.append(nuk)

    nuk += 10

fig, (ax0,ax1,ax2,ax3) = plt.subplots(nrows = 1, ncols = 4, figsize = (24,8))
ax0.plot(NK, RADAXI)
ax0.set_ylim(0, 1000)
ax0.set_xlabel(r"$\nu_K\,\mathrm{(Hz)}$")
ax0.set_ylabel(r"$\nu_r\,\mathrm{(Hz)}$")
ax0.grid(alpha = 0.3)


ax1.plot(NK, RADPR)
ax1.set_ylim(0, 1000)
ax1.set_xlabel(r"$\nu_K\,\mathrm{(Hz)}$")
ax1.set_ylabel(r"$\nu_{PR}\,\mathrm{(Hz)}$")
ax1.grid(alpha = 0.3)


ax2.plot(NK, VERTAXI)
ax2.set_ylim(0, 3000)
ax2.set_xlabel(r"$\nu_K\,\mathrm{(Hz)}$")
ax2.set_ylabel(r"$\nu_V\,\mathrm{(Hz)}$")
ax2.grid(alpha = 0.3)


ax3.plot(NK, VERTPR)
ax3.set_ylim(0, 1000)
ax3.set_xlabel(r"$\nu_K\,\mathrm{(Hz)}$")
ax3.set_ylabel(r"$\nu_{LT}\,\mathrm{(Hz)}$")
ax3.grid(alpha = 0.3)


fig.suptitle(r'$j = {{{0:.2f}}}, q = {{{1:.2f}}} j^2, \beta = {{{2:.2f}}}, M = {{{3:.2f}}} M_\odot $'.format(j,q_in,b,M))
fig.tight_layout()
fig.savefig('freqs.pdf',bbox_inches='tight')
#fig.show()
