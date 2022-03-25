import swiftsimio
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
import numpy as np
import copy
import unyt
import sys
import os
import csv

# Plot parameters
params = {
'axes.labelsize': 14,
'axes.titlesize': 14,
'font.size': 14,
'legend.fontsize': 14,
'xtick.labelsize': 12,
'ytick.labelsize': 12,
'xtick.direction': 'in',
'ytick.direction': 'in',
'xtick.top': True,
'ytick.right': True,
'xtick.major.width': 1.5,
'ytick.major.width': 1.5,
'axes.linewidth': 1.5,
'text.usetex': True,
'figure.figsize' : (10,4),
#'figure.figsize' : (9.90,3.25),
'figure.subplot.left'    : 0.045,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.05,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 1,
'lines.linewidth' : 2.,
#'text.latex.unicode': True
}
mpl.rcParams.update(params)
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Times']})

scatterplot_kwargs = {
    "alpha": 0.6,
    "s": 4,
    "marker": ".",
    "linewidth": 0.0,
    "facecolor": "blue",
}

# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True

snapshot_base = "output"



def get_TT1Dsolution():
    TT1D_runit = 5.4 * unyt.unyt_array(1.0,'kpc') #kpc
    xtt1dlist=np.array([])
    rtt1dlist=np.array([])
    with open('data/xTT1D_Stromgren100Myr.csv') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            rtt1dlist=np.append(rtt1dlist,float(row[0]))
            xtt1dlist=np.append(xtt1dlist,10**float(row[1]))
    rtt1dlist *= TT1D_runit

    Ttt1dlist=np.array([])
    rTtt1dlist=np.array([])
    with open('data/TTT1D_Stromgren100Myr.csv') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            rTtt1dlist=np.append(rTtt1dlist,float(row[0]))
            Ttt1dlist=np.append(Ttt1dlist,10**float(row[1]))
    rTtt1dlist *= TT1D_runit
    Ttt1dlist *= unyt.unyt_array(1.0, 'K')
    outdict = {}
    outdict['rtt1dlist'] = rtt1dlist
    outdict['xtt1dlist'] = xtt1dlist
    outdict['rTtt1dlist'] = rTtt1dlist
    outdict['Ttt1dlist'] = Ttt1dlist    
    return outdict

def mean_molecular_weight(XH0, XHp, XHe0, XHep, XHepp):
    """
    Determines the mean molecular weight for given 
    mass fractions of
        hydrogen:   XH0
        H+:         XHp
        He:         XHe0
        He+:        XHep
        He++:       XHepp

    returns:
        mu: mean molecular weight [in atomic mass units]
        NOTE: to get the actual mean mass, you still need
        to multiply it by m_u, as is tradition in the formulae
    """

    # 1/mu = sum_j X_j / A_j * (1 + E_j)
    # A_H    = 1, E_H    = 0
    # A_Hp   = 1, E_Hp   = 1
    # A_He   = 4, E_He   = 0
    # A_Hep  = 4, E_Hep  = 1
    # A_Hepp = 4, E_Hepp = 2
    one_over_mu = XH0 + 2 * XHp + 0.25 * XHe0 + 0.5 * XHep + 0.75 * XHepp

    return 1.0 / one_over_mu


def gas_temperature(u, mu, gamma):
    """
    Compute the gas temperature given the specific internal 
    energy u and the mean molecular weight mu
    """

    # Using u = 1 / (gamma - 1) * p / rho
    #   and p = N/V * kT = rho / (mu * m_u) * kT

    T = u * (gamma - 1) * mu * unyt.atomic_mass_unit / unyt.boltzmann_constant

    return T.to("K")


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list
    """

    snaplist = []

    if plot_all:
        dirlist = os.listdir()
        for f in dirlist:
            if f.startswith(snapshot_basename) and f.endswith("hdf5"):
                snaplist.append(f)

        snaplist = sorted(snaplist)

    else:
        fname = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"
        if not os.path.exists(fname):
            print("Didn't find file", fname)
            quit(1)
        snaplist.append(fname)

    return snaplist

def get_imf(scheme,data):
    """
    Get the ion mass fraction (imf) according to the scheme.
    """
    if scheme.startswith("GEAR M1closure"):
        imf = data.gas.ion_mass_fractions
    elif scheme.startswith("SPH M1closure"):
        # atomic mass
        mamu = {'e':0.0, 'HI':1.0, 'HII':1.0, 'HeI':4.0, 'HeII':4.0, 'HeIII':4.0}
        mass_function_hydrogen = data.gas.rt_element_mass_fractions.hydrogen
        imf = copy.deepcopy(data.gas.rt_species_abundances)
        named_columns = data.gas.rt_species_abundances.named_columns
        for column in named_columns:
            # abundance is in n_X/n_H unit. We convert it to mass fraction by multipling mass fraction of H
            mass_function = getattr(data.gas.rt_species_abundances,column) * mass_function_hydrogen * mamu[column]
            setattr(imf,column, mass_function)
    return imf


def trim_paramstr(paramstr):
    # clean string up
    if paramstr.startswith("["):
        paramstr = paramstr[1:]
    if paramstr.endswith("]"):
        paramstr = paramstr[:-1]

    # transform string values to floats with unyts
    params = paramstr.split(",")
    paramtrimmed = []
    for er in params:
        paramtrimmed.append(float(er))
    return paramtrimmed


def plot_compare(filename):
    # Read in data first
    print("working on", filename)
    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))
    gamma = meta.hydro_scheme["Adiabatic index"][0]

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))

    imf = get_imf(scheme,data)
    xHI = imf.HI/(imf.HI+imf.HII)

    mu = mean_molecular_weight(imf.HI, imf.HII, imf.HeI, imf.HeII, imf.HeIII)
    data.gas.T = gas_temperature(data.gas.internal_energies, mu, gamma)


    outdict = get_TT1Dsolution()

    fig, ax = plt.subplots(1,2)

    ax[0].scatter(r,xHI, **scatterplot_kwargs)
    ax[0].plot(outdict['rtt1dlist'], outdict['xtt1dlist'], color='k',lw=2.0,label="TT1D")
    ax[0].set_ylabel('Neutral Fraction')
    xlabel_units_str = meta.boxsize.units.latex_representation()
    ax[0].set_xlabel("r [$" + xlabel_units_str + "$]") 
    ax[0].set_yscale('log')
    ax[0].set_xlim([0,10])

    ax[1].scatter(r,data.gas.T, **scatterplot_kwargs)
    ax[1].plot(outdict['rTtt1dlist'], outdict['Ttt1dlist'], color='k',lw=2.0,label="TT1D")
    ax[1].set_ylabel('T [K]')
    xlabel_units_str = meta.boxsize.units.latex_representation()
    ax[1].set_xlabel("r [$" + xlabel_units_str + "$]") 
    ax[1].set_yscale('log')
    ax[1].set_xlim([0,10])

    plt.tight_layout()
    figname = filename[:-5]
    figname += "-Stromgren3DMF.png"
    plt.savefig(figname)
    plt.close()


if __name__ == "__main__":
    snaplist = get_snapshot_list(snapshot_base)
    for f in snaplist:
        plot_compare(f)
