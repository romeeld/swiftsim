import swiftsimio
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.integrate import odeint
import numpy as np
import copy
import unyt
import sys
import os


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
'figure.figsize' : (5,4),
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

# analytic solution
def neutralfraction3d(rfunc,nH,sigma,alphaB,dNinj,rini):
    def fn(xn, rn):
        """this is the rhs of the ODE to integrate, i.e. dx/drn=fn(x,r)=x*(1-x)/(1+x)*(2/rn+x)"""
        return xn*(1.0-xn)/(1.0+xn)*(2.0/rn+xn)
    xn0 = nH * alphaB * 4.0 * np.pi / sigma / dNinj * rini * rini 
    rnounit = rfunc*nH*sigma
    xn = odeint(fn, xn0, rnounit)
    return xn 

def get_analytic_solution(data):
    meta = data.metadata
    rho = data.gas.densities
    rini_value = 0.1 
    r_ana = np.linspace(rini_value,10.0,100) * unyt.unyt_array(1.0, 'kpc')
    rini = rini_value * unyt.unyt_array(1.0, 'kpc')
    nH = np.mean(rho.to('g/cm**3')/unyt.proton_mass)
    sigma_cross = trim_paramstr(meta.parameters['SPHM1RT:sigma_cross'].decode("utf-8"))*unyt.unyt_array(1.0,'cm**2')
    sigma = sigma_cross[0]
    alphaB = trim_paramstr(meta.parameters['SPHM1RT:alphaB'].decode("utf-8"))*unyt.unyt_array(1.0,'cm**3/s')
    unit_l_in_cgs = float(meta.parameters["Snapshots:UnitLength_in_cgs"])*unyt.unyt_array(1.0,'cm')
    unit_v_in_cgs = float(meta.parameters["Snapshots:UnitVelocity_in_cgs"])*unyt.unyt_array(1.0,'cm/s')
    unit_m_in_cgs = float(meta.parameters["Snapshots:UnitMass_in_cgs"])*unyt.unyt_array(1.0,'g')
    star_emission_rates = trim_paramstr(meta.parameters['SPHM1RT:star_emission_rates'].decode("utf-8"))* unit_m_in_cgs * unit_v_in_cgs**3 / unit_l_in_cgs
    ionizing_photon_energy_erg = trim_paramstr(meta.parameters['SPHM1RT:ionizing_photon_energy_erg'].decode("utf-8"))*unyt.unyt_array(1.0,'erg')
    dNinj = star_emission_rates[1]/ionizing_photon_energy_erg[0]
    xn = neutralfraction3d(r_ana,nH,sigma,alphaB,dNinj,rini)
    return r_ana, xn


def plot_analytic_compare(filename):
    # Read in data first
    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))

    imf = get_imf(scheme,data)
    xHI = imf.HI/(imf.HI+imf.HII)

    r_ana,xn = get_analytic_solution(data)
    plt.scatter(r,xHI, **scatterplot_kwargs)
    plt.plot(r_ana,xn)
    plt.ylabel('Neutral Fraction')
    xlabel_units_str = meta.boxsize.units.latex_representation()
    plt.xlabel("r [$" + xlabel_units_str + "$]") 
    plt.yscale('log')
    plt.xlim([0,10])
    plt.tight_layout()
    figname = filename[:-5]
    figname += "-Stromgren3Dsinglebin.png"
    plt.savefig(figname)
    plt.close()


if __name__ == "__main__":
    snaplist = get_snapshot_list(snapshot_base)
    for f in snaplist:
        plot_analytic_compare(f)
