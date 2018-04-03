#!/usr/bin/env python3

from pdspy.constants.physics import c, m_p, G
from pdspy.constants.physics import k as k_b
from pdspy.constants.astronomy import M_sun, AU
from matplotlib.backends.backend_pdf import PdfPages
import pdspy.modeling.mpi_pool
import pdspy.interferometry as uv
import pdspy.spectroscopy as sp
import pdspy.modeling as modeling
import pdspy.imaging as im
import pdspy.misc as misc
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects
from matplotlib.colors import LinearSegmentedColormap
import pdspy.table
import pdspy.dust as dust
import pdspy.gas as gas
import pdspy.mcmc as mc
import scipy.signal
import argparse
import numpy
import time
import sys
import os
import emcee
import corner
from mpi4py import MPI

comm = MPI.COMM_WORLD

################################################################################
#
# Parse command line arguments.
#
################################################################################

# Set the number of cpus to use.

ncpus = 4

# Get the source name and check that it has been set correctly.

source = "NSFSPG"

################################################################################
#
# Create a function which returns a model of the data.
#
################################################################################

def model(visibilities, params, parameters, plot=False):

    # Set the values of all of the parameters.

    p = {}
    for key in parameters:
        if parameters[key]["fixed"]:
            if parameters[key]["value"] in parameters.keys():
                if parameters[parameters[key]["value"]]["fixed"]:
                    value = parameters[parameters[key]["value"]]["value"]
                else:
                    value = params[parameters[key]["value"]]
            else:
                value = parameters[key]["value"]
        else:
            value = params[key]

        if key[0:3] == "log":
            p[key[3:]] = 10.**value
        else:
            p[key] = value

    # Make sure alpha and beta are defined.

    if p["disk_type"] == "exptaper":
        t_rdisk = p["T0"] * (p["R_disk"] / 1.)**-p["q"]
        p["h_0"] = ((k_b*(p["R_disk"]*AU)**3*t_rdisk) / (G*p["M_star"]*M_sun * \
                2.37*m_p))**0.5 / AU
    else:
        p["h_0"] = ((k_b * AU**3 * p["T0"]) / (G*p["M_star"]*M_sun * \
                2.37*m_p))**0.5 / AU
    p["beta"] = 0.5 * (3 - p["q"])
    p["alpha"] = p["gamma"] + p["beta"]

    # Set up the dust.

    dustopac = "pollack_new.hdf5"

    dust_gen = dust.DustGenerator(dust.__path__[0]+"/data/"+dustopac)

    ddust = dust_gen(p["a_max"] / 1e4, p["p"])
    edust = dust_gen(1.0e-4, 3.5)

    # Set up the gas.

    gases = []
    abundance = []

    index = 1
    while index > 0:
        if "gas_file"+str(index) in p:
            g = gas.Gas()
            g.set_properties_from_lambda(gas.__path__[0]+"/data/"+\
                    p["gas_file"+str(index)])

            gases.append(g)
            abundance.append(p["abundance"+str(index)])

            index += 1
        else:
            index = -1

    # Make sure we are in a temp directory to not overwrite anything.

    original_dir = os.environ["PWD"]
    os.mkdir("/tmp/temp_{1:s}_{0:d}".format(comm.Get_rank(), source))
    os.chdir("/tmp/temp_{1:s}_{0:d}".format(comm.Get_rank(), source))

    # Write the parameters to a text file so it is easy to keep track of them.

    f = open("params.txt","w")
    for key in p:
        f.write("{0:s} = {1}\n".format(key, p[key]))
    f.close()

    # Set up the model. 

    m = modeling.YSOModel()
    m.add_star(mass=p["M_star"], luminosity=p["L_star"],temperature=p["T_star"])

    if p["envelope_type"] == "ulrich":
        m.set_spherical_grid(p["R_in"], p["R_env"], 100, 51, 2, code="radmc3d")
    else:
        m.set_spherical_grid(p["R_in"], max(5*p["R_disk"],300), 100, 51, 2, \
                code="radmc3d")

    if p["disk_type"] == "exptaper":
        m.add_pringle_disk(mass=p["M_disk"], rmin=p["R_in"], rmax=p["R_disk"], \
                plrho=p["alpha"], h0=p["h_0"], plh=p["beta"], dust=ddust, \
                t0=p["T0"], plt=p["q"], gas=gases, abundance=abundance,\
                aturb=p["a_turb"])
    else:
        m.add_disk(mass=p["M_disk"], rmin=p["R_in"], rmax=p["R_disk"], \
                plrho=p["alpha"], h0=p["h_0"], plh=p["beta"], dust=ddust, \
                t0=p["T0"], plt=p["q"], gas=gases, abundance=abundance,\
                aturb=p["a_turb"])

    if p["envelope_type"] == "ulrich":
        m.add_ulrich_envelope(mass=p["M_env"], rmin=p["R_in"], rmax=p["R_env"],\
                cavpl=p["ksi"], cavrfact=p["f_cav"], dust=edust, \
                t0=p["T0_env"], tpl=p["q_env"], gas=gases, abundance=abundance,\
                aturb=p["a_turb_env"])
    else:
        pass

    m.grid.set_wavelength_grid(0.1,1.0e5,500,log=True)

    # Run the images/visibilities/SEDs.

    for j in range(len(visibilities["file"])):
        # Shift the wavelengths by the velocities.

        b = p["v_sys"]*1.0e5 / c
        lam = c / visibilities["data"][j].freq / 1.0e-4
        wave = lam * numpy.sqrt((1. - b) / (1. + b))

        # Set the wavelengths for RADMC3D to use.

        m.set_camera_wavelength(wave)

        if p["docontsub"]:
            m.run_image(name=visibilities["lam"][j], nphot=1e5, \
                    npix=visibilities["npix"][j], lam=None, \
                    pixelsize=visibilities["pixelsize"][j], tgas_eq_tdust=True,\
                    scattering_mode_max=0, incl_dust=True, incl_lines=True, \
                    loadlambda=True, incl=p["i"], pa=p["pa"], dpc=p["dpc"], \
                    code="radmc3d", verbose=False, writeimage_unformatted=True,\
                    setthreads=ncpus)

            m.run_image(name="cont", nphot=1e5, \
                    npix=visibilities["npix"][j], lam=None, \
                    pixelsize=visibilities["pixelsize"][j], tgas_eq_tdust=True,\
                    scattering_mode_max=0, incl_dust=True, incl_lines=False, \
                    loadlambda=True, incl=p["i"], pa=p["pa"], dpc=p["dpc"], \
                    code="radmc3d", verbose=False, writeimage_unformatted=True,\
                    setthreads=ncpus)

            m.images[visibilities["lam"][j]].image -= m.images["cont"].image
        else:
            m.run_image(name=visibilities["lam"][j], nphot=1e5, \
                    npix=visibilities["npix"][j], lam=None, \
                    pixelsize=visibilities["pixelsize"][j], tgas_eq_tdust=True,\
                    scattering_mode_max=0, incl_dust=False, incl_lines=True, \
                    loadlambda=True, incl=p["i"], pa=p["pa"], dpc=p["dpc"], \
                    code="radmc3d", verbose=False, writeimage_unformatted=True,\
                    setthreads=ncpus)

        m.visibilities[visibilities["lam"][j]] = uv.interpolate_model(\
                visibilities["data"][j].u, visibilities["data"][j].v, \
                visibilities["data"][j].freq, \
                m.images[visibilities["lam"][j]], dRA=-p["x0"], dDec=-p["y0"], \
                nthreads=ncpus)

        if plot:
            lam = c / visibilities["image"][j].freq / 1.0e-4
            wave = lam * numpy.sqrt((1. - b) / (1. + b))

            m.set_camera_wavelength(wave)

            if p["docontsub"]:
                m.run_image(name=visibilities["lam"][j], nphot=1e5, \
                        npix=visibilities["image_npix"][j], lam=None, \
                        pixelsize=visibilities["image_pixelsize"][j], \
                        tgas_eq_tdust=True, scattering_mode_max=0, \
                        incl_dust=True, incl_lines=True, loadlambda=True, \
                        incl=p["i"], pa=-p["pa"], dpc=p["dpc"], code="radmc3d",\
                        verbose=False, setthreads=ncpus)

                m.run_image(name="cont", nphot=1e5, \
                        npix=visibilities["image_npix"][j], lam=None, \
                        pixelsize=visibilities["image_pixelsize"][j], \
                        tgas_eq_tdust=True, scattering_mode_max=0, \
                        incl_dust=True, incl_lines=False, loadlambda=True, \
                        incl=p["i"], pa=-p["pa"], dpc=p["dpc"], code="radmc3d",\
                        verbose=False, setthreads=ncpus)

                m.images[visibilities["lam"][j]].image -= m.images["cont"].image
            else:
                m.run_image(name=visibilities["lam"][j], nphot=1e5, \
                        npix=visibilities["image_npix"][j], lam=None, \
                        pixelsize=visibilities["image_pixelsize"][j], \
                        tgas_eq_tdust=True, scattering_mode_max=0, \
                        incl_dust=False, incl_lines=True, loadlambda=True, \
                        incl=p["i"], pa=-p["pa"], dpc=p["dpc"], code="radmc3d",\
                        verbose=False, setthreads=ncpus)

            x, y = numpy.meshgrid(numpy.linspace(-256,255,512), \
                    numpy.linspace(-256,255,512))

            beam = misc.gaussian2d(x, y, 0., 0., \
                    visibilities["image"][j].header["BMAJ"]/2.355/\
                    visibilities["image"][j].header["CDELT2"], \
                    visibilities["image"][j].header["BMIN"]/2.355/\
                    visibilities["image"][j].header["CDELT2"], \
                    (90-visibilities["image"][j].header["BPA"])*\
                    numpy.pi/180., 1.0)

            for ind in range(len(wave)):
                m.images[visibilities["lam"][j]].image[:,:,ind,0] = \
                        scipy.signal.fftconvolve(\
                        m.images[visibilities["lam"][j]].image[:,:,ind,0], \
                        beam, mode="same")

    os.system("rm params.txt")
    os.chdir(original_dir)
    os.system("rmdir /tmp/temp_{1:s}_{0:d}".format(comm.Get_rank(), source))

    return m

# Define a useful class for plotting.

class Transform:
    def __init__(self, xmin, xmax, dx, fmt):
        self.xmin = xmin
        self.xmax = xmax
        self.dx = dx
        self.fmt = fmt

    def __call__(self, x, p):
        return self.fmt% ((x-(self.xmax-self.xmin+1)/2)*self.dx)

################################################################################
#
# In case we are restarting this from the same job submission, delete any
# temporary directories associated with this run.
#
################################################################################

os.system("rm -r /tmp/temp_{0:s}_{1:d}".format(source, comm.Get_rank()))

################################################################################
#
# Read in the data.
#
################################################################################

# Import the configuration file information.

sys.path.insert(0, '')

from config import *

# Set up the places where we will put all of the data.

visibilities["data"] = []
visibilities["data1d"] = []
visibilities["image"] = []

images["data"] = []

# Make sure "fmt" is in the visibilities dictionary.

if not "fmt" in visibilities:
    visibilities["fmt"] = ['4.1f' for i in range(len(visibilities["file"]))]

# Decide whether to use an exponentially tapered 

if not "disk_type" in parameters:
    parameters["disk_type"] = {"fixed":True, "value":"truncated", \
            "limits":[0.,0.]}

# Make sure the code doesn't break if envelope_type isn't specified.

if not "envelope_type" in parameters:
    parameters["envelope_type"] = {"fixed":True, "value":"none", \
            "limits":[0.,0.]}

# Decide whether to do continuum subtraction or not.

if not "docontsub" in parameters:
    parameters["docontsub"] = {"fixed":True, "value":False, "limits":[0.,0.]}

# Make sure the code is backwards compatible to a time when only a single gas
# file was being supplied.

if "gas_file" in parameters:
    parameters["gas_file1"] = parameters["gas_file"]
    parameters["logabundance1"] = parameters["logabundance"]

######################################
# Read in the millimeter visibilities.
######################################

for j in range(len(visibilities["file"])):
    # Read the raw data.

    data = uv.Visibilities()
    data.read(visibilities["file"][j])

    # Center the data. => need to update!

    data = uv.center(data, [visibilities["x0"][j], visibilities["y0"][j], 1.])

    # Add the data to the dictionary structure.

    visibilities["data"].append(data)

    # Scale the weights of the visibilities to force them to be fit well.

    visibilities["data"][j].weights *= visibilities["weight"][j]

    # Average the visibilities radially.

    visibilities["data1d"].append(uv.average(data, gridsize=20, radial=True, \
            log=True, logmin=data.uvdist[numpy.nonzero(data.uvdist)].min()*\
            0.95, logmax=data.uvdist.max()*1.05, mode="spectralline"))

    # Read in the image.

    visibilities["image"].append(im.readimfits(visibilities["image_file"][j]))

#####################
# Read in the images.
#####################

for j in range(len(images["file"])):
    images["data"].append(im.readimfits(images["file"][j]))

################################################################################
#
# Fit the model to the data.
#
################################################################################

# Set up the emcee run.

ndim = 0
for key in parameters:
    if not parameters[key]["fixed"]:
        ndim += 1

chain = numpy.load("chain.npy")

# Get the best fit parameters and uncertainties.

samples = chain[:,-100:,:].reshape((-1, ndim))

params = numpy.median(samples, axis=0)
sigma = samples.std(axis=0)

# Make a dictionary of the best fit parameters.

keys = []
for key in sorted(parameters.keys()):
    if not parameters[key]["fixed"]:
        keys.append(key)

params = dict(zip(keys, params))

############################################################################
#
# Plot the results.
#
############################################################################

# Create a high resolution model for averaging.

m = model(visibilities, params, parameters, plot=True)

# Loop through the visibilities and plot.

for j in range(len(visibilities["file"])):
    # Calculate the velocity for each image.

    v = c * (float(visibilities["freq"][j])*1.0e9 - \
            visibilities["image"][j].freq)/(float(visibilities["freq"][j])*\
            1.0e9)

    # Set the ticks.

    ticks = visibilities["image_ticks"][j]

    # Plot the best fit model over the data.

    fig, ax = plt.subplots(nrows=2, ncols=7, sharex=True, sharey=True)

    for n in range(2):
        if n == 0:
            plot_image = visibilities["image"][j]

            cont_image = images["data"][j]

            # The plot range.

            vmin = 0
            vmax = numpy.nanmax(visibilities["image"][j].image)

            levels = numpy.linspace(3, 30, 10)*numpy.nanstd(plot_image.image)

            # Get the centroid position.

            if "x0" in params:
                xmin, xmax = int(round(visibilities["image_npix"][j]/2 + \
                        visibilities["x0"][j]/\
                        visibilities["image_pixelsize"][j]+ \
                        params["x0"]/visibilities["image_pixelsize"][j]+ \
                        ticks[0]/visibilities["image_pixelsize"][j])), \
                        int(round(visibilities["image_npix"][j]/2+\
                        visibilities["x0"][j]/\
                        visibilities["image_pixelsize"][j]+ \
                        params["x0"]/visibilities["image_pixelsize"][j]+ \
                        ticks[-1]/visibilities["image_pixelsize"][j]))

                cont_xmin, cont_xmax = int(round(images["npix"][j]/2 + \
                        images["x0"][j]/images["pixelsize"][j]+ \
                        ticks[0]/images["pixelsize"][j])), \
                        int(round(images["npix"][j]/2+\
                        images["x0"][j]/images["pixelsize"][j]+\
                        ticks[-1]/images["pixelsize"][j]))
            else:
                xmin, xmax = int(round(visibilities["image_npix"][j]/2 + \
                        visibilities["x0"][j]/\
                        visibilities["image_pixelsize"][j]+ \
                        parameters["x0"]["value"]/\
                        visibilities["image_pixelsize"][j]+ \
                        ticks[0]/visibilities["image_pixelsize"][j])), \
                        int(round(visibilities["image_npix"][j]/2+\
                        visibilities["x0"][j]/\
                        visibilities["image_pixelsize"][j]+ \
                        parameters["x0"]["value"]/\
                        visibilities["image_pixelsize"][j]+ \
                        ticks[-1]/visibilities["image_pixelsize"][j]))
            if "y0" in params:
                ymin, ymax = int(round(visibilities["image_npix"][j]/2-\
                        visibilities["y0"][j]/\
                        visibilities["image_pixelsize"][j]- \
                        params["y0"]/visibilities["image_pixelsize"][j]+ \
                        ticks[0]/visibilities["image_pixelsize"][j])), \
                        int(round(visibilities["image_npix"][j]/2-\
                        visibilities["y0"][j]/\
                        visibilities["image_pixelsize"][j]- \
                        params["y0"]/visibilities["image_pixelsize"][j]+ \
                        ticks[-1]/visibilities["image_pixelsize"][j]))

                cont_ymin, cont_ymax = int(round(images["npix"][j]/2 - \
                        images["y0"][j]/images["pixelsize"][j]+ \
                        ticks[0]/images["pixelsize"][j])), \
                        int(round(images["npix"][j]/2-\
                        images["y0"][j]/images["pixelsize"][j]+\
                        ticks[-1]/images["pixelsize"][j]))
            else:
                ymin, ymax = int(round(visibilities["image_npix"][j]/2-\
                        visibilities["y0"][j]/\
                        visibilities["image_pixelsize"][j]- \
                        parameters["y0"]["value"]/\
                        visibilities["image_pixelsize"][j]+ \
                        ticks[0]/visibilities["image_pixelsize"][j])), \
                        int(round(visibilities["image_npix"][j]/2-\
                        visibilities["y0"][j]/\
                        visibilities["image_pixelsize"][j]- \
                        parameters["y0"]["value"]/\
                        visibilities["image_pixelsize"][j]+ \
                        ticks[-1]/visibilities["image_pixelsize"][j]))
        elif n == 1:
            # Create a model image.

            m.visibilities[visibilities["lam"][j]].weights = \
                    visibilities["data"][j].weights

            model_image = uv.invert(m.visibilities[visibilities["lam"][j]], \
                    imsize=visibilities["image_npix"][j], \
                    pixel_size=visibilities["image_pixelsize"][j], \
                    weighting="natural", convolution="pillbox", \
                    mfs=False, mode="spectralline")

            plot_image = model_image
            
            # Get the scaling for the image.

            temp = uv.invert(visibilities["data"][j], \
                    imsize=visibilities["image_npix"][j], \
                    pixel_size=visibilities["image_pixelsize"][j], \
                    weighting="natural", convolution="pillbox", \
                    mfs=False, mode="spectralline")

            vmin = 0
            vmax = numpy.nanmax(temp.image)

            levels = numpy.linspace(3, 30, 10) * numpy.nanstd(temp.image)
            
            # Now make the centroid the map center for the model.

            xmin, xmax = int(visibilities["image_npix"][j]/2+1 + \
                    ticks[0]/visibilities["image_pixelsize"][j] + 0.5), \
                    int(visibilities["image_npix"][j]/2+1 +\
                    ticks[-1]/visibilities["image_pixelsize"][j] + 0.5)
            ymin, ymax = int(visibilities["image_npix"][j]/2+1 + \
                    ticks[0]/visibilities["image_pixelsize"][j] + 0.5), \
                    int(visibilities["image_npix"][j]/2+1 + \
                    ticks[-1]/visibilities["image_pixelsize"][j] + 0.5)
        elif n == 2:
            # Create a residual image. Centroid and scaling remain the same.

            residuals = uv.Visibilities(visibilities["data"][j].u, \
                    visibilities["data"][j].v, visibilities["data"][j].freq, \
                    visibilities["data"][j].real.copy(), \
                    visibilities["data"][j].imag.copy(),\
                    visibilities["data"][j].weights)

            residuals.real -= m.visibilities[visibilities["lam"][j]].real
            residuals.imag -= m.visibilities[visibilities["lam"][j]].imag

            residual_image = uv.invert(residuals, \
                    imsize=visibilities["image_npix"][j], \
                    pixel_size=visibilities["image_pixelsize"][j], \
                    weighting="natural", convolution="pillbox", \
                    mfs=False, mode="spectralline")

            plot_image = residual_image

        # Plot the image.

        for l in range(7):
            ind = 5*l + visibilities["ind0"][j]

            if v[ind]/1e5 < params["v_sys"]:
                cdict1 = {'red':   ((0.0, 1.0, 1.0),
                                    (1.0, 0.0, 0.0)),
                          'green': ((0.0, 1.0, 1.0),
                                    (1.0, 0.0, 0.0)),
                          'blue':  ((0.0, 1.0, 1.0),
                                    (1.0, 1.0, 1.0))}
                blues = LinearSegmentedColormap('blues', cdict1)
                cmap = blues
            else:
                cdict2 = {'red':   ((0.0, 1.0, 1.0),
                                    (1.0, 1.0, 1.0)),
                          'green': ((0.0, 1.0, 1.0),
                                    (1.0, 0.0, 0.0)),
                          'blue':  ((0.0, 1.0, 1.0),
                                    (1.0, 0.0, 0.0))}
                reds = LinearSegmentedColormap('reds', cdict2)
                cmap = reds

            plot_image.image[plot_image.image < 0] = 0.

            scalefn = lambda x: (numpy.arctan(x*7-12)+1.3) / \
                    (numpy.arctan(1*7-1.5)+1.3)

            # Plot the image.

            ax[n,l].imshow(plot_image.image[ymin:ymax,xmin:xmax,ind,0] * \
                    scalefn(abs(v[ind]/1e5 - params["v_sys"])), \
                    origin="lower", interpolation="nearest", vmin=vmin, \
                    vmax=vmax, cmap=cmap)

            ax[n,l].contour(plot_image.image[ymin:ymax,xmin:xmax,ind,0], \
                    levels=levels, colors='k')

            if n == 0:
                x = numpy.arange(cont_image.image[cont_ymin:cont_ymax,\
                        cont_xmin:cont_xmax,0,0].shape[1]) * \
                        images["pixelsize"][j] / \
                        visibilities["image_pixelsize"][j]
                y = numpy.arange(cont_image.image[cont_ymin:cont_ymax,\
                        cont_xmin:cont_xmax,0,0].shape[0]) * \
                        images["pixelsize"][j] / \
                        visibilities["image_pixelsize"][j]

                cont_levels = (numpy.arange(2)+0.5)/2.*\
                        numpy.nanmax(cont_image.image)

                ax[n,l].contour(x, y, cont_image.image[cont_ymin:cont_ymax,\
                        cont_xmin:cont_xmax,0,0], colors="gray", \
                        levels=cont_levels)

            # Fix the axes labels.

            transformx = ticker.FuncFormatter(Transform(xmin, xmax, \
                    visibilities["image_pixelsize"][j], '%.1f"'))
            transformy = ticker.FuncFormatter(Transform(ymin, ymax, \
                    visibilities["image_pixelsize"][j], '%.1f"'))

            ax[n,l].set_xticks(visibilities["image_npix"][j]/2+\
                    ticks[1:-1]/visibilities["image_pixelsize"][j]-xmin)
            ax[n,l].set_yticks(visibilities["image_npix"][j]/2+\
                    ticks[1:-1]/visibilities["image_pixelsize"][j]-ymin)

            ax[n,l].get_xaxis().set_major_formatter(transformx)
            ax[n,l].get_yaxis().set_major_formatter(transformy)

            ax[n,l].set_adjustable('box-forced')

            ax[-1,l].set_xlabel("$\Delta$RA")

            # If this is the data, show the beam and the velocities.

            if n == 0:
                # Add the velocity to the map.

                txt = ax[n,l].annotate(r"$v=%{0:s}$ km s$^{{-1}}$".format(\
                        visibilities["fmt"][j]) % (v[ind]/1e5),\
                        xy=(0.1,0.8), xycoords='axes fraction')

                # Show the size of the beam.

                bmaj = visibilities["image"][j].header["BMAJ"] / \
                        abs(visibilities["image"][j].header["CDELT1"])
                bmin = visibilities["image"][j].header["BMIN"] / \
                        abs(visibilities["image"][j].header["CDELT1"])
                bpa = visibilities["image"][j].header["BPA"]

                ax[n,l].add_artist(patches.Ellipse(xy=(12.5,12.5), \
                        width=bmaj, height=bmin, angle=(bpa+90), \
                        facecolor="white", edgecolor="black"))

        ax[n,0].set_ylabel("$\Delta$Dec")

    # Adjust the plot and save it.

    fig.set_size_inches((10.,3.125))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.98, bottom=0.15, \
            wspace=0.0,hspace=0.0)

    # Adjust the figure and save.

    fig.savefig("HOPS-370_{0:s}.pdf".format(visibilities["lam"][j]))

    plt.clf()
