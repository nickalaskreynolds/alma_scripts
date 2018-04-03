#!/usr/bin/env python3

import pdspy.interferometry as uv
import pdspy.imaging as im
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy
import corner
import os

# Read in the data.

data = uv.Visibilities()
data.read("./L1448IRS3B.H13COp.ms.selfcal.contsub.concat.wide.lsrk.tavg.hdf5")

print(data.uvdist[data.uvdist > 0].min(),data.uvdist[data.uvdist > 0].max())

# Read in the image

image = im.readimfits("./L1448IRS3B.H13COp.triplet.fits")

# Average the data to a more manageable size.

vis = uv.average(data, gridsize=2048, binsize=11000)

# Average the visibilities radially.

data_1d = uv.average(data, gridsize=20, binsize=500000., radial=True, \
        log=True, logmin=data.uvdist[data.uvdist > 0].min()*0.95, \
        logmax=data.uvdist[data.uvdist > 0].max()*1.05)

# Fit the data with a gaussian model.

params, sigma, samples = uv.fit_model(vis, funct='gauss', nsteps=1e3)

# Write out the results.

f = open("gaussian_model/gaussian_fit.txt", "w")
f.write("Best fit to HOPS-370:\n\n")
f.write("x0 = {0:f} +/- {1:f}\n".format(params[0], sigma[0]))
f.write("y0 = {0:f} +/- {1:f}\n".format(params[1], sigma[1]))
f.write("r = {0:f} +/- {1:f}\n".format(params[2], sigma[2]))
f.write("i = {0:f} +/- {1:f}\n".format(params[3], sigma[3]))
f.write("P.A. = {0:f} +/- {1:f}\n".format(params[4], sigma[4]))
f.write("Flux = {0:f} +/- {1:f}\n\n".format(params[5], sigma[5]))
f.close()

print()
os.system("cat gaussian_model/gaussian_fit.txt")

# Plot histograms of the resulting parameters.

xlabels = ["x$_0$","y$_0$","$r$","$i$","P.A.",r"F$_{\nu}$"]

fig = corner.corner(samples, labels=xlabels, truths=params)

plt.savefig("gaussian_model/gaussian_fit.pdf")
plt.close(fig)

# Plot the best fit model over the data.

fig, ax = plt.subplots(nrows=1, ncols=2)

# Create a high resolution model for averaging.

u, v = numpy.meshgrid( numpy.linspace(-5000,5000,10000) * 1000., \
        numpy.linspace(-5000,5000,10000) * 1000.)
u = u.reshape((u.size,))
v = v.reshape((v.size,))

model = uv.model(u, v, params, return_type="data", funct="gauss")
model1d = uv.average(model, gridsize=10000, binsize=1000, radial=True)

# Plot the visibilities.

ax[0].errorbar(data_1d.uvdist/1000, data_1d.amp, \
        yerr=numpy.sqrt(1./data_1d.weights),\
        fmt="bo", markersize=8, markeredgecolor="b")

# Plot the best fit model

ax[0].plot(model1d.uvdist/1000, model1d.amp, "g-")

# Create a model image to contour over the image.

model = uv.model(data.u, data.v, params, return_type="data", funct="gauss")
model.weights = data.weights

model_image = uv.invert(model, imsize=1024, pixel_size=0.01, \
        centering=params)

# Plot the image.

xmin, xmax = 512+int(params[0]/0.01)-128, 512+int(params[0]/0.01)+128
ymin, ymax = 512-int(params[1]/0.01)-128, 512-int(params[1]/0.01)+128

ax[1].imshow(image.image[ymin:ymax,xmin:xmax,0,0], origin="lower", \
        interpolation="none")

xmin, xmax = 512-128, 512+128
ymin, ymax = 512-128, 512+128

ax[1].contour(model_image.image[xmin:xmax,ymin:ymax,0,0])

transform = ticker.FuncFormatter(lambda x,p : '%.1f"'%\
        ((x-(256.-1)/2)*0.01))

ax[1].set_xticks([27.5,77.5,127.5,177.5,227.5])
ax[1].set_yticks([27.5,77.5,127.5,177.5,227.5])
ax[1].get_xaxis().set_major_formatter(transform)
ax[1].get_yaxis().set_major_formatter(transform)

# Adjust the plot and save it.

ax[0].axis([10,5000,0,data_1d.amp.max()*1.1])

ax[0].set_xscale("log", nonposx='clip')
#ax[0].set_yscale("log", nonposx='clip')

#x0, x1 = ax[0].get_xlim()
#y0, y1 = ax[0].get_ylim()
#ax[0].set_aspect((x1-x0)/(y1-y0))

ax[0].set_xlabel("U-V Distance [k$\lambda$]")
ax[0].set_ylabel("Amplitude [Jy]")

ax[1].set_xlabel("$\Delta$RA")
ax[1].set_ylabel("$\Delta$Dec")

fig.set_size_inches((9,4))
fig.subplots_adjust(left=0.1, right=0.99, top=0.95, bottom=0.15, wspace=0.3)

# Adjust the figure and save.

fig.savefig("gaussian_model/gaussian_model.pdf")

# Close the figure.

plt.close(fig)
plt.clf()
