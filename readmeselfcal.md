README.md
===========

Self Calibration
----------------

## This README.md is for walking you through reducing data from ALMA, VLA,etc. I will separate between the different facilities and the processes used. THis is not meant to be holistic but is meant to hey guide you through your analysis

# ALMA
## Main parts are `Preliminary`, `Self Calibration`, `Imaging Molecules and Continuum`, `Moment Maps`, `Removing sources`, `Deprojection`, `Pretty Plotting`, `Analyzing Ellipses of Images`

### Preliminary
##### I would read up as to how casa handles data and data sets. Everything is store within tables or key value pairs.This is extremely important as you could be applying work onto the model or the corrected data columns mistakenly. There are numerous columns just make sure to use the correct one.
##### Need to backup the original data file(s)
##### If needed combine all relevant data sets
##### I would first check the observation log and see if any antennas were having issues
##### Also make sure you plot the weather pattern, the antenna position, etc to get all theinformation about the observation. This can save you a lot of time later if you need to flag out an antenna due to errors or maintenance before you calibrae and realize something is off and have to redo the whole thing.
##### Use listobs to determine the Spectral Windows and the fields for observations
##### Split out the separate molecular observations based on spectral window without keeping the flags
###### This will save the split data sets to new data sets, Isuggest calling the new ms an intelligble name
##### Flag out the spectral windows and channels where molecular emission was present
###### This will construct a clean continuum ms
##### Now split out all data sets with the previous flaged data preserved. This is your new clean continuum ms.
###### You can bin the data using `width=''` which will reduce the total channel size into a more reasonable form
### Now time for full Self Calibration
##### There are several imagers modes, however I have found that the CS (Cotton Schwab) cleaning algorithm works the best for strong point sources and extended sources combined. (<http://www.cv.nrao.edu/~abridle/deconvol/node10.html>. This deals with minor and mahor sequences of subtracting out the FFT of the integrated intensity of the convolved synthetic beam
##### There are numerous weighing schemes but briggs weighing is good as you can specify a robust range between -2 and 2 from uniform to natural (balance between high res low amp to low res high amp)
##### Call the clean sequence with the proper commands specified
##### I would first call it without cleaning through the parameters and denote this as the `noselfcal` version in order to sanity check yourself
##### You will want to do 2-5 cleaning sequences depending on the brightness of the source.
##### Starting at the 0th sequence, lightly clean and move up to aggressive cleaning as you better calibrate
#### To walk through the full sequence
##### Call Clean with no cleaning to denote base image
##### call first clean with light cleaning
##### create your first phase calibration table with gaincal and plot this to verify good fit
##### then apply this phase calibration ont to the visibility file
##### Just for safety check, I would then go through the clean1 process again but cleaning to the noise floor. Call it something intelligble`p0_fullclean`. This can help tell you your progress during the self-calibration sequence
##### start the 2nd sequence with more moderate cleaning
##### rinse and repeat
##### This might seem very ambiguous, but this process is more of an art form and these is no good copy-paste recipe to follow. Everything depends on the extent, amplitude, and the quality of the observations
##### Once you reach the end of the sequence and you believe you are not increasing the amplitude of the observation anymore by correcting for the phase calibration, now time to move onto amplitude calibration
##### Now you are going to correct through calibration of the amplitude. Similar process as the phase solution howeverm you want to make sure to apply the phase and the amplitude correction at the end.
##### Generate an image with strong cleaning of the amplitude and finalized phased calibrated image.
##### I would generate and image with both 0.5, and -0.5 robust parameter to give a good baseline for comparing observations later.
##### You can now split out the data from the corrected column into a new image file.
## Imaging Molecules and continuum
### Continuum scriptforContinuumImaging.py
##### The continuum is very straightforward
##### You will generate images similar to the processes in self cal, iterating through numerous robust parameters and weighing schemes. Make sure to clean to noise floor
### Molecules scriptForLineSelfcal.py
#### General lineselfcal
##### You will want to apply your phase and amplitude calibration into each of the separate lines
##### The subtract your continuum from the line files. Careful uvcontsub natively accepts line free regions, make sure you are subtracting the right channels
##### it will be very obvious if it is correct
#### Specific line imaging
##### I would generate a separate file per molecule you are imaging
##### Now you will iterate through, similar to the continuum imaging, you are iterating through different weighing schemes and robust parameters
##### You will also want to change the tapering and uv ranges to exclude small and large baselines if certain antenna have worse sampling errors
## Moment Maps
#####
