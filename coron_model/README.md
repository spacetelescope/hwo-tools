# coronagraph

A Python noise model for directly imaging exoplanets with a space based coronagraph.  

[![Build Status](https://travis-ci.org/jlustigy/coronagraph.svg?branch=master)](https://travis-ci.org/jlustigy/coronagraph) 
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/jlustigy/coronagraph/blob/master/LICENSE)

If you use this model in your own research please cite [Robinson et al (2016)](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1507.00777) and include the following acknowledgement: "This work made use of the Python coronagraph noise model, developed by J. Lustig-Yaeger and available at https://github.com/jlustigy/coronagraph/".

## Install

* Clone this github repository:
```bash
git clone git@github.com:jlustigy/coronagraph.git
```
* (optional) Add to python path to use `coronagraph` in any dirctory. This is ideal if you want to be able to access the core coronagraph functions from different project directories. If you're on a mac, then open up the `.bash_profile` file and add the following lines:
```bash
# create python path if it isn't already an environment variable
export PYTHONPATH=$PATH
# add this repo to the python path
export PYTHONPATH=${PYTHONPATH}:/Users/Name/Folder/coronagraph/
```
If you use the `.cshrc` file it might look something like this:
```shell
# For .cshrc (I know, terrible...)
setenv PYTHONPATH ${PYTHONPATH}:/Users/Name/Folder/coronagraph/
```

## Examples

#### Running `luvoir_demo.py`
You can start to get a feel for how the coronagraph noise model works by examining and running `luvoir_demo.py`. From the command line type:
```bash
python luvoir_demo.py
```
You'll get an output file that looks like this:

<img src="https://github.com/jlustigy/coronagraph/blob/master/examples/plots/luvoir_demo_.png" width="60%" height="60%" align="middle" />

Now you can create your own scripts like this to play with all the tunable parameters!

#### Simulate observation with the Integral Field Spectrograph (IFS)

```python
# Import coronagraph package
import coronagraph as cg

# Initialize Telescope, Planet, and Star objects
telescope = cg.Telescope()
planet = cg.Planet()
star = cg.Star()

# Read-in wavelength, reflectance model
model = np.loadtxt('planets/earth_quadrature_radiance_refl.dat', skiprows=8)
lam = model[:,0]            # wavelength (microns)
refl = np.pi * model[:,3]   # geometric albedo
solhr = model[:,2]          # solar flux

# Specify telescope integration time in hours
integration_time = 10.0

# Observe!
lam, dlam, Cratio, spec, sig, SNR = \
      cg.generate_observation(lam, refl, solhr, integration_time, telescope, planet, star)

```
<img src="https://github.com/jlustigy/coronagraph/blob/master/examples/plots/earth_quad_R70.png" width="100%" height="100%" align="middle" />

#### Simulate observation with the Imaging camera

```python
# Set telescope to 'Imaging' mode
telescope.mode = 'Imaging'

# Load Filter Wheel for obsevation (the default filters are the Johnson-Counsins UBVRI filters)
landsat = cg.filters.landsat()
jc = cg.filters.johnson_cousins2()

# Add Filter Wheel to Telescope
telescope.filter_wheel = jc

# Observe!
lam, spec, sig = cg.generate_observation(lam, refl, integration_time, telescope, planet, star)
```
<img src="https://github.com/jlustigy/coronagraph/blob/master/examples/plots/earth_quad_jc.png" width="100%" height="100%" align="middle" />

## Notes

* Check out the [SVO Filter Profile Service](http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=2MASS/2MASS.J&&mode=browse&gname=2MASS&gname2=2MASS#filter)
