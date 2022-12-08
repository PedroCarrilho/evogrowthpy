# evogrowthpy (version 1.0.0)
This repository contains the source code of evogrowthpy, a python wrapper for a fast code to compute the scale-independent growth factor for wACDM models of interacting dark energy and for the nDGP model of modified gravity. The main computations are performed in C++ using a code based on part of [ReACT](https://github.com/nebblu/ReACT).

Authors:   Pedro Carrilho, Ben Bose<br/>
Date of last update:    December 2022<br/>

<b>Contact information:</b> If you have any questions, contact Pedro Carrilho (pedromgcarrilhoATgmail.com) or post an issue on this page.

## Quick start
### Prerequisites
In any case you need:
 * C++11 or later
 * GNU Scientific Library version 2.5 or higher (GSL; see https://www.gnu.org/software/gsl/)
 * g++ version 4.9.1 or higher

For the python wrapper, you need python with the following packages:
 * cython, numpy, scipy (essential)
 * jupyter, matplotlib (for plotting in the example notebook)

### Building and installation
The current version of evogrowthpy is meant to be used only via the python wrapper.

To install it, clone this repo, cd into the evogrowthpy directory and execute

```
pip install .
```

To check if the installation was successful, just open python and do

```
import evogrowthpy
```

A python notebook (test_evogrowthpy.ipynb) is also included, which includes an example for running the code in the most general setup. If it runs correctly, the installation was successful.

### Usage

To run the python wrapper, first import the package via
```
import evogrowthpy as evgr
```
and use the function `get_growth_wrap`
```
D, f = evgr.get_growth_wrap(Omm, hubble, [z, w0, wa, xi, Omega_rc], accuracy=1e-3, om_cb=0)
```
where `Omm` is the total matter density parameter, `hubble` is the Hubble constant divided by 100, `z` is the redshift, `w0` is the EOS parameter of dark energy and `wa` is its derivative w.r.t. `a`. `xi` is the interaction parameter of the dark scattering model and `Omega_rc` is the parameter of nDGP. The `accuracy` is set to `1e-3` which should be sufficient in most cases. The parameter `om_cb` is set to `0` by default and should only be non-zero when massive neutrinos are present (i.e. when `Om_cb != Om_m`) and the user wishes to calculate the growth rate on small scales.

The output of this function is always composed of the growth factor `D` and the growth rate `f`, with `D` normalised to be equal to the scale factor `a` at early times.

The code does not yet allow for a full calculation of the scale-dependent growth of massive neutrinos, but is sufficient to compute the growth in two limits. For the large-scale limit, the user should use the default `om_cb=0`, while for the small scale limit, the user should use instead `om_cb = Omm - Om_nu`

See the python notebook (test_evogrowthpy.ipynb) for an example of a full run and more details.

## License
evogrowthpy is free software, distributed under the GNU General Public License. This implies that you may freely distribute and copy the software. You may also modify it as you wish, and distribute these modified versions. You must always indicate prominently any changes you made in the original code and leave the copyright notices, and the no-warranty notice intact. Please read the General Public License for more details.
