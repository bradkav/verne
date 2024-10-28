# Verne

## Usage

### Heavy Dark Matter

The core of the standard code is in the module `verne.py`. Import in python with `import verne`.

Then you need to do some initialisation:

`verne.loadFFcorrections(m_x)` - this tabulates the correction factors due to Nuclear form factors. It must be reloaded for each new DM mass you want to calculate for.

Note that there is a hard-coded flag `NEGLECT_FF` which you can set to true to neglect corrections due to form factors (which should be valid for low-mass DM). 


### Light Dark Matter

For sub-GeV Dark Matter, the formalism is a little different. In that case, you should load the module `verne_light.py`.

In this case, you should run:

`verne_light.Tabulate_Column_Density(depth, target)` - where `depth` is the depth of the detector in metres and `target = "earth", "atmos" or "full"` depending on whether you want to include interaction with the Earth, Atmosphere or both. You should re-tabulate whenever you want to study detectors at a different depth.

**Warning:** the `verne_light` module is only valid for use with a Standard Maxwell Boltzmann distribution (as defined in `MaxwellBoltzmann.py`). It should be fine to alter the parameters for the Earth velocity, escape velocity and velocity dispersion defined in that module, but the code may not work as expected if you attempt to implement a different form for the input velocity distribution. This is because the grids of `theta` defined internally in the code are defined assuming a Maxwell Boltzmann distribution.


### Scripts

For *heavy DM*, calculate the speed distribution at a detector, for a range of values of gamma and v, use:

```python
python3 CalcVelDist.py -m_x M_X -sigma_p SIGMA_P -loc LOCATION -int INTERACTION
```
where
   - `M_X` is the DM mass in GeV  
  - `SIGMA_P` is the DM-nucleon cross section in cm^2  
  - `LOCATION` is the detector location, which can be "surface" or a number of other underground locations (e.g. "SUF" for Underground Facility, "MPI" for Max Planck Institute Munich, "MOD" for Modane, ...)
  - `INTERACTION` is the type of the type of interaction, either "SI" (Spin-independent),  "SD" (spin-dependent)

Note that this will save a file into `results/veldists` and will probably take a few minutes on a single core. The integration parameters which control the speed and precision of the calculations (`TOL` and `NSTEP`) can be updated near the top of the `verne.py` module. The script `CalcVelDist.py` is actually just a command-line wrapper for the function in `CalcVelDist_function.py`, which may be useful if you want to call this from another python script.

For *light DM*, use:
```python
python3 CalcVelDist_light.py -m_x M_X -sigma_p SIGMA_P -loc LOCATION -int INTERACTION
```
where in this case, the possible interactions are "SI" (spin-independent), "hm" (heavy dark photon mediator) and "ulm" (ultra-light dark photon mediator).

It should be straightforward to dig into the internals of the `CalcVelDist_function.py` script to adjust the grid of gamma and v that you want to use.

### Other examples

A general example (which performs an example velocity distribution calculation and event rate calculation) can be run with:
```
python3 Test.py
```
Some more specific examples, calculating the number of signal events at specific detectors, can be run with:
```python
python3 CalcRate_nucleus.py -m_x M_X -sigma_p SIGMA_P
```
or 
```python
python3 CalcRate_CDMS.py -m_x M_X -sigma_p SIGMA_P
```
depending on the detector. This will read in the relevant file from `results/veldists` and calculate the corresponding event rate and number of signal events, appending to a file in `results/Nevents`. 

The file `RunMPI_verne.py` is also provided, which helps with paralellising these calculations. It essentially launches a number of copies `CalcVelDist.py` in parallel using MPI, so that the grid scan can be performed more quickly.

## General code structure

### Isotopes

The isotopes to be considered are identified by identifiers `isoID` from `0` to `9`:

| `isoID` | Element |
| ------- | ------- |
| 0		  |	O		|
| 1		  | Si		|
| 2		  | Mg      |
| 3       | Fe      |
| 4		  | Ca		|
| 5		  | Na		|
| 6		  | S		|
| 7		  | Al		|
| 8		  | O (Atmos.) |
| 9	      | N (Atmos.) |

In addition, we calculate the form factor corrections for Lead and Copper, as these are required in calculating the shielding effects.

### Coordinate systems

**Gamma** is the angle between the mean DM velocity and the position vector of the Detector (measured from the center of the Earth). So gamma = 0 corresponds to DM particles travelling though most of the Earth before reaching the detector. While gamma = pi corresponds to DM particles coming mostly from 'overhead'.

**Theta** is the incoming direction of a particular DM particle, as measured from the position vector of the detector. Theta = 0 corresponds to particles coming directly from *below*.

### Velocity distribution

The calculation of the unperturbed velocity distribution is performed by `MaxwellBoltzmann.py`, often loaded as `MB`. In here, you can adjust the values of the escape speed and SHM speed dispersion, which are stored as module-scope variables.
