## Dependencies

The code is compatible with Python2.7. Requires [numpy](http://www.numpy.org) and [scipy](https://www.scipy.org).

## To-Do

- Make sure that all cross section calculations are the full cross section (not the high mass limits)... also check other parts which might depend on the DM mass...  
- Check gamma values off-axis  
- Check all `BJK` in the code - these are comments to me...  
- Check that the notation mx, mX, etc. is consistent  
- Check out CRESST sapphire
- For the HPC code, make sure to copy the parameters first maybe!!!

## General code structure

`verne.py` is where most of the heavy-duty calculations take place.

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


### Mapping Velocities

`VA.calcVfinal_full` calculates the final velocity of the DM particle for a given initial velocity, incident angle and depth (taking into acocunt either the atmosphere, or the Earth, or both). 

`VA.calcVfinal_shield` takes into account attenuation from shielding (a single straight line attentuation).

### Velocity distribution

The calculation of the unperturbed velocity distribution is performed by `MaxwellBoltzmann.py`, often loaded as `MB`. In here, you can adjust the values of the escape speed and SHM speed dispersion, which are stored as module-scope variables.

### Calculating the velocity distribution

`VA.CalcFullF` calculates the speed distribution at the detector, assuming only a single isotope (`ID`).

`VA.CalcFullF_full` calculates the speed distribution, taking into account all isotopes. You can specify `target` = `atmos` or `earth` (or anything else just includes both). This does not currently include the contribution from the detector shielding. This should be added separately later, once everything is running stably.

The calculation of the velocity distribution at speed v\_final is done by evaluating the unperturbed Maxwell-Boltzmann distribution at v\_initial, the speed which a particle must have had initially in order to reach your detector as a speed v\_final.

**NB:** Need to be careful about the theta-cut procedure when you're only looking at scattering off the atmosphere. In particular, the maximum speed at the detector is *not* a monotonic function of the incoming DM angle (because the distance travelled is not a monotonic function of the angle either).

## Calculation scripts

**CalcVelDist\_grid\_single.py** - Calculate the speed distribution at the detector, including on the Earth stopping and assuming a single element in the Earth (usually Silicon...).

**CalcVmap.py** - Calculate tabulated values of the final DM velocity, as a function of the initial DM velocity and incoming direction.

**CalcVelDist.py** - Calculate the velocity distribution (taking into account all effects).

**CalcMissingDM.py** - Test script for working out how to calculate the amount of DM that is 'missing' from the distribution. That is, particles which end up being too slow to be above the analysis cut off (v\_cut ~ 1 km/s).

## Initialisation

Need to run the following:

`VA.loadPhiInterp()` - this loads in some numerical results which are required for some intermediate integrals over the velocity distribution. *(This is now done inside the MaxwellBoltzmann.py module)*

`VA.loadIsotopes()` - this loads the isotopes and density profiles for the Earth and atmosphere.

`VA.loadFFcorrections(m_x)` - this tabulates the correction factors due to Nuclear form factors. It must be reloaded for each new DM mass you want to calculate for.

`VA.calcDeff_interp(depth, isoID)` - *(deprecated)* this tabulates the path length to the detector as a function of the incoming direction (for a fixed Earth element, identified by isoID). This is only used by `VA.CalcFullF` and other parts of the code which assumes a single Earth isotope.

`VA.LoadVelocityFunctions()` - *(deprecated)* this tabulates the Velocity Functions, which map the initial to final velocity for a given isotope. This framework only makes sense in the case when the Form Factor corrections are *not* taken into account. However, the Form Factor corrections depend on the DM speed, so the full ODE must be solved and the simple interpolation functions used here cannot be applied.




## Plotting scripts

**PlotCorrectionFactors.py** - Plot correction factors (due to nuclear form factors) as a function of v (for a couple of example nuclei).

**PlotAlpha.py** - Plot alpha (the latitude of mean incoming DM) as a function of time.

**PlotD.py** - Plot D, the distance travelled underground (or through the atmosphere) for a detector at a particular depth.

**PlotGamma\_STAN.py**, **PlotGamma\_SURF.py** - Plot gamma, the mean DM incidence angle for a given lab - either Stanford, or SURF.

**PlotVelocityTransfer\_all.py** - Plot the velocity transfer function (i.e. v\_final as a function of v\_initial), decomposed sequentially in terms of the atmosphere, the Earth and the shielding.

**PlotVelocityTranfer\_elements.py** - Plot the velocity transfer function in the Earth, assuming individual elements. *Doesn't really make sense, since the ODE is highly non-linear.*

**PlotVelDists.py** - Plot the DM speed distributions at the detector, from tables of final velocities, calculated using **CalcVmap.py**.

**PlotCutOffAngle.py** - Checking and plotting the final speed which is obtained by a particle of initial speed vesc. This is done as a function of incoming angle theta, in order to determine which angle of theta we should cut off at.