# Verne

## Usage

The core code for calculations is in the module `verne.py`. Import in python with `import verne`.


### Initialisation

Then you need to do some initialisation:

`verne.loadIsotopes()` - this loads the isotopes and density profiles for the Earth and atmosphere.

`verne.loadFFcorrections(m_x)` - this tabulates the correction factors due to Nuclear form factors. It must be reloaded for each new DM mass you want to calculate for.

### Scripts

To calculate the speed distribution at a detector, for a range of values of gamma and v, use:

```python
python CalcVelDist.py -m_x M_X -sigma_p SIGMA_P -loc LOCATION
```
where
`M_X` is the DM mass in GeV  
`SIGMA_P` is the DM-nucleon cross section in cm^2  
`LOCATION` is the detector location, "MPI" or "SUF".

Note that this will save a file into `results/veldists` and will probably take about 1 hour on a single core.

The number of signal events at a given detector is then calculated with
```python
python CalcRate_nucleus.py -m_x M_X -sigma_p SIGMA_P
```
or 
```python
python CalcRate_CDMS.py -m_x M_X -sigma_p SIGMA_P
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
