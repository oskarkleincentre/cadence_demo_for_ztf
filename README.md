# Cadence Demonstration for the Zwicky Transient Facility (ZTF)

A piece of code to demonstrate the ztf cadence by:
 -  generating images of the entire sky at different times
 - showing background positions like the rough extent of the Milky Way
 - Areas of the sky that are desirable for pointing at the given instance for SN science
 - The ZTF field of view, the rough patch it covers and the filter used at a viewing if taken from a scheduler.
 - Showing simulated SNIa in the background (and potentially showing the ZTF capabiltiy in detecting and characterizing them) 

Animations can then be created from these images.

## Current Status:

Currently, we have an implementation which uses outputs form external programs to show the entire sky, a statistical sample of SNIa in the low z (0.0-0.2) universe and their changes in brightness, an approximate region of the Milky Way, the variation of 
 position of the ZTF field of view and filter, and a region of the sky expected to be ideal for studying SN at a particlar time. A brief of explanation of how different components can being currently used here is explained below:

-  __Motion of the ZTF Field of View and Filter Changes__   we are using examples of ZTF schedules (from E. Bellm) to illustrate the motion of the ZTF field of view and changing of filters.
- __The Milky Way Region__ This is used in conjuction with an approximate outline of the region occupied by the Milky Way, which obscures the extragalctic universe. This is a region between +20 and -20 degrees from the disk of the Milky Way. 
- __Visible Fields__  Ideal locations in the sky to study supernovae are also affected by factors like airmass, location of the moon, location of the Milky Way Disc etc.  The location of such fields on the basis of such conditions are computed using a separate program (Goobar, Nugent), and the outputs from this can be conveniently read in and converted to Healpixels of sizes similar to the field of view. These can then be displayed as polygons.
- __SNIa__ A supernovae Type Ia  catalog for the area of ZTF over a redshift range of 0.0-0.2 (higher than the max redshift of the planned ZTF cosmology sample) is used (but their diversity is neglected) and their shown with a color corresponding to redshift, and point sizes related to the apparent magnitude at the time of the observation.
 
Sets of components above may be chosen as the sets to be shown in a particular visualization.

## Using this pacakge
In order to use this package, please install the software in the [`requirements file`](./master/requirements.md) along with their dependencies. Installing this package after that involves running
```python setup.py install```
from the top level directory.
 In order to show the SN, and the visible fields from external sources, you will also need to have the SN catalog, and a directory containing the computed visible fields.

To get started, please use the demo jupyter notebooks in `examples`.

Some animations made with this package are listed in [animations](./data/animations_released.md)
## Future Plans :

We would like to integrate futher functionality into this program:
- Calculate the visible fields independently in this program from a set of conditions, and track the position of the moon.
- Introduce additional methods and potentially subplots to track SN that meet certain conditions cumulatively. Examples of such conditions may be being detected, or well characterized. 

Comments and suggestions on these are welcome through the issues page.

## Contributors
- R. Biswas
- U. Feindt
- E. Bellm
- A. Goobar
- M. Bulla
- S. Dhawan
- L. Hangard
