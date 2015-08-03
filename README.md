# mkipp
Kippenhahn plotter for MESA. BEWARE: Currently undergoing big changes, this README might be out of date!

This repository provides two files:
- mesa_data.py: a simple parser of MESA data files that only parses specific columns to save time
- kipp_data.py: process data for Kippenhan plot
- mkipp.py: the Kippenhahn plotter
 
Although lower-level functions are provided that allow plotting directly into a matplotlib axis object, a high level function produces a fully decorated plot. The repository comes with some sample MESA data for a 20 Msun star, which if used with the following python script

```python
import mkipp
mkipp.kipp_plot(mkipp.Kipp_Args())
```
results in the following plot

![kipp](Kippenhan.png)

Switching to xaxis = "star_age" can be used to create a plot with age as the xaxis. The time_units options which can have the values "yr", "1000 yr", "Myr" and "Gyr" controls the unit of time.
