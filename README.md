# mkipp
Kippenhahn plotter for MESA. BEWARE: Currently undergoing big changes, this README might be out of date!

This repository provides two files:
- mesa_data.py: a simple parser of MESA data files that only parses specific columns to save time
- mkipp.py: the Kippenhahn plotter
 
Although lower-level functions are provided that allow plotting directly into a matplotlib axis object, a high level function produces a fully decorated plot. The repository comes with some sample MESA data for a 20 Msun star, which if used with the following python script

```python
import mkipp
mkipp.decorated_kipp_plot(logs_dir = "LOGS", profile_numbers = range(1,45,1),\
   contour_plots = ["eps_nuc"], core_masses = ["He","C","O"], \
   xaxis = "model_number", time_units = "Myr", \
   save_file = False, mass_tolerance = 1e-10)
```
results in the following plot

![kipp](example_kipp.png)

Switching to xaxis = "star_age" can be used to create a plot with age as the xaxis. The time_units options which can have the values "yr", "1000 yr", "Myr" and "Gyr" controls the unit of time.
