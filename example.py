#!/usr/bin/env python
import mkipp
import matplotlib.pyplot as plt
mkipp.decorated_kipp_plot(logs_dir = "LOGS", profile_numbers = range(1,45,1), contour_plots = ["eps_nuc"], core_masses = ["He","C","O"], \
        xaxis = "model_number", time_units = "Myr", save_file = False, mass_tolerance = 1e-10)
