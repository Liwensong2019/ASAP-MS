# ASAP-MS Advion ASAP-MS Data Preprocessing

This project contains scripts and notebooks for preprocessing mass spectrometry (MS) data obtained using the Advion ASAP-MS system. It is designed for batch processing of multiple CSV files, background subtraction, peak detection, and visualisation.

# Advion ASAP-MS datapreprocessing.ipynb
  A Jupyter Notebook demonstrating the step-by-step usage of the data preprocessing functions. It includes file path handling, visualization, and example outputs.
  
# dataprocessing-bin1.py
  A Python script containing reusable functions for:
    Reading and cleaning raw ASAP-MS CSV data
    Background noise removal
    Total ion count (TIC) calculation
    Peak detection using scipy.signal.find_peaks

# Reuired packages:
pandas
numpy
matplotlib
scipy

# Notes:
The preprocessing pipeline assumes that the input .csv files are raw csv files directly exported by Advion Mass Express Software with bin = 1.
You need to adjust the dataprocessing-bin1.py when you use different bins
The notebook includes parameter tuning sections such as background frame size (bg_time) and TIC selection, which should be adjusted based on your experimental design.
