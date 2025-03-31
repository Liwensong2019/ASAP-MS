#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


# In[ ]:


def read_raw_data(file_name,bg_time):
    total_ion_count = []
    ms_data = []
    # read csvfile
    df = pd.read_csv(file_name,sep=',',header=1, low_memory=False)
    # remove experimental method
    idx = df[df.applymap(lambda x: 'Experiment Method' in str(x)).any(1)].index[0]
    df = df.iloc[:idx]
    # remove first row and first column
    df_data=df.iloc[0:,1:]
    df_data = df_data.apply(pd.to_numeric, errors='coerce')
    # calculate background spectra
    tic = df_data.sum(axis=1)
    # obtain the tic of the highest peak please adjust this based on your experiment
    select = df_data.iloc[:,360]
    total_ion_count.append(tic)
    #calculatebackgroundspectra
    bgstd = df_data.iloc[:bg_time,:].std(axis = 0)
    bg=df_data.iloc[:bg_time,:].mean(axis=0) + 3*bgstd
    file = file_name.split('/')[-1]
    return df_data, tic, bg, file,select

# In[ ]:

def find_start_points(select,find_shreshold):
    # Find peaks in the select data using find_peaks function
    peaks, _ = find_peaks(select, prominence=find_shreshold, distance=10)
    # Find the start of each peak by looking for the point where the gradient changes from positive to negative
    start_points = []
    for peak in peaks:
        start_point = peak
        while start_point > 0 and np.gradient(select)[start_point - 1] > 0:
            start_point -= 1
        start_points.append(start_point)
    # Plot the select data, peaks, and start points
#     plt.figure().set_figwidth(150)
    plt.figure(figsize=(10, 4))
    plt.plot(select)
    plt.plot(peaks, select[peaks], 'ro', label='Peaks',color='blue')
    plt.plot(start_points, select[start_points], 'rx', label='Start Points')
    for peak in start_points:
        plt.text(peak, select[peak], f'{peak}', verticalalignment='center')
    plt.xlabel('Scan Number')
    plt.ylabel('Total ion count')
    plt.legend()
    plt.show()
    return start_points,peaks


# In[ ]:

def find_and_extend_start_points(select, find_threshold,reps):
    # Step 1: Find the start points in the select data
    start_points, peaks = find_start_points(select, find_threshold)
    
    # Step 2: Keep halving the find_threshold until more than reps/2 start points are found
    while len(start_points) <= reps/2:
        find_threshold /= 1.5
        start_points = find_start_points(select, find_threshold)

    # Step 3: Ensure start_points has exactly reps elements
    if len(start_points) < reps:
        # Append the first element until the list has reps elements
        start_points += [start_points[0]] * (reps - len(start_points))
    elif len(start_points) > reps:
        # Get the select values corresponding to the start points
        select_values = [select[peak] for peak in peaks]
        # Sort start points based on the corresponding select values in descending order
        start_points = [x for _, x in sorted(zip(select_values, start_points), reverse=True)]
        # Keep only the top reps start points
        start_points = start_points[:reps]

    # Return the start_points if it has exactly 5 elements
    if len(start_points) == reps:
        return start_points
# In[ ]:
def get_user_input_insert_points(start_points):
    print("Potential insert points: ", start_points)
    user_input = input("Enter . if there is no problem or enter the desired insert points (separated by commas): ")
    if user_input == '.':
        insert_points = start_points
    else:
        insert_points = [int(x.strip()) for x in user_input.split(",")]
    return insert_points


# In[ ]:


def plot_insert_points(insert_points,tic,fig_path,file,sp_time):
    plot_points = [x + i for x in insert_points for i in range(0, sp_time)]
    plt.figure(figsize=(10, 4))
    plt.plot(tic[plot_points],'rx',color='r',label='Selected data points')
    plt.plot(tic)
    for peak in insert_points:
        plt.text(peak, tic[peak], f'{peak:}', verticalalignment='bottom')
    plt.xlabel('Scan Number')
    plt.ylabel('Total ion count')
#     plt.xlim(0,1150)
    plt.legend()
    plt.savefig(os.path.join(fig_path, file.replace('.csv','_repdetect.png')),dpi=400)


# In[ ]:


def data_process(df_data,insert_points,bg,reps,sp_time,data_path,fig_path,file):
    # calculate the raw ms data and substract background and perform normalisation
    # save the normalised data
    raw_ms_list = []
    ms_list = []
    for i in range(reps):
        sp = df_data.iloc[insert_points[i]:insert_points[i]+sp_time,:].mean(axis=0)
        ms = sp - bg
        ms [ms<0] = 0
        ms[np.isnan(ms)] = 0
        ms_list.append(ms)
    # save ms_list to csv file 
    ms_list_filename = os.path.join(data_path, f"{file.replace('.csv', '')}_rep.csv")
    pd.concat(ms_list, axis=1).to_csv(ms_list_filename, index=False, header=False)
    print('The ms_list csv file has been saved')
    # plot all reps in the same figure
    fig, axs = plt.subplots(reps, sharex=True, sharey =False, figsize = (6,12))
    for i in range(reps):
        axs[i].plot(ms_list[i].index.values.astype(float), ms_list[i].values, label = f"rep {i+1}")
        axs[i].set_xlabel('m/z')
        axs[i].set_ylabel('Ion Count')
        axs[i].legend(loc='upper right')
    plt.tight_layout()
    figname = os.path.join(fig_path, f"{file.replace('.csv', '')}_reps_sum.png")
    plt.savefig(figname, dpi=100)
    print('The plot of all reps in the same figure saved')
    # average spectrum
    avg_spectrum = sum(ms_list) / reps
#     avg_spectrum = (avg_spectrum/(avg_spectrum.sum())) *20
    # plot average spectrum
    fig, ax = plt.subplots()
    ax.plot(avg_spectrum.index.values.astype(float), avg_spectrum.values, color='red')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    ax.set_title('Average Spectrum')
    # save the figure
    figname = os.path.join(fig_path, f"{file.replace('.csv', '')}_avg.png")
    plt.savefig(figname, dpi=100)
    print('The plot of avg ms saved')
    # save average spectrum
    avg_sp_filename = os.path.join(data_path, f"{file.replace('.csv', '')}_avg.csv")
    avg_spectrum.to_csv(avg_sp_filename, index=False, header=False)
    print('avg ms csv file saved')


# In[ ]:


# Combine data in one csv file 
def combine_csv_files_by_column(input_folder, output_file,start,end):
    csv_files = [f for f in os.listdir(input_folder) if f.startswith(start)and f.endswith(end)]
    combined_data = None

    for file in csv_files:
        file_path = os.path.join(input_folder, file)
        data = pd.read_csv(file_path,header=None)

        if combined_data is None:
            combined_data = data
        else:
            combined_data = pd.concat([combined_data, data], axis=1)

    combined_data.to_csv(output_file, index=False,header=None)

