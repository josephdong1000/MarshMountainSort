#!/usr/bin/env python
# coding: utf-8

# # Mountain Sort - Marsh Version

# Remember to install all of the packages used in this script. I installed all of them locally in my home directory '/~'
# 
# To run without crashing, try the following:
# - Export as a .py file and run on HPC with `python marshmountainsort.py`
# - Comment out blocks that render large areas of signal, since this may overflow RAM and crash the kernel
# - Allocate more memory for the kernel
#     - Generate a config file with `jupyter lab --generate-config`
#     - Open up the config file in `~/.jupyter/jupyter_lab_config.py`. You may need to show hidden files with Ctrl+H
#     - Uncomment and edit the line `c.ServerApp.max_buffer_size = 1000000000`. This raises the buffer size to 1GB. If there are still issues, keep raising buffer size.
#     - Restart Jupyter

# ## Run Code on HPC

# ### Run 1 File

# ### Generate All Binaries

# ### Generate All Figures

# ## Analysis Constants

# In[2]:


# If True, uses rec_folders to get file names and run scripts. Code will ignore rec_folder and tetrode_name.
# Set True if generating all .bins or figures.
# Remember to reset filecount.pk! See below scripts
GENERATEALL = True

# If you want to convert .txt into .bin (Warning: time consuming!!! ~4hrs)
# Set True when you want to create/regenerate .bins. If so, run on HPC to save time
PARSENEW = True
# Saves converted .txt files to .npy.gz instead of .bin (to save space). Should always be True
SAVECOMP = True
# If a .npy.gz is present alongside .bin, use the .npy.gz
FAVORCOMP = True

# Make figures. Turn off if you just want to convert TXTs for instance
MAKEFIGS = False
# Save all figures to folder /output (may overwrite)
SAVEFIGS = False

# Determines if sorting 1 tetrode or all tetrodes at once. Should always be True.
SINGLETETRODE = True


# In[3]:


from pathlib import Path

# Name of the TXT folder containing 1 .txt file with tetrodes. Ignored if GENERATEALL
# rec_folder = '1133_WT_4 half turns + 3-8th of full turn_right after turning_01-25-17' # For visualization testing
rec_folder = "1176_WT Bl6_5 half turns + 3-8th of full turn_right after turning_02-09-17"

# Choose between analyzing ca3, ca1s, or ca1o. Ignored if GENERATEALL
# tetrode_name = ''
tetrode_name = ['ca3', 'ca1s', 'ca1o'][2]


# GENERATEALL: Name of every binary/TXT file. Can be pulled with findnewtxt.ipynb
rec_folders = ["1238_Exp_3 half turn after turning_5_3_17",
                "766_WT_Final Protocol Morning 6th half turn +1_4 _recorded_ after turning_7-28-2017",
                "766_WT_5th half turn_recorded after turning_7-26-2017",
                "1269_EXP_4th half turn_recorded after recovery_7-31-2017",
                "766_WT_Final Protocol Night 6th half turn +1_4 _recorded after turning_7-28-2017",
                "766_WT_Final Protocol Morning 6th half turn_recorded after turning_7-27-2017",
                "1244_Bl6_5th half turn_recorded after turning_ 7-02-2017",
                "1269_EXP_Final Protocol 6th half + 1_4th quarter turn Morning_recorded after recovery_8-3-2017",
                "1236_Exp_6 half turn after turning_Day1 night 5_5_17",
                "766_WT_3rd half turn_recorded after turning_7-24-2017",
                "1244_Bl6_Final Protocol_6th half turn_recorded after turning_Day1 morning 7-03-2017",
                "1269_EXP_5th half turn_recorded after recovery_8-1-2017",
                "1269_EXP_3rd half turn_recorded after recovery_7-28-2017",
                "1269_ Exp_2nd 12 hour recording_8-9-17_Spike Wave-0 - Copy",
                "1236_Exp_3 half turn after turning_5_3_17",
                "766_WT_Final Protocol Night 6th half turn+2_4_recorded after turning_7-31-2017",
                "1269_EXP_2nd half turn_recorded after recovery_7-27-2017",
                "1244_Bl6_1st half turn_recorded after turning_6-28-2017",
                "1236_Exp_4 half turn after turning_5_4_17",
                "766_WT_Final Protocol Morning 6th half turn+3_4_recorded after turning_8-02-2017",
                "1269_EXP_Final Protocol 6th half turn Morning_recorded after recovery_8-2-2017",
                "1236_Expanded_2 half turns_after turning_05-02-17",
                "1269_EXP_Final Protocol 6th half + 1_4th quarter turn Night_recorded after recovery_8-3-2017",
                "1244_Bl6_2nd half turn_recorded after turning_ 6-29-2017",
                "1244_Bl6_4th half turn_recorded after turning_ 7-01-2017",
                "1236_Exp_6 half turn after turning_Day1 morning 5_5_17",
                "1238_Exp_5 half turn after turning_5_5_17",
                "1238_Exp_4 half turn after turning_5_4_17",
                "766_WT_Final Protocol Night 6th half turn_recorded after turning_7-27-2017",
                "1269_EXP_Final Protocol 6th half turn Night_recorded after recovery_8-2-2017",
                "1236_Exp_5 half turn after turning_5_5_17",
                "766_WT_Final Protocol Morning 6th half turn+2_4_recorded after turning_7-31-2017",
                "1238_Exp_6 half turn after turning_ Day 1 morning 5_5_17",
                "1269_ Exp_3rd 12 hour recording_8-10-17_Spike Wave-0 - Copy",
                "1238_Expanded_2 half turns_after turning_05-02-17",
                "1244_Bl6_3rd half turn_recorded after turning_ 6-30-2017",
                "1238_Exp_6 half turn after turning_ Day 1 night 5_5_17",
                "766_WT_4th half turn_recorded after turning_7-25-2017",
                "766_WT_2nd half turn_recorded after turning_7-21-2017"]
print(len(rec_folders))

# GENERATEALL: Name of every tetrode position
tetrode_names = ['ca3', 'ca1s', 'ca1o']


# Define Paths and folders used in analysis.
base_folder_path = Path(f'/mnt/isilon/marsh_single_unit/MarshMountainSort')
raw_txts_folder = 'rawtxts'
binary_folder = 'bins'
temporary_folder = 'temp'
output_folder = 'output'
filecount_name = 'filecount.pk'


# ## Import Dependencies

# In[4]:


# %matplotlib widget
# import pexpect
from tempfile import TemporaryDirectory
import numpy as np
import pkg_resources
pkg_resources.require('matplotlib==3.8.3')
import matplotlib.pyplot as plt
import os
import shutil
import warnings
import time
import csv
import sys
import pickle
import tempfile
import gzip
# from IPython import get_ipython

import mountainsort5 as ms5
from mountainsort5.util import create_cached_recording
import spikeinterface.core as si
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre
# import spikeinterface.sorters as ss
import spikeinterface.postprocessing as spost
# import spikeinterface.qualitymetrics as sqm
# import spikeinterface.exporters as sexp
# import spikeinterface.comparison as scmp
# import spikeinterface.curation as scur
import spikeinterface.sortingcomponents as sc
import spikeinterface.widgets as sw
# import pyedflib
import sonpy
import zugbruecke
import probeinterface as pi
from probeinterface.plotting import plot_probe_group, plot_probe
import ffmpeg
# plt.rcParams(['animation.ffmpeg_path']) = ''
plt.rcParams['figure.dpi'] = 150

# import matplotlib.animation as manimation
# manimation.writers.list()

# import sortingview
# import pyvips
# import ephyviewer


# ## `GENERATEALL`: Get Filename and Step Filecount 

# In[5]:


fcount_file_path = base_folder_path / filecount_name
print(f'Path for filecount:\n\t{fcount_file_path}')

def reset_filecount():
    with open(fcount_file_path, 'wb') as file:
        pickle.dump(0, file, protocol=pickle.HIGHEST_PROTOCOL)

def step_filecount():
    with open(fcount_file_path, 'rb') as file:
        out = pickle.load(file)
    out += 1
    with open(fcount_file_path, 'wb') as file:
        pickle.dump(out, file, protocol=pickle.HIGHEST_PROTOCOL)
        
def get_filecount():
    with open(fcount_file_path, 'rb') as file:
        out = pickle.load(file)
        return out

def get_file(n):
    return rec_folders[n // 3], tetrode_names[n % 3]
    


# In[125]:


##### DO NOT run this cell after resetting filecount. This will step filecount by 1! #####
# If generating all figures, read folder and tetrode name and step filecount
print("=" * 50)
if GENERATEALL:
    rec_folder, tetrode_name = get_file(get_filecount())
    print(f"\tFILECOUNT: {get_filecount()}")
    if PARSENEW:
        for i in range(len(tetrode_names)):
            step_filecount() # Skip 3 to next file, ignoring tetrodes
    else:
        step_filecount()
        

print(f"\tREC_FOLDER: {rec_folder}\n\tTETRODE_NAME: {tetrode_name}")
print("=" * 50)


# ## Load .txt and Convert to .bin

# In[9]:


txt_folder_path = base_folder_path / raw_txts_folder / rec_folder
bin_folder_path = base_folder_path / binary_folder / rec_folder
temp_folder_path = base_folder_path / temporary_folder
output_folder_path = base_folder_path / output_folder / rec_folder

txt_file = None
for file in os.listdir(txt_folder_path):
    if file.endswith(".txt"):
        txt_file = Path(file)
        break
if txt_file is None:
    warnings.warn("No .TXT file found!")
    txt_file = ""

txt_file_path = txt_folder_path / txt_file
bin_file_path = bin_folder_path / f"{rec_folder}.bin"
bin_file_path_ca3 = bin_folder_path / f"{rec_folder}_ca3.bin"
bin_file_path_ca1s = bin_folder_path / f"{rec_folder}_ca1s.bin"
bin_file_path_ca1o = bin_folder_path / f"{rec_folder}_ca1o.bin"

# Report the PLX file
print(f'TXT to be used:\n\t{txt_file_path}')
print(f'Path for BIN to be saved:\n\t{bin_file_path}')
print(f'Path for cached files to be saved:\n\t{temp_folder_path}')
print(f'Path for image outputs:\n\t{output_folder_path}')


# In[45]:


def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

if PARSENEW:
    tstart = time.time()
    with open(txt_file_path, "r", encoding="utf-8", errors='ignore') as f:
        nlines = sum(bl.count("\n") for bl in blocks(f))
    tend = time.time()
    print(f"Time elapsed: {tend - tstart} seconds")


# In[46]:


# Read TXT file and convert to Binary

ncol = 16 # Number of channels in full recording, not including time column
ncol_tet = 4 # Number of channels in 1 tetrode
nheader = 12 # Number of lines in Dwave block header
f_s = 25000

if PARSENEW:
    ncol_arr = np.arange(ncol + 1) # including time column
    nrblocks = nlines / (f_s + nheader) # 1 recording block has header + 1s of recording data
    print(f"Number of recording blocks: {nrblocks}")
    nrblocks = round(nrblocks)


# In[47]:


if PARSENEW:
    txtblocks = []
    tstart = time.time()
    
    for i in range(nrblocks):
        txtblocks.append(np.loadtxt(txt_file_path, dtype=np.float32, delimiter='\t', usecols=ncol_arr, skiprows=(f_s+nheader)*i+nheader, max_rows=f_s))
        if i % 50 == 0:
            print(f"Processing block {i}..")
    
    txtdata = np.concatenate(txtblocks, axis=0)
    print(txtdata.shape)

    tend = time.time()
    print(f"Time elapsed: {tend - tstart} seconds")


# In[48]:


print(bin_folder_path)
if PARSENEW:
    # if os.path.exists(bin_folder_path):
    #     shutil.rmtree(bin_folder_path, ignore_errors=True)
    os.makedirs(bin_folder_path, exist_ok=True)


# In[49]:


def save_comp(arr, filepath):
    fcomp_name = filepath.parent / f"{filepath.stem}.npy.gz"
    with gzip.GzipFile(fcomp_name, "w") as fcomp:
        np.save(file=fcomp, arr=arr)
    print(f"Compressed: {fcomp_name.name}")


if PARSENEW:
    tstart = time.time()
    txtdata_notime = txtdata[:, 1:] # Discard time

    # Reorder data because of weird Datawave order convention.
    # Datawave sorts alphabetically, so 10-15 are before 2 and after 1
    # permutation = [0, 1, 10, 11, 12, 13, 14, 15, 2, 3, 4, 5, 6, 7, 8, 9]
    # idx = np.empty_like(permutation)
    # idx[permutation] = np.arange(len(permutation))
    # txtdata_notime[:] = txtdata_notime[:, idx]

    print(txtdata_notime.shape)
    print("Saving files..")
    # txtdata_notime = np.transpose(txtdata_notime)
    
    txtdata_ca3 = txtdata_notime[:, 0:4]
    txtdata_ca1o = txtdata_notime[:, 4:8]
    txtdata_ca1s = txtdata_notime[:, 10:14]

    if SAVECOMP:
        save_comp(txtdata_notime, bin_file_path)
        save_comp(txtdata_ca3, bin_file_path_ca3)
        save_comp(txtdata_ca1o, bin_file_path_ca1o)
        save_comp(txtdata_ca1s, bin_file_path_ca1s)
    else:
        txtdata_notime.tofile(bin_file_path) # Will save in row major order 'C'
        txtdata_ca3.tofile(bin_file_path_ca3)
        txtdata_ca1o.tofile(bin_file_path_ca1s)
        txtdata_ca1s.tofile(bin_file_path_ca1o)

    tend = time.time()
    print(f"Binary files saved. Time elapsed: {tend - tstart} seconds")
    print("=" * 50)


# In[50]:


if not MAKEFIGS:
    print("MAKEFIGS is false. Finishing program.")
    time.sleep(3)
    sys.exit()


# ## Define Probe Groups

# In[51]:


pg = pi.ProbeGroup()

tetrode_x = np.array([-2500, 2000, -2000, 2500 ]) # um
tetrode_y = np.array([-480, 0, 0, -480]) # um
tetrode_z = np.array([-1000, 0, 0, -1000]) # um
global_device_channel_indices = [0, 1, 2, 3, # corresponds to Datawave channel
                                 4, 5, 6, 7,
                                 10, 11, 12, 13,
                                 -1, -1, -1, -1]
plot_xlim = np.array([-3000, 3000])
plot_ylim = np.array([-1000, 1000])
plot_zlim = np.array([-1500, 500])

### TEST MULTIPLIER. VISUALIZATION PURPOSES ONLY
testmult = None
if testmult is not None:
    warnmessage = f'testmult is not None. Distances are {testmult}x smaller!'
    warnings.warn(warnmessage)
    tetrode_x = tetrode_x * testmult
    tetrode_y = tetrode_y * testmult
    tetrode_z = tetrode_z * testmult
    plot_xlim = plot_xlim * testmult
    plot_ylim = plot_ylim * testmult
    plot_zlim = plot_zlim * testmult
### TEST MULTIPLIER

fig, ax1 = plt.subplots(1, 1, subplot_kw={'projection': '3d'})
fig, ax2 = plt.subplots(1, 1)
for i in range(4):
    tetrode = pi.generate_tetrode()

    if not SINGLETETRODE:
        tetrode_3d = tetrode.to_3d('xy')
        tetrode_3d.move([tetrode_x[i], tetrode_y[i], tetrode_z[i]])
        tetrode_3d.set_device_channel_indices(global_device_channel_indices[i*4:i*4+4])
        tetrode_3d.set_contact_ids(np.arange(i*4, i*4+4, 1))
        tet_add = tetrode_3d
        
        tetrode_3dproj = tetrode_3d.to_2d('xy')
        tetrode_3dproj.set_device_channel_indices(global_device_channel_indices[i*4:i*4+4])
        tetrode_3dproj.set_contact_ids(np.arange(i*4, i*4+4, 1))
        
        plot_probe(tet_add, ax=ax1)
        plot_probe(tetrode_3dproj, ax=ax2, with_device_index=True, with_contact_id=True)
    else:
        tetrode.set_device_channel_indices(np.arange(4))
        tetrode.set_contact_ids(np.arange(4))
        tet_add = tetrode
        
        plot_probe(tet_add, ax=ax2, with_device_index=True, with_contact_id=True)

    pg.add_probe(tet_add)

    if SINGLETETRODE:
        break

print(pg.get_global_device_channel_indices())
if not SINGLETETRODE:
    # ax1.set_xlim(plot_xlim)
    # ax1.set_ylim(plot_ylim)
    # ax1.set_zlim(plot_zlim)
    ax2.set_xlim(plot_xlim)
    ax2.set_ylim(plot_ylim)
else:
    ax1.remove()
plt.show()


# ## Load .bin

# In[52]:


# Load binaries. If SINGLETETRODE = False, load all tetrodes regardless of tetrode_name
if not SINGLETETRODE:
    full_raw_rec = se.read_binary(bin_file_path, sampling_frequency=f_s, dtype=np.float32, num_channels=ncol, gain_to_uV=1000)
else:
    if tetrode_name == 'ca3':
        tet_bin_path = bin_file_path_ca3
    elif tetrode_name == 'ca1s':
        tet_bin_path = bin_file_path_ca1s
    elif tetrode_name == 'ca1o':
        tet_bin_path = bin_file_path_ca1o
    else:
        raise Exception(f'Unrecognized tetrode_name: {tetrode_name}')

    tet_fcomp_path = tet_bin_path.parent / f"{tet_bin_path.stem}.npy.gz"
    if tet_fcomp_path.is_file() and FAVORCOMP:
        print("Using .NPY.GZ file..")
        with tempfile.NamedTemporaryFile(dir=tet_bin_path.parent) as tmp:
            fcomp = gzip.GzipFile(tet_fcomp_path, "r")
            bin_arr_recov = np.load(fcomp)
            bin_arr_recov.tofile(tmp)

            full_raw_rec = se.read_binary(tmp.name, sampling_frequency=f_s, dtype=np.float32, num_channels=ncol_tet, gain_to_uV=1000)
            full_raw_rec = full_raw_rec.set_probes(pg)
    else:
        print("Using .BIN file..")
        full_raw_rec = se.read_binary(tet_bin_path, sampling_frequency=f_s, dtype=np.float32, num_channels=ncol_tet, gain_to_uV=1000)
        full_raw_rec = full_raw_rec.set_probes(pg)

print(full_raw_rec)
print(full_raw_rec.get_probegroup())
full_raw_rec.get_probegroup().to_dataframe(complete=True)


# In[53]:


recording = full_raw_rec # load your recording using SpikeInterface
freq_min = 400 # Hz
freq_max = 8064 # Hz

# Preprocess recording
recording_prep = recording
recording_prep = spre.common_reference(recording_prep, dtype=np.float32)
recording_prep = spre.scale(recording_prep, gain=10, dtype=np.float32) # Scaling for whitening to work properly
recording_prep = spre.whiten(recording_prep, dtype=np.float32)
recording_prep = spre.bandpass_filter(recording_prep, freq_min=freq_min, freq_max=freq_max, dtype=np.float32)
recording_preprocessed = recording_prep

# Separately, downsample recording for plotting
downsamp_factor = 25
recording_downsamp = recording
recording_downsamp = spre.bandpass_filter(recording_downsamp, freq_min=freq_min, freq_max=freq_max)
recording_downsamp = spre.resample(recording_downsamp, resample_rate=f_s // downsamp_factor, dtype=np.float32)
print(f"Downsampled f_s for plotting: {f_s // downsamp_factor} Hz")


# ## Sort Spikes

# In[54]:


print(temp_folder_path)

if os.path.exists(temp_folder_path):
    shutil.rmtree(temp_folder_path, ignore_errors=True)
os.makedirs(temp_folder_path, exist_ok=True)

try:
    with TemporaryDirectory(dir=temp_folder_path, prefix='temp-') as tmpdir:
        print(tmpdir)
        
        # cache the recording to a temporary directory for efficient reading
        recording_cached = create_cached_recording(recording_preprocessed, folder=tmpdir)
    
        # use scheme 1
        # sorting = ms5.sorting_scheme1(
        #     recording=recording_cached,
        #     sorting_parameters=ms5.Scheme1SortingParameters(detect_channel_radius=100)
        # )
    
        # or use scheme 2
        sorting = ms5.sorting_scheme2(
            recording=recording_cached,
            sorting_parameters=ms5.Scheme2SortingParameters(
                phase1_detect_channel_radius=100,
                detect_channel_radius=100,
                # detect_sign=0,
            )
        )
    
        # # or use scheme 3
        # sorting = ms5.sorting_scheme3(
        #     recording=recording_cached,
        #     sorting_parameters=ms5.Scheme3SortingParameters(...)
        # )
    
        shutil.rmtree(tmpdir, ignore_errors=True)
    
    # Now you have a sorting object that you can save to disk or use for further analysis

except OSError: # Weird bug with clearing the temporary directory
    pass


# ## Plot and Save Figures

# In[55]:


if sorting.get_num_units() == 0:
    warnings.warn("No units found!")
    if SAVEFIGS:
        os.makedirs(output_folder_path / f"{tetrode_name}_NO_UNITS", exist_ok=True)
    sys.exit()


# In[ ]:


print(sorting)
print(output_folder_path)

num_units = sorting.get_num_units()

if SAVEFIGS:
    # if os.path.exists(output_folder_path):
    #     shutil.rmtree(output_folder_path, ignore_errors=True)
    os.makedirs(output_folder_path, exist_ok=True)

time_plotsample = 0 # s
if SAVEFIGS:
    time_plotsample_range = recording_preprocessed.get_duration()
else:
    time_plotsample_range = 50 # s

print(f"Plotting from {time_plotsample} s to {time_plotsample_range} s")


# In[ ]:


# Create filtered recording for waveform extraction
recording_temps = recording
recording_temps = spre.bandpass_filter(recording_temps, freq_min=400, freq_max=8064, dtype=np.float32)

# Preview traces to be plotted
sw.plot_traces(recording_downsamp, # Plot downsampled recording to not crash kernel 
               time_range=[0, 20], 
               show_channel_ids=True,
              )


# In[ ]:


# Extract waveforms from recording based on sorting
we = si.extract_waveforms(recording=recording_temps, # Extract from filtered signals
                          sorting=sorting, 
                          folder=base_folder_path / "waveforms_dense", 
                          sparse=False, 
                          overwrite=True,
                          allow_unfiltered=True)
print(we)

# Calculate metrics for plotting
corr = spost.compute_correlograms(we)
pcs = spost.compute_principal_components(we)
amps = spost.compute_spike_amplitudes(we)
sloc = spost.compute_spike_locations(we, method="monopolar_triangulation")
tmet = spost.compute_template_metrics(we, include_multi_channel_metrics=False)
sim = spost.compute_template_similarity(we)
uloc = spost.compute_unit_locations(we, method="monopolar_triangulation")


# In[ ]:


# Set a unit coloring scheme
unit_cmap = sw.get_unit_colors(sorting, 
                               map_name='gnuplot2',
                               # color_engine='matplotlib',
                               # color_engine='colorsys',
                               # color_engine='distinctipy',
                              )


# In[ ]:


fig, ax = plt.subplots(2, 1, figsize=(10,5), sharex=True, gridspec_kw={'height_ratios':[1, 0.3]})
plt.subplots_adjust(hspace=0)
backend = "matplotlib"

tstart = time.time()

sw.plot_traces(recording_downsamp, # Plot downsampled recording to not crash kernel 
               time_range=[time_plotsample, time_plotsample + time_plotsample_range], 
               show_channel_ids=True,
               backend=backend,
               # color_groups=True,
               # color='C0',
               ax=ax[0],
              )
ax[0].set_ylabel("Datawave Channel")
ax[0].get_legend().remove()
sw.plot_rasters(sorting, 
                time_range=[time_plotsample, time_plotsample + time_plotsample_range], 
                backend=backend,
                color='black',
                ax=ax[1],
               )
ax[1].set_ylabel("Unit Index")
ax[1].set_xlabel("Time (s)")
fig.suptitle("Raster Plot", y=0.95)

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_raster', bbox_inches='tight')
plt.show()

tend = time.time()
print(f"Time elapsed: {tend - tstart} seconds")


# In[ ]:


# widg = sw.plot_spikes_on_traces(we,
#                                 # time_range=[time_plotsample, time_plotsample + time_plotsample_range],
#                                 time_range=[216, 217],
#                                 unit_colors=unit_cmap,
#                                 show_channel_ids=True,
#                                 figsize=(10,5),
#                                )
# widg.ax.set_title("Spikes on Traces")
# widg.ax.set_ylabel("Channel ID")

# if SAVEFIGS:
#     plt.savefig(output_folder_path / f'{tetrode_name}_spiketraces', bbox_inches='tight')
# plt.show()


# In[ ]:


widg = sw.plot_template_similarity(we, display_diagonal_values=True, figsize=(1.5 * num_units, 1.5 * num_units))

for (j, i), label in np.ndenumerate(sim):
    widg.ax.text(i, j, round(label, 3), ha='center', va='center', bbox=dict(facecolor='white'))
widg.ax.set_title("Template Similarity")

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_tempsim', bbox_inches='tight')
plt.show()


# In[ ]:


widg = sw.plot_amplitudes(we, plot_histograms=True, unit_colors=unit_cmap, figsize=(10, 5))
widg.axes[0][0].set_title("Spike Amplitudes")
widg.axes[0][1].set_title("Amplitude Histogram")

widg.axes[0][0].patch.set_edgecolor('black')
widg.axes[0][0].patch.set_linewidth(1)
widg.axes[0][1].patch.set_edgecolor('black')
widg.axes[0][1].patch.set_linewidth(1)
# widg.axes[0][1].spines["right"].set_visible(True)

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_spikeamps', bbox_inches='tight')
plt.show()


# In[ ]:


widg = sw.plot_unit_depths(we, figsize=(10,3), depth_axis=0, unit_colors=unit_cmap)
# widg.figure.set_figheight(10)
widg.ax.set_ylim(-50, 50)
widg.ax.set_title("Unit Depth and Amplitude")

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_depthamp', bbox_inches='tight')
plt.show()


# In[ ]:


# fig, ax = plt.subplots(1, 1, figsize=(8,5))
# print(we.channel_ids)

widg = sw.plot_unit_templates(we, 
                       # channel_ids=np.array([[0]]),
                       # channel_ids=[1],
                       # unit_ids=[1],
                       unit_colors=unit_cmap,
                       scale=3e15,
                       # backend="matplotlib", 
                       # max_spikes_per_unit=100,
                       templates_percentile_shading=None,
                       # same_axis=True,
                       # plot_channels=True,
                       alpha_waveforms=0.25,
                       figsize=(4.5 * min(num_units, 5), 4 * (num_units // 5 + 1) ),
                       axis_equal=True,
                      )

widg.figure.suptitle("Unit Templates", y=1.03, fontsize='x-large')

# ax.set_xlim(2100, 1900)
# ax.set_ylim(-100, 100)

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_unittemp', bbox_inches='tight')
plt.show()


# In[ ]:


widg = sw.plot_unit_locations(we, unit_colors=unit_cmap, with_channel_ids=True, plot_legend=True, figsize=(6, 6))
widg.ax.set_title("Unit Locations", y=1.2)
widg.ax.set_xlim(-40, 40)
widg.ax.set_ylim(-40, 40)

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_unitloc', bbox_inches='tight')
plt.show()


# In[ ]:


widg = sw.plot_unit_probe_map(we, figsize=(5 * num_units, 4), animated=True)

# plt.subplots_adjust(wspace=0.5)
# widg.animation.to_html5_video()
# widg.animation.to_html5_video()
if SAVEFIGS:
    widg.animation.save(output_folder_path / f'{tetrode_name}_animation.gif', fps=5, )

# plt.show()

# if SAVEFIGS:


# In[ ]:


widg = sw.plot_unit_presence(sorting,
                      time_range=[time_plotsample, time_plotsample + time_plotsample_range],
                      smooth_sigma=120,
                      figsize=(12, 1 * num_units))

widg.ax.set_title("Unit Presence", y=1.2)

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_unitpresence', bbox_inches='tight')
plt.show()


# In[ ]:


for i in range(num_units):
    sw.plot_unit_summary(we, i+1, unit_colors=unit_cmap, figsize=(12, 5))
    
    if SAVEFIGS:
        plt.savefig(output_folder_path / f'{tetrode_name}_unitsummary{i+1}', bbox_inches='tight')
    plt.show()


# In[ ]:


widg = sw.plot_unit_waveforms_density_map(we, figsize=(10, 2.5 * num_units))
widg.figure.suptitle("Unit Waveform Density Map")

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_densitymap', bbox_inches='tight')
plt.show()


# In[ ]:


widg = sw.plot_spike_locations(we, with_channel_ids=True, unit_colors=unit_cmap, plot_legend=True, figsize=(6, 6))
widg.ax.set_title("Spike Locations", y=1.2)
# widg.legend.set_loc('upper right')
# widg.legend.set_bbox_to_anchor((0.5, 0.9))
# print(widg.legend.get_bbox_to_anchor())
widg.ax.set_xlim((-40, 40))
widg.ax.set_ylim((-40, 40))

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_spikelocs2', bbox_inches='tight')
plt.show()


# In[ ]:


# Cross-Correlograms and Autocorrelograms
widg = sw.plot_crosscorrelograms(we, unit_colors=unit_cmap, figsize=(3 * num_units, 2 * num_units))
if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_ccg', bbox_inches='tight')
plt.show()


# In[ ]:


widg = sw.plot_isi_distribution(sorting, figsize=(4 * num_units, 3))
for ax in widg.axes[-1]:
    ax.set_xlabel("Time (ms)")
plt.subplots_adjust(wspace = 0.3)

if SAVEFIGS:
    plt.savefig(output_folder_path / f'{tetrode_name}_isi', bbox_inches='tight')
plt.show()

