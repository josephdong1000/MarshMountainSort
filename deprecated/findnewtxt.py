#!/usr/bin/env python
# coding: utf-8

# # Find New TXTs
# 
# This will find any TXT files in the raw_txts folder that have not been processed, print out their names, and move them into their individual folders

# In[11]:


# Create new folders and move all raw .TXTs into individual folders
PUTINDIRECTORY = True

# Show converted files. Set to false if there are too many converted/outputted files
SHOWCONVERTED = False

# Show files that have figures generated. Set to false if there are too many files.
SHOWOUTPUTTED = False

# Compress all present .BIN files
COMPRESSBINS = False

# Move output files into subdirectories. You can also do this manually
MOVEOUTPUT = False


# In[12]:


import glob
import shutil
import os
import gzip
import numpy as np
from pathlib import Path


# ## Define Constants

# In[13]:


base_folder_path = Path(f'/mnt/isilon/marsh_single_unit/MarshMountainSort')
raw_txts_folder = 'rawtxts'
binary_folder = 'bins'
output_folder = 'output'
tetrode_names = ['ca3', 'ca1s', 'ca1o']

txt_folder_path = base_folder_path / raw_txts_folder
bin_folder_path = base_folder_path / binary_folder
output_folder_path = base_folder_path / output_folder

print(txt_folder_path)
print(bin_folder_path)
print(output_folder_path)


# ## Find Unconverted .TXTs

# In[15]:


# Quickly view all .TXT files in the raw_txts folder

# Find all files
txt_file_list = glob.glob(f'{txt_folder_path}/*.txt')
txt_in_folder_list = glob.glob(f'{txt_folder_path}/*/*.txt')
bin_file_list = glob.glob(f'{bin_folder_path}/*/*.bin')
gzip_file_list = glob.glob(f'{bin_folder_path}/*/*.npy.gz')

txt_file_names = [Path(e).stem for e in txt_file_list]
txt_in_folder_names = [Path(e).parent.stem for e in txt_in_folder_list]
bin_file_names = [Path(e).stem for e in bin_file_list]
gzip_file_names = [Path(e).name.partition('.')[0] for e in gzip_file_list]

# Print out files
print('===UNCONVERTED TXT FILES===')

print('NOT IN FOLDER:')
for e in txt_file_names:
    if e not in bin_file_names + gzip_file_names:
        print(f'"{e}",')
    else:
        if SHOWCONVERTED:
            print(f'\t\t[CONVERTED] {e}')
print()
print('IN FOLDER:')
for e in txt_in_folder_names:
    if e not in bin_file_names + gzip_file_names:
        print(f'"{e}",')
    else:
        if SHOWCONVERTED:
            print(f'\t\t[CONVERTED] {e}')

# 1185_Exp_4 half turns + 3-8th of a full turn_after turning_1-27-17


# ## Move .TXTs into Folders

# In[7]:


# Put exposed TXT files into individual directories
if PUTINDIRECTORY:
    for e in txt_file_names:
        fpath = txt_folder_path / e / f'{e}.txt'
        opath = txt_folder_path / f'{e}.txt'
        # print(txt_folder_path / f'{e}.txt')
        
        os.makedirs(fpath.parent, exist_ok=True)
        shutil.move(opath, fpath)


# ## Compress .BINs

# In[8]:


if COMPRESSBINS:
    for i, e in enumerate(bin_file_list):
        fcomp_name = Path(e).parent / f"{Path(e).stem}.npy.gz"
        
        # If the compressed file exists, skip
        if fcomp_name.is_file():
            continue
        
        bin_arr = np.fromfile(e, dtype=np.float32)
        with gzip.GzipFile(fcomp_name, "w") as fcomp:
            np.save(file=fcomp, arr=bin_arr)
        print(f"Compressed: {fcomp_name.name}")


# ## Find Unoutputted .BINs

# In[9]:


# Find unoutputted binaries and print them out
bin_folder_list = glob.glob(f'{bin_folder_path}/*')
bin_folder_names = [Path(e).stem for e in bin_folder_list]

conv_folder_parent_list = glob.glob(f'{output_folder_path}/*')
conv_folder_parent_names = [Path(e).stem for e in conv_folder_parent_list]

print('===UNOUTPUTTED BINARY FILES===')
for e in bin_folder_names:
    # print(e)
    if e not in conv_folder_parent_names:
        print(f'"{e}", ')


# ## Move Outputs into Folders

# In[17]:


# Move all output files into respective folders if applicable.
# If an output folder exists and is not empty, do not move and report the conflict
# Otherwise, create folder and move outputs into folder

if MOVEOUTPUT:
    print('===OUTPUT FOLDERS===')
    
    for e in conv_folder_parent_names:
        print(e)
        for i, region in enumerate(tetrode_names):
            region_folder_path = Path(output_folder_path / e / region)
    
            region_files = glob.glob(f'{region_folder_path.parent}/{region}*.*')
            region_present = len(region_files) != 0
            # print(region_present)
            
            # print(region_folder_path)
            
            if os.path.exists(region_folder_path):
                if len(os.listdir(region_folder_path)) != 0 and region_present:
                    print(f'\tCONFLICT {i} : {region} folder not empty')
                    continue
            elif region_present:
                os.makedirs(region_folder_path)
    
            
            for file in region_files:
                new_file = region_folder_path / Path(file).name
                print(f'\tMOVED FILE : {file}')
                shutil.move(file, new_file)
                
            


# ## Find Unoutputted Tetrode .BINs

# In[18]:


# Report which tetrodes have not been processed, based on presence of /ca3 /ca1s /ca1o subfolders
conv_folder_list = glob.glob(f'{output_folder_path}/*/*')
conv_folder_names = [Path(e).parent.stem for e in conv_folder_list]

print('===UNPROCESSED TETRODE BINARY FILES===')
for e in bin_file_names:
    if e not in txt_file_names + txt_in_folder_names:
        continue
    
    if e not in conv_folder_names:
        print(e)
    else:
        hasall = True
        msg = ''
        for i, region in enumerate(tetrode_names):
            if f'{output_folder_path}/{e}/{region}' not in conv_folder_list and f'{output_folder_path}/{e}/{region}_NO_UNITS' not in conv_folder_list:
                msg += f'\t\tMISSING OUTPUT {i} : {region}\n'
                hasall = False
            else:
                if SHOWOUTPUTTED:
                    msg += f'\t\tOK {i} : {region}\n'
        if not hasall:
            msg = f'\t[INCOMPLETE] {e}\n' + msg
        else:
            if SHOWOUTPUTTED:
                msg = f'\t[COMPLETE] {e}\n' + msg
        
        if msg != '':
            print(msg)


# In[ ]:




