#!/usr/bin/env python
# coding: utf-8

# # Constants

# In[2]:


from typing import Literal

GENOTYPE_ALIASES = {'WT' : ['Wt', 'wildtype', 'Wildtype', 'Bl6'],
                    'Exp' : ['EXP', 'exp', 'Expanded', 'GCG']}
REGIONS = ['ca3', 'ca3o', 'ca1s', 'ca1o']
REGION_TO_INTAN_CHANNEL = {'ca3o': [8, 9, 10, 11],
                            'ca1o': [12, 13, 14, 15],
                            'ca1s': [16, 17, 18, 19],
                            'ca3': [20, 21, 22, 23]}
REGION_TO_DATAWAVE_CHANNEL = {'ca3': [0, 1, 2, 3],
                              'ca1o': [4, 5, 6, 7],
                              'ca1s': [10, 11, 12, 13]}
GLOBAL_F_S = 25000
GLOBAL_BANDPASS = [400, 8064]
DEPTHPLOT_TYPES = Literal['heatmap', 'location', 'waveform', 'template']
FAST_JOB_KWARGS = dict(n_jobs=4, chunk_duration='5s', progress_bar=True)

