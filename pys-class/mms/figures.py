#!/usr/bin/env python
# coding: utf-8

# # Make Stats And Figs

# A dedicated notebook for generating statistics and figures from sorting data.
# 
# Run after executing core modules. Will pull information from cached sorting and recording files in `sortings/`

# In[141]:


import glob
import shutil
import os
import csv
import gzip
from pathlib import Path
import datetime
from warnings import warn
from abc import ABC, abstractmethod
import re
import dateutil.parser as dparser
from textwrap import wrap
from functools import reduce
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# import probeinterface as pi
import pandas as pd
# import openpyxl
from statannotations.Annotator import Annotator
# import starbars
# from probeinterface.plotting import plot_probe
# import tempfile
# from tempfile import TemporaryDirectory
from scipy.optimize import curve_fit
from scipy.stats import skewnorm, norm, zscore
from scipy.signal import savgol_filter, savgol_coeffs, convolve, windows
import scipy.cluster.hierarchy as shc
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import euclidean_distances, manhattan_distances

import mountainsort5 as ms5
from mountainsort5.util import create_cached_recording
import spikeinterface.core as si
import spikeinterface.extractors as se
# import spikeinterface.preprocessing as spre
# import spikeinterface.sorters as ss
import spikeinterface.postprocessing as spost
import spikeinterface.qualitymetrics as sqm
# import spikeinterface.exporters as sexp
# import spikeinterface.comparison as scmp
# import spikeinterface.curation as scur
# import spikeinterface.sortingcomponents as sc
import spikeinterface.widgets as sw

from mms import core
# from mms.constants import GENOTYPE_ALIASES, REGIONS, FAST_JOB_KWARGS
from mms import constants
from mms.parser import FolderPathParser


# In[142]:


class DepthSheetReader():
    # REVIEW this can kind of just be a static class with states
    # But it is used as an instance class

    df_depths:pd.DataFrame = None
    df_optdepth:pd.DataFrame = None

    def __init__(self, base_folder: str, sortdir_name: str = 'sortings', truncate: bool = False, verbose: bool = True, omit: list[str] = []) -> None:
        self.base_folder = Path(base_folder)
        self.sortdir_name = sortdir_name
        self.sort_folder = self.base_folder / sortdir_name
        self.truncate = truncate
        self.verbose = verbose
        self.omit = omit

        self.sort_folders = self._glob_all_folders(self.sortdir_name)
        self.sort_all_ids = self._get_all_ids()

    # Theres some overlap in function implementation with core.IAnimalAnalyer,
    # but functionality > maintainability
    def _glob_all_folders(self, base_subdir_name:str):
        searchstring = self.base_folder / base_subdir_name / '*'
        subfolders = glob.glob(str(searchstring))
        subfolders = list(set(subfolders))
        subfolders = [x for x in subfolders if os.path.isdir(x)]
        subfolders = [x for x in subfolders if not any(substring in x for substring in self.omit)]
        subfolders.sort()
        return subfolders

    def make_filenametodepth_xlsx(self, xlsx_name:str='filename-to-depth.xlsx', xlsx_dir_name:str='xlsxs', read_old:bool=True):
        # self.sort_folders = self._glob_all_folders(self.sortdir_name)
        names = [FolderPathParser.parse_name(f) for f in self.sort_folders]
        ids = [FolderPathParser.parse_id(f, as_int=True) for f in self.sort_folders]
        genotypes = [FolderPathParser.parse_genotype(f) for f in self.sort_folders]
        dates = [FolderPathParser.parse_date(f) for f in self.sort_folders]
        df = pd.DataFrame({'folder' : self.sort_folders,
                           'name' : names,
                           'id' : ids,
                           'genotype' : genotypes,
                           'date' : dates,
                           'depth' : np.nan * len(self.sort_folders)})
        # df['id'] = df['id'].astype(str)
        today = datetime.date.today().strftime(r'%m-%d-%Y')
        out_path = self.base_folder / xlsx_dir_name / xlsx_name

        if read_old:
            # Merge all old spreadsheets
            cols = df.columns.drop(['depth']).tolist() # Without dropping depth this would just be concat
            dfs_old = pd.read_excel(out_path, sheet_name=None, index_col=0)
            dfs_old = list(dfs_old.values())
            dfs = [df] + dfs_old
            df_merge = reduce(lambda left,right: pd.merge(left, right, on=cols, how='outer', suffixes=('', '_y')), 
                              dfs)
            # Pull last value updated in the spreadsheet for each recording
            mapping = {'depth':'depth', 'depth_y':'depth'}
            mapping.update({k:v for k,v in zip(cols, cols)})
            df_last:pd.DataFrame = df_merge.set_index('name').T.groupby(mapping).last()
            df_last = df_last.T.reset_index().fillna(value=np.NaN)
            df_out = df_last[df.columns]
            
            df = df_out
            
        with pd.ExcelWriter(out_path, engine='openpyxl', mode='a', if_sheet_exists='new') as writer:
            df.to_excel(writer, sheet_name=today)
    
    def read_filenametodepth_xlsx(self, xlsx_name:str='filename-to-depth.xlsx', xlsx_dir_name:str='xlsxs', sheet_name=-1):
        sheet_path = self.base_folder / xlsx_dir_name / xlsx_name
        if sheet_name is None:
            df = pd.read_excel(sheet_path, sheet_name=-1, index_col=0)
        else:
            df = pd.read_excel(sheet_path, sheet_name=sheet_name, index_col=0)
        DepthSheetReader.df_depths = df
        return df
    
    @staticmethod
    def get_row_from_depth_df(identifier:int, region:str, depth:float=None, name:str=None):
        # def _load_sortinganalyzer(self, ):
        if depth is None and name is None:
            raise TypeError("Cannot have both depth and name as None")
        mask = DepthSheetReader.df_depths['id'] == identifier
        if depth is not None and depth == depth:
            mask = mask & (DepthSheetReader.df_depths['depth'] == depth)
        if name is not None and name == name:
            mask = mask & (DepthSheetReader.df_depths['name'] == name)

        df_mask = DepthSheetReader.df_depths[mask]
        if len(df_mask.index) == 0:
            warn("No match found for criteria")
            return None
        elif len(df_mask.index) >= 2:
            raise ValueError("Too many matches found for criteria")

        return df_mask.iloc[0]

    @staticmethod
    def __validate_and_get_dfdepth() -> pd.DataFrame:
        if DepthSheetReader.df_depths is None:
            raise TypeError('df_depths was not generated yet')
        return DepthSheetReader.df_depths

    @staticmethod
    def parse_depth(path:str):
        df = DepthSheetReader.__validate_and_get_dfdepth()
        name = FolderPathParser.parse_name(path)
        value = df.loc[df['name'] == name, 'depth']
        if len(value.index) == 0:
            raise IndexError('Path not found: {path}')
        elif len(value.index) >= 2:
            raise ValueError('Too many depths found for path: {path}')
        return value.iloc[0]

    def plot_all_depths_regions(self, identifier:int, output_subdir:str, sort_by_depth=True, save_figs=True, 
                                location_method: Literal['monopolar_triangulation', 'center_of_mass'] = 'monopolar_triangulation', **kwargs):
        df = DepthSheetReader.__validate_and_get_dfdepth()
        df = df[df['id'] == identifier]
        if sort_by_depth:
            if len(df.index) == 0:
                raise KeyError(f'ID {identifier} not found in depth short spreadsheet. sort_by_depth=True')
            df = df.sort_values(['depth', 'date', 'name'], ascending=[True, True, True])
            order = df['name'].tolist()
        else:
            order = None
        
        asl = core.AnimalSortLoader(self.base_folder, identifier=identifier, sortdir_name=self.sortdir_name)    
        
        for region in constants.REGIONS:
            asl.load_sortings_df()
            asl.extract_analyzers_df()

            asl.compute_extension_df(extension='random_spikes', region=region, max_spikes_per_unit=1000, **kwargs)
            asl.compute_extension_df(extension='waveforms', region=region, ms_before=1.0, ms_after=2.0, **kwargs)
            asl.compute_extension_df(extension='templates', region=region, ms_before=1.0, ms_after=2.0, **kwargs)
            asl.compute_extension_df(extension='spike_locations', region=region, ms_before=1.0, ms_after=2.0, method=location_method, **kwargs)
            # asl.compute_extension_df(extension='spike_locations', region=region, ms_before=0.5, ms_after=0.5, method='center_of_mass')

            tempsavepath = asl.base_folder / Path(output_subdir) / str(asl.identifier) / region
            if sort_by_depth:
                titles = df.apply(self.__get_title_from_row, axis=1, region=region).tolist()
            else:
                titles = None
            plot_params = {'region' : region,
                           'max_cols' : 10,
                           'order' : order,
                           'titles' : titles}
            asl.plot_units_oneregion(xybound_location=(-30, 30), plot_type='location', save_path=tempsavepath / 'location' if save_figs else None, **plot_params)
            asl.plot_units_oneregion(ybound_heatmap=(-200, 100), plot_type='heatmap', save_path=tempsavepath / 'heatmap' if save_figs else None, **plot_params)
            asl.plot_units_oneregion(ybound_wavetemp=(-40, 40), plot_type='waveform', save_path=tempsavepath / 'waveform' if save_figs else None, **plot_params)
            asl.plot_units_oneregion(ybound_wavetemp=(-40, 40), plot_type='template', save_path=tempsavepath / 'template' if save_figs else None, **plot_params)

            
    def __get_title_from_row(self, x:pd.Series, region:str):
        name = x['name']
        name = '\n'.join(wrap(name, 40))
        return f"{name}\n{x.id} | {x.genotype} | {x.depth} turns | {x.date.strftime(r'%m/%d/%Y')} | {region}"

    def _get_all_ids(self):
        ids = [FolderPathParser.parse_id(f, as_int=True) for f in self.sort_folders]
        ids = list(set(ids))
        ids.sort()
        return ids

    def make_animaltobestdepth_xlsx(self, xlsx_name:str='animal-to-bestdepth.xlsx', xlsx_dir_name:str='xlsxs', read_old:bool=True):
        
        sort_regs = sorted(constants.REGIONS)
        df = pd.DataFrame(
            {'id' : self.sort_all_ids,
            'region' : [sort_regs] * len(self.sort_all_ids),
            'best_depth' : [None] * len(self.sort_all_ids),
            'name' : [None] * len(self.sort_all_ids)}
        )
        df = df.explode(column='region')
        df = df.reset_index(drop=True)
        today = datetime.date.today().strftime(r'%m-%d-%Y')
        out_path = self.base_folder / xlsx_dir_name / xlsx_name

        if read_old:
            # Merge all old spreadsheets
            cols = df.columns.drop(['best_depth', 'name']).tolist()
            dfs_old = pd.read_excel(out_path, sheet_name=None, index_col=0)
            dfs_old = list(dfs_old.values())
            dfs = [df] + dfs_old
            df_merge = reduce(lambda left,right: pd.merge(left, right, on=cols, how='outer', suffixes=('', '_y')), 
                              dfs)

            # Pull last value updated in the spreadsheet for each recording
            mapping = {'best_depth':'best_depth', 'best_depth_y':'best_depth',
                       'name':'name', 'name_y':'name'}
            mapping.update({k:v for k,v in zip(cols, cols)})
            df_last:pd.DataFrame = df_merge.set_index(['id', 'region']).T.groupby(mapping).last()
            df_last = df_last.T.reset_index().fillna(value=np.NaN)
            df_out = df_last[df.columns]
            
            df = df_out

        with pd.ExcelWriter(out_path, engine='openpyxl', mode='a', if_sheet_exists='new') as writer:
            df.to_excel(writer, sheet_name=today)
    

    def read_animaltobestdepth_xlsx(self, xlsx_name:str='animal-to-bestdepth.xlsx', xlsx_dir_name:str='xlsxs', sheet_name=-1, drop_na:bool=True):
        sheet_path = self.base_folder / xlsx_dir_name / xlsx_name
        if sheet_name is None:
            df = pd.read_excel(sheet_path, sheet_name=-1, index_col=0)
        else:
            df = pd.read_excel(sheet_path, sheet_name=sheet_name, index_col=0)
            
        if drop_na:
            df = df.dropna(how='all', subset=['best_depth', 'name'])
            df = df.reset_index(drop=True)
        
        DepthSheetReader.df_optdepth = df
        return df
    


# In[144]:


class IExperimentAnalyzer(ABC):
    EXTENSION_PARAMS = {
        'random_spikes' : {'max_spikes_per_unit':1000},
        # 'waveforms' : dict({'ms_before':1.0, 'ms_after':2.0}, **constants.FAST_JOB_KWARGS),
        'waveforms' : {'ms_before':1.0, 'ms_after':2.0},
        'templates' : {'ms_before':1.0, 'ms_after':2.0},
        # 'spike_locations' : {'ms_before':1.0, 'ms_after':2.0, 'method':'center_of_mass'},
        # 'spike_locations' : {'ms_before':1.0, 'ms_after':2.0, 'method':'monopolar_triangulation'},
        # 'principal_components' : {'mode':'concatenated'},
        'principal_components' : {'mode':'by_channel_local'},
        # 'spike_amplitudes' : constants.FAST_JOB_KWARGS,
        'spike_amplitudes' : {},
        'noise_levels' : {},
    }
    # QUALITY_METRIC_FEATURES = sqm.get_quality_metric_list() + sqm.get_quality_pca_metric_list()
    # TEMPLATE_METRIC_FEATURES = spost.get_template_metric_names()
    QUALITY_METRIC_FEATURES = ['num_spikes', 'firing_rate', 'presence_ratio', 'snr',
       'isi_violations_ratio', 'isi_violations_count', 'rp_contamination',
       'rp_violations', 'sliding_rp_violation', 'amplitude_cutoff',
       'amplitude_median', 'amplitude_cv_median', 'amplitude_cv_range',
       'sync_spike_2', 'sync_spike_4', 'sync_spike_8', 'firing_range',
       'sd_ratio', 'isolation_distance', 'l_ratio', 'd_prime', 'silhouette',
       'nn_hit_rate', 'nn_miss_rate']
    TEMPLATE_METRIC_FEATURES = ['peak_to_valley', 'peak_trough_ratio',
       'half_width', 'repolarization_slope', 'recovery_slope',
       'num_positive_peaks', 'num_negative_peaks', 'velocity_above',
       'velocity_below', 'exp_decay', 'spread']
    DO_NOT_EXPLODE_FEATURES = ['id', 'region', 'best_depth', 'name', 'sa', 'n_units', 'waveform']
    UNIT_UNIQUE_FEATURES = ['id', 'region', 'best_depth', 'name', 'sa', 'n_units', 'unit', 'waveform']
    BAD_FEATURES = ['sliding_rp_violation', 'velocity_above', 'velocity_below',
                        'amplitude_cutoff', 'amplitude_cv_median', 'amplitude_cv_range']
    
    @abstractmethod
    def __init__(self):
        super().__init__()
        self.df_all_units: pd.DataFrame = None
    
    def get_df(self, drop_nan_rows=False, drop_nan_ignorecols:list[str]=[], 
               drop_bad_cols=True, drop_unit_unique_cols=False, drop_extra_cols:list[str]=None, drop_extra_cols_substr:list[str]=None,
               only_get_cols:list[str]=None):
        df = self.df_all_units.copy()
        if only_get_cols is not None:
            df = df.loc[:, only_get_cols]
        else:
            dropcols = []
            if drop_bad_cols:
                dropcols.extend(self.BAD_FEATURES)
            if drop_unit_unique_cols:
                dropcols.extend(self.UNIT_UNIQUE_FEATURES)
            if drop_extra_cols is not None:
                dropcols.extend(drop_extra_cols)
            if drop_extra_cols_substr is not None:
                matches = [[x for x in df.columns if substr in x] for substr in drop_extra_cols_substr]
                matches = [x for xs in matches for x in xs]
                # matches = np.array(matches)
                # matches = matches.flatten('C')
                dropcols.extend(matches)

            df = df.drop(dropcols, axis=1, errors='ignore')

        if drop_nan_rows:
            df = df.dropna(axis=0, subset=df.columns.difference(self.UNIT_UNIQUE_FEATURES + drop_nan_ignorecols))

        return df


class ExperimentFeatureExtractor(core.IAnimalAnalyzer, IExperimentAnalyzer):

    # FIXME sliding_rp_violation seems like a good stat, but returns many many np.Nans. See documentation to avoid this
    # FIXME i have a feeling these other stats are also good, but the spikeinterface settings need changing    

    def __init__(self, base_folder: str, dsr:DepthSheetReader, sortdir_name: str = 'sortings', 
                 truncate: bool = False, verbose: bool = True, omit: list[str] = ...,
                 **job_kwargs) -> None:
        
        super().__init__(base_folder, '', '', sortdir_name, truncate, verbose, omit)
        self._dsr = dsr
        self.__filename_to_depth_df = self._dsr.read_filenametodepth_xlsx()
        self.__animal_to_best_depth_df = self._dsr.read_animaltobestdepth_xlsx()

        self.df_all_units = self.__animal_to_best_depth_df.copy()
        if truncate:
            self.df_all_units = self.df_all_units.sample(5, random_state=42)

        self.df_all_units['sa'] = \
            self.df_all_units.apply(
                lambda row: self._load_sortinganalyzer(identifier=row.id,
                                                       region=row.region,
                                                       depth=row.best_depth,
                                                       name=row['name']), 
                axis=1
            )
        for k,kwargs in ExperimentFeatureExtractor.EXTENSION_PARAMS.items():
            self.df_all_units.apply(
                lambda row: row['sa'].compute_one_extension(extension_name=k, **(kwargs | job_kwargs)) if row['sa'] is not None else None,
                axis=1
            )
        self.df_all_units['n_units'] = self.df_all_units.apply(
            lambda row: row['sa'].get_num_units() if row['sa'] is not None else np.nan, axis=1)
        self.df_all_units['unit'] = self.df_all_units.apply(
            lambda row: list(range(1, row['sa'].get_num_units() + 1)) if row['sa'] is not None else np.nan, axis=1)
        
    def _load_sortinganalyzer(self, identifier:int, region:str, depth:float=None, name:str=None):
        row = DepthSheetReader.get_row_from_depth_df(identifier=identifier, region=region, depth=depth, name=name)
        return core.AnimalSortLoader.load_sortinganalyzer(row.folder, region=region)
    
    def remove_no_sa_rows(self, in_place=True, df=None):
        df = self.df_all_units if df is None else df
        df = df.dropna(subset=['sa'])
        if in_place:
            self.df_all_units = df
        return df

    def compute_quality_metrics(self, in_place=True, df=None):
        df = self.df_all_units if df is None else df
        df = df.apply(self.__compute_qm_row, axis=1)
        if in_place:
            self.df_all_units = df
        return df

    def __compute_qm_row(self, row: pd.Series):
        sa:si.SortingAnalyzer = row['sa']
        qm_ext = sa.compute('quality_metrics', skip_pc_metrics=False, **constants.FAST_JOB_KWARGS)
        metrics = qm_ext.get_data()
        for col in metrics.columns:
            row[col] = metrics[col]
        return row
    
    def compute_template_metrics(self, in_place=True, df=None):
        df = self.df_all_units if df is None else df
        df = df.apply(self.__compute_tm_row, axis=1)
        if in_place:
            self.df_all_units = df
        return df

    def __compute_tm_row(self, row:pd.Series):
        sa:si.SortingAnalyzer = row['sa']
        tm_ext = sa.compute('template_metrics', 
                            include_multi_channel_metrics=True, 
                            **constants.FAST_JOB_KWARGS)
        metrics = tm_ext.get_data()
        for col in metrics.columns:
            row[col] = metrics[col]
        return row
    
    def compute_mds_wave_features(self, n_components=4, in_place=True, df=None):
        df = self.df_all_units if df is None else df
        df = df.apply(self.__compute_mds_wave_row, axis=1)
        
        # Apply MDS
        X_waves = df['waveform']
        X_waves = np.array(X_waves.tolist())
        dist_waves = euclidean_distances(X_waves)
        mds = MDS(n_components=n_components, dissimilarity='precomputed')
        X_transform = mds.fit_transform(dist_waves)

        X_dict = {f'mds_wave_{i}' : X_transform[:, i] for i in range(X_transform.shape[1])}
        df = df.assign(**X_dict)
        
        if in_place:
            self.df_all_units = df
        return df
    
    def __compute_mds_wave_row(self, row: pd.Series):
        sa: si.SortingAnalyzer = row['sa']
        if sa is None:
            row['waveform'] = None
        else:
            unit = row['unit']
            ext_tmp: si.AnalyzerExtension = sa.get_extension('templates')
            # Kind of a hacky solution, unfortunately SA doesn't have a get_template function
            data = ext_tmp.get_data()
            extremum_ch = si.get_template_extremum_channel(sa, outputs='index')[unit]
            wave = data[unit - 1, :, extremum_ch]
            wave = (wave - wave.mean()) / wave.std() # NOTE standardizing the waveform
            row['waveform'] = wave
        return row

    # NOTE run this after quality or template metrics. You need to set the output to df_all_units manually
    def explode_dataframe(self, in_place=True, df=None):
        # explodecols = ExperimentFeatureExtractor.TEMPLATE_METRIC_FEATURES + \
        #     ExperimentFeatureExtractor.QUALITY_METRIC_FEATURES + ['unit']
        df = self.df_all_units if df is None else df
        explodecols = df.columns.difference(self.DO_NOT_EXPLODE_FEATURES, sort=False)
        explodecols = [x for x in explodecols if x in df.columns]
        df = df.explode(explodecols)
        if in_place:
            self.df_all_units = df
        return df

    def convert_datatypes(self, in_place=True, df=None):
        df = self.df_all_units if df is None else df
        datacols = df.columns.difference(self.DO_NOT_EXPLODE_FEATURES, sort=False)
        replace_dict = {'<NA>':np.nan, 'NA':np.nan, pd.NA:np.nan, None:np.nan}
        replace_dict = {x : replace_dict for x in datacols}
        df = df.replace(replace_dict)
        df = df.convert_dtypes()
        if in_place:
            self.df_all_units = df
        return df
    
    def remove_outliers(self, mode: Literal['iqr', 'std'] = 'iqr', iqr_mult:float = 1.5, std_mult:float = 3, 
                        ignore_cols: list[str]=['num_spikes', 'isi_violations_count', 'rp_violations'], 
                        ignore_cols_substr: list[str]=['num_', 'mds_wave_'], 
                        plot_rejects:bool = True,
                        in_place=True, df: pd.DataFrame=None):
        
        # unit: 0 outliers
        # num_spikes: 2 outliers
        # firing_rate: 0 outliers
        # presence_ratio: 1 outliers
        # snr: 1 outliers # REVIEW should this be taken out?
        # isi_violations_ratio: 1 outliers
        # isi_violations_count: 1 outliers
        # rp_contamination: 1 outliers
        # rp_violations: 0 outliers
        # amplitude_median: 1 outliers
        # sync_spike_2: 0 outliers
        # sync_spike_4: 0 outliers
        # sync_spike_8: 0 outliers
        # firing_range: 0 outliers
        # sd_ratio: 0 outliers
        # isolation_distance: 0 outliers
        # l_ratio: 0 outliers
        # d_prime: 0 outliers
        # silhouette: 0 outliers
        # nn_hit_rate: 0 outliers
        # nn_miss_rate: 0 outliers
        # peak_to_valley: 0 outliers
        # peak_trough_ratio: 0 outliers
        # half_width: 1 outliers
        # repolarization_slope: 0 outliers
        # recovery_slope: 0 outliers
        # num_positive_peaks: 4 outliers # REVIEW should this be considered for outlier rejection?
        # num_negative_peaks: 0 outliers # REVIEW same here
        # exp_decay: 0 outliers
        # spread: 0 outliers
        # mds_wave_0: 0 outliers
        # mds_wave_1: 0 outliers
        # mds_wave_2: 0 outliers
        # mds_wave_3: 0 outliers

        df = self.df_all_units if df is None else df
        
        ignore_cols = [] if ignore_cols is None else ignore_cols
        # datacols = df.columns.difference(self.DO_NOT_EXPLODE_FEATURES + self.BAD_FEATURES + ignore_cols, sort=False)
        # print(datacols)
        datacols = self.get_df(drop_unit_unique_cols=True, drop_extra_cols=ignore_cols, drop_extra_cols_substr=ignore_cols_substr).columns

        match mode:
            case 'iqr':
                dfiqr = df[datacols].quantile([0.25, 0.75])
                dfiqr.loc['iqr'] = dfiqr.iloc[1] - dfiqr.iloc[0]
                dfiqr.loc['min_bound'] = dfiqr.iloc[0] - iqr_mult * dfiqr.loc['iqr']
                dfiqr.loc['max_bound'] = dfiqr.iloc[1] + iqr_mult * dfiqr.loc['iqr']
                dfrange = dfiqr
            case 'std':
                std = df[datacols].std()
                mean = df[datacols].mean()
                dfrange = pd.concat([std, mean], axis=1)
                dfrange.loc['min_bound'] = dfrange['mean'] - std_mult * dfrange['std']
                dfrange.loc['max_bound'] = dfrange['mean'] + std_mult * dfrange['std']
            case _:
                raise ValueError(f'Invalid mode {mode}')
            
        df['reject'] = False
        for col in datacols:
            min_val = dfrange.loc['min_bound'][col]
            max_val = dfrange.loc['max_bound'][col]
            
            prev_reject = df['reject']
            df['reject'] = df['reject'] | (df[col] < min_val) | (df[col] > max_val)
            n_new_reject = df['reject'].sum() - prev_reject.sum()
            print(f'{col}: {n_new_reject} outliers')

            # ANCHOR fixing this
            if plot_rejects and n_new_reject > 0:
                sas = df.loc[df['reject'] != prev_reject, ['sa', 'unit']]
                fig, ax = plt.subplots(1, n_new_reject, figsize=(n_new_reject, 1), squeeze=False)
                for i, (index, row) in enumerate(sas.iterrows()):
                    sw.plot_unit_waveforms_density_map(row['sa'], use_max_channel=True, unit_ids=[row['unit']], ax=ax[0, i])
                    ax[0, i].set_facecolor('black')
                    ax[0, i].set_ylabel('')
                    ax[0, i].set_ylabel("uV")
                    # ax[i].set_ybound(ybound_heatmap[0], ybound_heatmap[1])
                fig.suptitle(col, y=1.05)
            
        if in_place:
            self.df_all_units = df
        return df


    def compute_clusters(self, feature_substrs: list[str] = None, threshold=0.7, plot_tree=True, in_place=True, df=None):
        df = self.df_all_units if df is None else df
        
        features = df.columns.difference(self.BAD_FEATURES)
        if feature_substrs is not None:
            allfeats = []
            for f in feature_substrs:
                allfeats.extend([x for x in features if f in x])
            print(f"Features used for clustering: {allfeats}")
            features = allfeats + self.UNIT_UNIQUE_FEATURES + ['reject']
        df2 = df.drop(df.columns.difference(features), axis=1, errors='ignore')
        df2 = df2.dropna(axis=0, subset=df2.columns.difference(self.DO_NOT_EXPLODE_FEATURES))
        if 'reject' in df2.columns:
            df2 = df2[~df2['reject']]

        df3 = df2.drop(self.UNIT_UNIQUE_FEATURES, axis=1, errors='ignore')
        df3_std = StandardScaler().fit_transform(df3)

        clusters = shc.linkage(df3_std, method='ward', metric='euclidean', optimal_ordering=True)
        self.Z_linkage = clusters
        threshold = threshold * max(clusters[:,2])
        cluster_labels = shc.fcluster(Z=clusters, t=threshold, criterion='distance')

        df2 = df2.assign(cluster=cluster_labels)
        df2['cluster'] = pd.to_numeric(df2['cluster'])
        df2 = df2.drop(df2.columns.difference(list(set(self.UNIT_UNIQUE_FEATURES) - set(['waveform'])) + ['cluster']), axis=1)
        df = df.merge(df2, how='left', on=list(set(self.UNIT_UNIQUE_FEATURES) - set(['waveform'])))
        
        if plot_tree:
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))
            shc.dendrogram(Z=clusters, ax=ax, color_threshold=threshold)
            ax.axhline(threshold, c='black', ls='--')
            plt.show()

        if in_place:
            self.df_all_units = df
        return df



# In[145]:


class ExperimentPlotter(core.IAnimalAnalyzer, IExperimentAnalyzer):

    def __init__(self, base_folder:str, efe:ExperimentFeatureExtractor = None,
                 dsr:DepthSheetReader = None, sortdir_name:str = 'sortings', 
                 truncate: bool = False, verbose: bool = True):
        super().__init__(base_folder, '', '', sortdir_name, truncate, verbose, omit=None)
        self.df_all_units = efe.get_df(drop_bad_cols=False)
        self._dsr = dsr if dsr is not None else efe._dsr
        self.__filename_to_depth_df = self._dsr.read_filenametodepth_xlsx()
        self.__animal_to_best_depth_df = self._dsr.read_animaltobestdepth_xlsx()
        if truncate:
            self.df_all_units = self.df_all_units.sample(10)


    def __get_clusterdf_clusterids_colorvec(self, features:list[str]=None, random_state:int=None):
        df = self.get_df(drop_nan_rows=True, drop_unit_unique_cols=True, drop_extra_cols=['cluster'])
        if features is not None:
            df = df.drop(df.columns.difference(features), axis=1)
        if 'reject' in df.columns:
            df = df[~df['reject']]
            
        df_clusters = self.get_df(drop_nan_rows=True)
        cluster_ids = df_clusters["cluster"].to_numpy(int)
        colorvec = [f'C{int(x)}' for x in df_clusters['cluster']]
        return df, cluster_ids, colorvec

    def __aggregate_and_show_legend(self, ax):
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=[1.05, 0.5], loc='center left')

    def plot_pca(self, features:list[str]=None):
        df_pca, cluster_ids, colorvec = self.__get_clusterdf_clusterids_colorvec()

        df_pca = StandardScaler().fit_transform(df_pca)
        pca = PCA(n_components=2)
        X_transform = pca.fit_transform(df_pca)

        fig, ax = plt.subplots(1, 1)
        for i in range(X_transform.shape[0]):
            ax.scatter(X_transform[i, 0], X_transform[i, 1], c=colorvec[i], label=f'Cluster {cluster_ids[i]}')

        self.__aggregate_and_show_legend(ax)

        plt.title('PCA Plot')
        plt.show()

    def plot_mds(self, features:list[str]=None, random_state:int=None):
        df_mds, cluster_ids, colorvec = self.__get_clusterdf_clusterids_colorvec()

        df_mds = StandardScaler().fit_transform(df_mds)
        mds = MDS(n_components=2, random_state=random_state)
        X_transform = mds.fit_transform(df_mds)

        fig, ax = plt.subplots(1, 1)
        for i in range(X_transform.shape[0]):
            ax.scatter(X_transform[i, 0], X_transform[i, 1], c=colorvec[i], label=f'Cluster {cluster_ids[i]}')

        self.__aggregate_and_show_legend(ax)

        plt.title('MDS Plot')
        plt.show()

    def plot_mds_waveforms(self, features:list[str]=None, width=2, height_mult=0.1, random_state:int=None):
        df_mds, cluster_ids, colorvec = self.__get_clusterdf_clusterids_colorvec()

        waveforms = np.array(self.get_df(drop_nan_rows=True)['waveform'].tolist())
        waveforms *= height_mult
        xvals = np.linspace(-width/2, width/2, waveforms.shape[1])

        df_mds = StandardScaler().fit_transform(df_mds)
        mds = MDS(n_components=2, random_state=random_state)
        X_transform = mds.fit_transform(df_mds)

        fig, ax = plt.subplots(1, 1)
        for i in range(X_transform.shape[0]):
            ax.plot(xvals + X_transform[i, 0], waveforms[i, :] + X_transform[i, 1], c=colorvec[i], label=f'Cluster {cluster_ids[i]}')
            
        self.__aggregate_and_show_legend(ax)

        plt.title('MDS Plot with Waveforms')
        plt.show()

    def plot_cluster_waveforms(self):
        df_waves = self.get_df(drop_nan_rows=True, drop_nan_ignorecols=['cluster'])
        df_waves.groupby('cluster', dropna=False).apply(self.__df_cluster_groupby_func)

    def __df_cluster_groupby_func(self, df: pd.DataFrame):
        n = df.index.size
        n_col = 10
        n_row = n // n_col + 1
        fig, axes = plt.subplots(n_row, n_col, figsize=(n_col*2, n_row*2), sharex=True, sharey=True, squeeze=False)
        axes[0, 0].set_aspect('equal')
        axes[0, 0].set_ylim(-20, 20)
        for i, (index, row) in enumerate(df.iterrows()):
            ax = axes[i // n_col, i % n_col]
            sa: si.SortingAnalyzer = row['sa']
            sw.plot_unit_waveforms(sa,
                                    unit_ids=[row['unit']],
                                    templates_percentile_shading=None,
                                    alpha_waveforms=0.1,
                                    set_title=False,
                                    scale=1e16,
                                    same_axis=True,
                                    plot_legend=False,
                                    ax=ax)
            if row['reject']:
                ax.set_facecolor('xkcd:red orange')
            ax.set_title(f"{row['id']} {row['region']} #{row['unit']}")
            ax.text(0.05, 0.05, f"n={row['num_spikes']} / {round(row['num_spikes']/sa.get_total_duration(), 2)} Hz", 
                    transform=ax.transAxes, 
                    ha='left', va='bottom', c='black', fontsize='small', fontweight='bold')
        fig.suptitle(f"Cluster {df.name} waveforms", y=1.05, fontsize='xx-large')
        plt.show()


    def plot_feature_boxplots(self, by: Literal['cluster', 'region']='cluster', features:list[str]=None):
        
        df_box = self.get_df(drop_nan_rows=True, drop_unit_unique_cols=True)
        if features is not None:
            df_box = df_box.drop(df_box.columns.difference(features), axis=1)

        df_box.boxplot(by=by, **{'sharey':False, 'figsize':(20, 20)})
        plt.show()



    # def plot_region_boxplots(self):
    #     # Group by region, plot the values of each column
    #     pass
    


# ## Read in Optimal Depths

# ## Load Saved Recordings and Sortings

# ## Calculate Spike Statistics

# ## Calculate Unit Quality Metrics

# ## Remove Low Quality Units

# ## Visualize All High Quality Units

# ## Plot Spike Statistics
