#!/usr/bin/env python
# coding: utf-8

# # Parser - Class

# In[3]:


from pathlib import Path
# import datetime
# from warnings import warn
# from abc import ABC, abstractmethod
import re
import dateutil.parser as dparser

from mms import constants
# import pandas as pd


# In[4]:


class FolderPathParser(object):

    @staticmethod
    def parse_name(path:str):
        return Path(path).name

    @staticmethod
    def parse_id(path:str, as_int=False):
        name = FolderPathParser.parse_name(path)
        idnum = re.split('[^0-9]', name)[0]
        if not idnum:
            raise ValueError('No ID found in name')
        if as_int:
            idnum = int(idnum)
        return idnum

    # @staticmethod
    # def parse_depth(path:str):
    #     return DepthSheetReader.parse_depth(path)
    
    @staticmethod
    def parse_genotype(path:str):
        name = FolderPathParser.parse_name(path)
        for key in constants.GENOTYPE_ALIASES.keys():
            if key in name:
                return key
        # Now check aliases
        for i, (k,v) in enumerate(constants.GENOTYPE_ALIASES.items()):
            if any([alias in name for alias in v]):
                return k
        raise KeyError(f'No valid genotype found in {name}')

    @staticmethod
    def parse_date(path:str, pattern=r'\d{1,2}[-_]\d{1,2}[-_]\d{2}(?:\d{2})?'):
        name = FolderPathParser.parse_name(path)
        name = name.removeprefix(FolderPathParser.parse_id(path))
        matches = re.findall(pattern, string=name)
        if len(matches) >= 2:
            raise ValueError(f'Too many date-like strings found in name: {matches}\nPath: {path}')
        elif len(matches) == 0:
            raise ValueError(f'No date-like strings found in name: {name}\nPath: {path}')
        
        name = matches[0]
        # Parse fuzzy with monthfirst always
        return dparser.parse(name, dayfirst=False, yearfirst=False, fuzzy=True)

