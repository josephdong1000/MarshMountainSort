# import sys
from timeit import default_timer as timer
import argparse
from tqdm import tqdm
from functools import partialmethod

from mms import core, figures, constants

def main():
    args = __get_parsed_args()
    print("Arguments:")
    for k,v in vars(args).items():
        print(f"\t{k}\t: {v}")
    
    tqdm.__init__ = partialmethod(tqdm.__init__, disable=not args.showtqdm)
    core.set_tempdir('/scr1/users/dongjp')
    # constants.GLOBAL_F_S = 25000 # global sampling frequency can be set like so
    # constants.GLOBAL_BANDPASS = [400, 1900] # same with global bandpass

    start = timer()

    match args.mode:
        case 'sort':
            ao = core.AnimalSorter(
                args.basefolder,  
                identifier=args.id,
                datadir_name=args.datadir,
                sortdir_name=args.sortdir,
                truncate=args.truncate,
                verbose=args.verbose,
                omit=args.omit,
                intan_port=args.intanport, 
                all_intan_ports=args.allintanports
            )
            ao.delete_sorting_folders()
            ao.sort_all()
        case 'depth':
            dsr = figures.DepthSheetReader(
                args.basefolder,
                sortdir_name=args.sortdir,
                truncate=args.truncate,
                verbose=args.verbose,
                omit=args.omit
            )
            dsr.read_filenametodepth_xlsx()
            dsr.plot_all_depths_regions(int(args.id), 
                                        output_subdir='output-depthsort/New Anims 1-16-25',
                                        # output_subdir='output-depthsort/New Anims 12-2-24',
                                        # output_subdir='output-depthsort/New Anims oldsort 12-17-24',
                                        sort_by_depth=not args.nosort,
                                        location_method='monopolar_triangulation',
                                        # location_method='center_of_mass',
                                        **constants.FAST_JOB_KWARGS)

        case 'figure':
            pass
        case _:
            raise ValueError('Invalid mode type for analysis')


    end = timer()
    print(f"Time elapsed: {end - start} s")

def __get_parsed_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-B', '--basefolder', help="BaseFolder for analysis", required=True)
    parser.add_argument('-i', '--id', help="Animal identifier", required=True)
    parser.add_argument('-D', '--datadir', default='rhds', help="Directory/filetype of data, either bins or rhds", choices=['bins', 'pyeegbins', 'rhds'])
    parser.add_argument('-S', '--sortdir', default='sortings', help="Directory/filetype of sorting files")
    parser.add_argument('-M', '--mode', help="Select mode for analysis", choices=['sort', 'depth', 'figure'], default='sort')
    parser.add_argument('-t', '--truncate', default=False, help="Only process the first few files. For testing only, default False")
    parser.add_argument('-v', '--verbose', default=True, help="Verbosity of script, default True")
    parser.add_argument('-o', '--omit', default=[], nargs='+', help="Substrings of filenames to omit, default empty")
    parser.add_argument('-ns', '--nosort', help="If flag is present, do not sort/label depth plots by depth. Set if incomplete depth turn sheet", action=argparse.BooleanOptionalAction)
    parser.add_argument('-I', '--intanport', help="Intan Port to read if datadir=rhds")
    parser.add_argument('-aI', '--allintanports', nargs='+', help="All available Intan Ports if datadir=rhds")
    parser.add_argument('-q', '--showtqdm', default=False, help="Show tqdm progress bars, defualt False")

    return parser.parse_args()


if __name__ == '__main__':
    main()