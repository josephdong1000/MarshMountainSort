from mms import constants, core, figures

def main():
    dsr = figures.DepthSheetReader('/mnt/isilon/marsh_single_unit/MarshMountainSort')
    # dsr.make_filenametodepth_xlsx()
    dsr.make_animaltobestdepth_xlsx()
    print("Done")

if __name__ == '__main__':
    main()