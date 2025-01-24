import subprocess

def main():

    defaultparams = "--basefolder /mnt/isilon/marsh_single_unit/MarshMountainSort/"
    batchparams = [
        # NOTE old sortings (Markus, Donald)
        "--id 766 --datadir pyeegbins",
        "--id 1131 --datadir pyeegbins",
        "--id 1132 --datadir pyeegbins",
        "--id 1133 --datadir pyeegbins",
        "--id 1155 --datadir pyeegbins",
        "--id 1158 --datadir pyeegbins",
        "--id 1176 --datadir pyeegbins",
        "--id 1177 --datadir pyeegbins",
        "--id 1178 --datadir pyeegbins",
        "--id 1185 --datadir pyeegbins",
        "--id 1203 --datadir pyeegbins",
        "--id 1211 --datadir pyeegbins",
        "--id 1214 --datadir pyeegbins",
        "--id 1225 --datadir pyeegbins",
        "--id 1226 --datadir pyeegbins",
        "--id 1227 --datadir pyeegbins",
        "--id 1231 --datadir pyeegbins",
        "--id 1233 --datadir pyeegbins",
        "--id 1234 --datadir pyeegbins",
        "--id 1236 --datadir pyeegbins",
        "--id 1238 --datadir pyeegbins",
        "--id 1244 --datadir pyeegbins",
        "--id 1269 --datadir pyeegbins",
        # NOTE new sortings (Joseph)
        "--id 501-502 --datadir rhds --intanport A -aI A B",
        "--id 501-502 --datadir rhds --intanport B -aI A B",
        "--id 513-514 --datadir rhds --intanport A -aI A B -o 14htL",
        "--id 513-514 --datadir rhds --intanport B -aI A B",
        "--id 529-530 --datadir rhds --intanport A -aI A B",
        "--id 529-530 --datadir rhds --intanport B -aI A B",
        "--id 537 --datadir rhds",
        "--id 556-557 --datadir rhds --intanport A -aI A B",
        "--id 556-557 --datadir rhds --intanport B -aI A B",
        "--id 579-580 --datadir rhds --intanport A -aI A B -o 'overnight 1-17-25'",
        "--id 579-580 --datadir rhds --intanport B -aI A B -o 'overnight 1-17-25'",
        
        # NOTE sortings to depthplots
        # "--mode depth --id 766 --datadir bins",
        # "--mode depth --id 1131 --datadir bins",
        # "--mode depth --id 1132 --datadir bins",
        # "--mode depth --id 1133", # REVIEW testing
        # "--mode depth --id 1155", # REVIEW testing
        # "--mode depth --id 1158", # REVIEW testing
        # "--mode depth --id 1176", # REVIEW testing
        # "--mode depth --id 1177 --datadir bins",
        # "--mode depth --id 1178 --datadir bins",
        # "--mode depth --id 1185 --datadir bins",
        # "--mode depth --id 1203 --datadir bins",
        # "--mode depth --id 1211 --datadir bins",
        # "--mode depth --id 1214 --datadir bins",
        # "--mode depth --id 1225 --datadir bins",
        # "--mode depth --id 1226 --datadir bins",
        # "--mode depth --id 1227 --datadir bins",
        # "--mode depth --id 1231 --datadir bins",
        # "--mode depth --id 1233 --datadir bins",
        # "--mode depth --id 1234 --datadir bins",
        # "--mode depth --id 1236 --datadir bins",
        # "--mode depth --id 1238 --datadir bins",
        # "--mode depth --id 1244 --datadir bins",
        # "--mode depth --id 1269 --datadir bins",
        # "--mode depth --id 501 --datadir bins",
        # "--mode depth --id 502 --datadir bins",
        # "--mode depth --id 513 --datadir bins",
        # "--mode depth --id 514 --datadir bins",
        # "--mode depth --id 529 --datadir bins",
        # "--mode depth --id 530 --datadir bins",
        # "--mode depth --id 537 --datadir bins",
        # "--mode depth --id 556 --datadir bins",
        # "--mode depth --id 557 --datadir bins",
        # "--mode depth --id 579",
        # "--mode depth --id 580",
    ]

    for bp in batchparams:
        params = ' '.join([defaultparams, bp])
        __run_main(params)


    print("sbatcher.py finished")

def __run_main(params:str):
    with open('/mnt/isilon/marsh_single_unit/MarshMountainSort/pys-class/batcher.sh', 'w') as bsh:
        bsh.write(f"""\
#!/bin/sh
module load Python/3.10.8-GCCcore-12.2.0.lua
cd /mnt/isilon/marsh_single_unit/MarshMountainSort
source .venv/bin/activate
cd /mnt/isilon/marsh_single_unit/MarshMountainSort/pys-class
python -u main.py {params}
echo "sbatcher.py subprocess finished"
""")

    cmd = [
        'cd', '/mnt/isilon/marsh_single_unit/MarshMountainSort/;',
        'sbatch', 
        '--mem', '50G', # NOTE if not sorting 513-514, 20G. Otherwise 50G
        '-t', '24:00:00', 
        '/mnt/isilon/marsh_single_unit/MarshMountainSort/pys-class/batcher.sh']
    cmd = ' '.join(cmd)
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    main()