import subprocess

def main():

    defaultparams = "--basefolder /mnt/isilon/marsh_single_unit/MarshMountainSort/"
    batchparams = [
        # NOTE old sortings (Markus, Donald)
        # "--id 766 --datadir pyeegbins",
        # "--id 1131 --datadir pyeegbins",
        # "--id 1132 --datadir pyeegbins",
        # "--id 1133 --datadir pyeegbins",
        # "--id 1155 --datadir pyeegbins",
        # "--id 1158 --datadir pyeegbins",
        # "--id 1176 --datadir pyeegbins",
        # "--id 1177 --datadir pyeegbins",
        # "--id 1178 --datadir pyeegbins",
        # "--id 1185 --datadir pyeegbins",
        # "--id 1203 --datadir pyeegbins",
        # "--id 1211 --datadir pyeegbins",
        # "--id 1214 --datadir pyeegbins",
        # "--id 1225 --datadir pyeegbins",
        # "--id 1226 --datadir pyeegbins",
        # "--id 1227 --datadir pyeegbins",
        # "--id 1231 --datadir pyeegbins",
        # "--id 1233 --datadir pyeegbins",
        # "--id 1234 --datadir pyeegbins",
        # "--id 1236 --datadir pyeegbins",
        # "--id 1238 --datadir pyeegbins",
        # "--id 1244 --datadir pyeegbins",
        # "--id 1269 --datadir pyeegbins -o '12 hour recording'",
        # NOTE new sortings (Joseph)
        # "--id 501-502 --datadir rhds --intanport A -aI A B",
        # "--id 501-502 --datadir rhds --intanport B -aI A B",
        # "--id 513-514 --datadir rhds --intanport A -aI A B -o 14htL",
        # "--id 513-514 --datadir rhds --intanport B -aI A B -o 14htL",
        # "--id 529-530 --datadir rhds --intanport A -aI A B",
        # "--id 529-530 --datadir rhds --intanport B -aI A B",
        # "--id 537 --datadir rhds",
        # "--id 556-557 --datadir rhds --intanport A -aI A B",
        # "--id 556-557 --datadir rhds --intanport B -aI A B",
        # "--id 579-580 --datadir rhds --intanport A -aI A B -o 'overnight 1-17-25'",
        # "--id 579-580 --datadir rhds --intanport B -aI A B -o 'overnight 1-17-25'",
        # "--id 619-613 --datadir rhds --intanport A -aI A B",
        # "--id 619-613 --datadir rhds --intanport B -aI A B",

        # NOTE convert old sortings (Markus, Donald)
        *[f"--mode convert --sortdir recordings-raw --id {id} --datadir pyeegbins" for id in [
            766, 1131, 1132, 1133, 1155, 1158, 1176, 1177, 1178, 1185,
            1203, 1211, 1214, 1225, 1226, 1227, 1231, 1233, 1234, 1236,
            1238, 1244
        ]],
        "--mode convert --sortdir recordings-raw --id 1269 --datadir pyeegbins -o '12 hour recording'",
        # # NOTE convert new sortings (Joseph)
        *[f"--mode convert --sortdir recordings-raw --id {id} --datadir rhds --intanport {port} -aI A B" for id in [
            "501-502", "513-514", "529-530", "556-557", "579-580", "619-613"
        ] for port in ["A", "B"]],
        "--mode convert --sortdir recordings-raw --id 537 --datadir rhds",



        # NOTE sortings to depthplots
        # "--mode depth --id 766",
        # "--mode depth --id 1131",
        # "--mode depth --id 1132",
        # "--mode depth --id 1133",
        # "--mode depth --id 1155",
        # "--mode depth --id 1158",
        # "--mode depth --id 1176",
        # "--mode depth --id 1177",
        # "--mode depth --id 1178",
        # "--mode depth --id 1185",
        # "--mode depth --id 1203",
        # "--mode depth --id 1211",
        # "--mode depth --id 1214",
        # "--mode depth --id 1225",
        # "--mode depth --id 1226",
        # "--mode depth --id 1227",
        # "--mode depth --id 1231",
        # "--mode depth --id 1233",
        # "--mode depth --id 1234",
        # "--mode depth --id 1236",
        # "--mode depth --id 1238",
        # "--mode depth --id 1244",
        # "--mode depth --id 1269", # TODO retry
        # "--mode depth --id 501",
        # "--mode depth --id 502",
        # "--mode depth --id 513",
        # "--mode depth --id 514",
        # "--mode depth --id 529",
        # "--mode depth --id 530",
        # "--mode depth --id 537",
        # "--mode depth --id 556",
        # "--mode depth --id 557",
        # "--mode depth --id 579",
        # "--mode depth --id 580",
        # "--mode depth --id 613",
        # "--mode depth --id 619",
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
        '--mem', '50G',
        '-c', '8',  # Request 8 CPUs per task
        '-t', '24:00:00',
        '/mnt/isilon/marsh_single_unit/MarshMountainSort/pys-class/batcher.sh']
    cmd = ' '.join(cmd)
    subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    main()