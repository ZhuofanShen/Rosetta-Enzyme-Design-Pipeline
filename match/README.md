**Generate position files**
> cd CPG2

Make staples across chain A and chain B
> python generate_position_files.py CPG2_relaxed.pdb -chains A B -wl 2000

or make staples within chain A
> python generate_position_files.py CPG2_relaxed.pdb -chain A -wl 2000

> ls

CPG2-AB (or CPG2-A)
> ls CPG2-AB

CPG2-AB_1 CPG2-AB_2 CPG2-AB_3
> ls CPG2-A

CPG2-A_1 CPG2-A_2 CPG2-A_3 CPG2-A_4 CPG2-A_5
> cd ..

**Match**
> mkdir CPG2_pAaF

> cd CPG2_pAaF

Make staples across chain A and chain B
> ../scripts/match.sh -pos ../CPG2/CPG2-AB -linker ../pAaF/pAaF-intermediate-R -mem 8000

> ../scripts/match.sh -pos ../CPG2/CPG2-AB -linker ../pAaF/pAaF-intermediate-S -mem 8000

> ../scripts/match.sh -pos ../CPG2/CPG2-AB -linker ../pAaF/pAaF-product -mem 8000

or make staples within chain A
> ../scripts/match.sh -pos ../CPG2/CPG2-A -linker ../pAaF/pAaF-intermediate-R -mem 8000

> ../scripts/match.sh -pos ../CPG2/CPG2-A -linker ../pAaF/pAaF-intermediate-S -mem 8000

> ../scripts/match.sh -pos ../CPG2/CPG2-A -linker ../pAaF/pAaF-product -mem 8000

**Convert output CloudPDB files into FastDesign input files**

Install Anaconda3 and PyMOL
> conda config --set auto_activate_base false

> conda init

> eval "$(~/anaconda3/bin/conda shell.bash hook)"

Make sure that your Anaconda3 has the pymol module installed.
> source ~/anaconda3/bin/activate

Make staples across chain A and chain B
> ../scripts/generate_fast_design_input.sh -pos ../CPG2/CPG2-AB -dup true

or make staples within chain A
> ../scripts/generate_fast_design_input.sh -pos ../CPG2/CPG2-A

Unload Anaconda
> conda deactivate
