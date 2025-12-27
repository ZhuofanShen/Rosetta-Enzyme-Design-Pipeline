```bash
sbatch scripts-match/generate_uaa_position_files.yCD.pbs

sbatch scripts-match/generate_uaa_position_files.antibody.pbs


sbatch scripts-match/match.yCD.pbs -nchains 1 -lig ligands/O2beY-CYS -pos *_*_C* -mem 10000

sbatch scripts-match/match.antibody.pbs -nchains 2 -lig ligands/O2beY-CYS -pos *_*_C* -mem 10000

sbatch scripts-match/match.antibody.pbs -nchains 2 -lig ligands/O2beY-HIS -pos *_*_H* -mem 10000

sbatch scripts-match/match.antibody.pbs -nchains 2 -lig ligands/O2beY-ASP -pos *_*_D* -mem 10000

sbatch scripts-match/match.antibody.pbs -nchains 2 -lig ligands/O2beY-GLU -pos *_*_E* -mem 10000

sbatch scripts-match/match.antibody.pbs -nchains 2 -lig ligands/O2beY-LYS -pos *_*_K* -mem 10000


conda activate py311

sbatch scripts-match/generate_fast_design_input.yCD.pbs -nchains 1 -lig ligands/O2beY-CYS -pos *_*_C*

sbatch scripts-match/generate_fast_design_input.antibody.pbs -nchains 2 -lig ligands/O2beY-CYS -pos *_*_C*


python scripts-design/fast_design_2_mpnn.py outyCD --process_variants 1ysb-A_D132_Q10 1ysb-A_L21_R134 1ysb-A_S6_E126 1ysb-A_Y24_E145 1ysb-A_D132_D9 1ysb-A_V105_D9 1ysb-A_V130_A5 1ysb-A_F52_E152 1ysb-A_L138_L21 -o mpnn -od mpnn_30

python scripts-design/fast_design_2_mpnn.py outantibody -o csv --process_variants 7l6v-A-B_E520_C166


conda activate aifold

sbatch scripts-design/mpnn.pbs outyCD mpnn_30


python scripts-design/esm_scoring.py outyCD --process_variants 1ysb-A_D132_Q10 1ysb-A_L21_R134 1ysb-A_S6_E126 1ysb-A_Y24_E145 1ysb-A_D132_D9 1ysb-A_V105_D9 1ysb-A_V130_A5 1ysb-A_F52_E152 1ysb-A_L138_L21 -od mpnn_30 -esm_n 10


conda activate py311

python scripts-design/mpnn_2_fast_relax.py outyCD --process_variants 1ysb-A_D132_Q10 1ysb-A_L21_R134 1ysb-A_S6_E126 1ysb-A_Y24_E145 1ysb-A_D132_D9 1ysb-A_V105_D9 1ysb-A_V130_A5 1ysb-A_F52_E152 1ysb-A_L138_L21 -od mpnn_30 -esm_n 10 -s 1

python scripts-design/mpnn_2_fast_relax.py outyCD --process_variants 1ysb-A_D132_Q10 1ysb-A_L21_R134 1ysb-A_S6_E126 1ysb-A_Y24_E145 1ysb-A_D132_D9 1ysb-A_V105_D9 1ysb-A_V130_A5 1ysb-A_F52_E152 1ysb-A_L138_L21 -od mpnn_30 -s 2


sbatch scripts-design/relax.pbs outyCD mpnn_30 relax.sh

sbatch scripts-design/relax.pbs outyCD mpnn_30 relax_WT.sh
```
