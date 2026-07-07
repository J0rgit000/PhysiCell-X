#Set of commands to run:
make cleanup-data
make clean
make reset
make prostate
python3 scripts/Create_Basic_Config_Files.py --cell_lines <cell_lines> #This will personalize the config files to the cell_line.
python3 scripts/physiboss_drugsim <project> <cell_line> <drugs> <drug_rest> <mode> <drug_concs> #This will create the final config files
#Example: python3 scripts/physiboss_drugsim.py -p prostate --cell_line LNCaP -d Ipatasertib,Afatinib,Ulixertinib --mode both --drug-concs IC10 IC90
bash scripts/Prepare_Simulations.sh <Num_Sim> <project>
#Example: bash scripts/Prepare_Simulations.sh 2 prostate
