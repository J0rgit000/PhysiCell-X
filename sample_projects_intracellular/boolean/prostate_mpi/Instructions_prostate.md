#Possible cell-lines: 22Rv1, BPH1, DU145, LNCaP, NCIH660, PC3, PWR1E, VCap
#Possible drugs: Ipatasertib, Afuresertib, Afatinib, Erlotinib, Ulixertinib, Luminespib, Trametinib, Selumetinib, Pictilisib, Alpelisib, BIBR1532
#Possible mode: single, double, both
#Set of commands to run:
make cleanup-data
make clean
make reset
make prostate
python3 scripts/Create_Basic_Config_Files.py --cell_lines <cell_lines> #This will personalize the config files to the cell_line.
#Example: python3 scripts/Create_Basic_Config_Files.py --cell_lines 22Rv1
python3 scripts/physiboss_drugsim.py <project> <cell_line> <drugs> <drug_rest> <mode> <drug_concs> #This will create the final config files
#Example: python3 scripts/physiboss_drugsim.py -p prostate --cell_line 22Rv1 -d Ipatasertib --mode single --drug-concs IC90
bash scripts/Prepare_Simulations.sh <Num_Sim> <project>
#Example: bash scripts/Prepare_Simulations.sh 1 prostate_mpi
