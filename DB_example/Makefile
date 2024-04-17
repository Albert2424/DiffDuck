.PHONY: run help plot rank cdb 
MAKE_COMMAND:=diffduck
help:
	@echo ''
	@echo -e '######################## HELP ##########################'
	@echo ''
	@echo -e '>>> Usage: diffduck [options]'
	@echo ''
	@echo -e '   \033[1m\033[31mrun\033[0m: Provided that you have placed the analyzed' 
	@echo -e '\tresults of the DiffDock run (whether each sample'
	@echo -e '\tattached to chain 1 or 0, i.e. results.csv) in'
	@echo -e '\tthe DD/ directory, generates all the outputs and'
	@echo -e '\tthe plots'
	@echo -e '   \033[1m\033[31mcdb\033[0m: Cleans the database and generates the necessary'
	@echo -e '\tfiles for generating a DiffDock input.'
	@echo -e '   \033[1m\033[31mfile_prep\033[0m: Generates the files and PDB necessary to'
	@echo -e '\trun the docking program.'
	@echo -e '   \033[1m\033[31mprep_ligands\033[0m: Uses the docking program selected to'
	@echo -e '\tpredict the ligand posses for the protein complex given.'
	@echo -e '   \033[1m\033[31manalysis\033[0m: Analyzes the results of the docking program.'
	@echo -e '   \033[1m\033[31mrank\033[0m: Ranks the results after the DiffDock run.'
	@echo -e '   \033[1m\033[31mplot\033[0m: Plots the results after ranking the results.'
	@echo ''
	@echo -e '########################################################'
	@echo ''

run:
	@./codes/run_files/run.sh

plot:
	@./codes/run_files/plot.sh
	
rank:
	@./codes/run_files/rank.sh

cdb:
	@./codes/run_files/cdb.sh

file_prep:
	@./codes/run_files/file_prep.sh

pred_ligands:
	@./codes/run_files/pred_ligands.sh

analysis:
	@./codes/run_files/analysis.sh