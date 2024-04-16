.PHONY: run help plot rank cdb 
MAKE_COMMAND:=diffduck
help:
	@echo ''
	@echo -e '######################## HELP ##########################'
	@echo ''
	@echo -e '>>> Usage: diffduck [options]'
	@echo ''
	@echo -e '   run: Provided that you have placed the analyzed' 
	@echo -e '\tresults of the DiffDock run (whether each sample'
	@echo -e '\tattached to chain 1 or 0, i.e. results.csv) in'
	@echo -e '\tthe DD/ directory, generates all the outputs and'
	@echo -e '\tthe plots'
	@echo -e '  plot: Plots the results after ranking the results.'
	@echo -e '  rank: Ranks the results after the DiffDock run.'
	@echo -e '   cdb: Cleans the database and generates the necessary'
	@echo -e '\tfiles for generating a DiffDock input.'
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