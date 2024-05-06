.PHONY: help clean ddi

help:
	@echo ''
	@echo -e '######################## HELP ##########################'
	@echo ''
	@echo -e '>>> Usage: peligand [\033[1m\033[31moptions\033[0m]'
	@echo ''
	@echo -e '    \033[1m\033[31mddi\033[0m: generates an input for the performance of the docking.' 
	@echo ''
	@echo -e '    \033[1m\033[31mrank\033[0m: ranks the results of the docking simulations.' 
	@echo ''
	@echo -e '    \033[1m\033[31mmutate\033[0m: Performs n_mutations single aa mutations to the protein.' 
	@echo ''
	@echo -e '########################################################'
	@echo ''

ddi:
	@./codes/run_files/c_i.sh

rank:
	@./codes/run_files/rank.sh

mutate:
	@./codes/run_files/mutate.sh	

clean: 
	@rm -rf proteins/paired_structures/*
	@rm -rf results/*
	@rm -f inputs/*
