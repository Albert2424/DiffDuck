.PHONY: help clean ddi

help:
	@echo ''
	@echo -e '######################## HELP ##########################'
	@echo ''
	@echo -e '>>> Usage: peligand [\033[1m\033[31moptions\033[0m]'
	@echo ''
	@echo -e '    \033[1m\033[31mddi\033[0m: generates an input for DD.' 

	@echo ''
	@echo -e '########################################################'
	@echo ''

ddi:
	@./codes/run_files/c_i.sh
clean: 
	@rm -rf proteins/paired_structures/*
	@rm -rf results/*
	@rm -f inputs/*
