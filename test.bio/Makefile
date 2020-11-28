
SHELL := /bin/bash

all-tests := $(addsuffix .test, $(basename $(wildcard *.awk)))

.PHONY : test all %.test

bioawk := ../bioawk

all: test 
	@echo "Success, all tests passed."
	
test: $(all-tests)

fastx_%.test : fastx_%.awk fastx_%.in fastx_%.expect $(bioawk)
	@$(bioawk) -c fastx -f $< $(word 2, $?) 2>&1 | diff -q $(word 3, $?) - >/dev/null || \
	(echo "Test $@ failed" && exit 1)

gff_%.test : gff_%.awk gff_%.in gff_%.expect $(bioawk)
	@$(bioawk) -c gff -f $< $(word 2, $?) 2>&1 | diff -q $(word 3, $?) - >/dev/null || \
	(echo "Test $@ failed" && exit 1)

%.test : %.awk %.in %.expect $(bioawk)
	@$(bioawk) -f $< $(word 2, $?) | diff -q $(word 3, $?) - >/dev/null || \
	(echo "Test $@ failed" && exit 1)