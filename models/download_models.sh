#!/bin/bash

for (( model=1; model<=9; model++)); do
    wget "http://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD000000000$model" 
	mkdir BIOMD000000000$model
	mv download?mid=BIOMD000000000$model BIOMD000000000$model/sbml.xml
  done

for (( model=10; model<=99; model++)); do
    wget "http://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD00000000$model" 
	mkdir BIOMD00000000$model
	mv download?mid=BIOMD00000000$model BIOMD00000000$model/sbml.xml
  done

for (( model=100; model<=706; model++)); do
    wget "http://www.ebi.ac.uk/biomodels-main/download?mid=BIOMD0000000$model" 
	mkdir BIOMD0000000$model
	mv download?mid=BIOMD0000000$model BIOMD0000000$model/sbml.xml
  done
