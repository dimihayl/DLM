#!/bin/bash

#Make sure you set the correct paths
PATH_TO_CATS='../CATS'
PATH_TO_DLMCPPTOOLS='../DLM_CppTools'
PATH_TO_GSL_INCLUDE='/usr/include/gsl'
PATH_TO_GSL_LIB='/usr/lib'

rm -rf ./bin/*
rm -rf ./obj/*

g++ -Wall -fexceptions -g -I${PATH_TO_GSL_INCLUDE} -I${PATH_TO_CATS} -I${PATH_TO_DLMCPPTOOLS} -c ${PATH_TO_CATS}/CATS.cpp -o obj/CATS.o
g++ -Wall -fexceptions -g -I${PATH_TO_GSL_INCLUDE} -I${PATH_TO_CATS} -I${PATH_TO_DLMCPPTOOLS} -c ${PATH_TO_CATS}/CATStools.cpp -o obj/CATStools.o
g++ -Wall -fexceptions -g -I${PATH_TO_GSL_INCLUDE} -I${PATH_TO_CATS} -I${PATH_TO_DLMCPPTOOLS} -c ${PATH_TO_DLMCPPTOOLS}/DLM_CppTools.cpp -o obj/DLM_CppTools.o
g++ -Wall -fexceptions -g -I${PATH_TO_GSL_INCLUDE} -I${PATH_TO_CATS} -I${PATH_TO_DLMCPPTOOLS} -c main.cpp -o obj/main.o
g++ -o bin/Executable obj/CATS.o obj/CATStools.o obj/DLM_CppTools.o obj/main.o ${PATH_TO_GSL_LIB}/libgsl.a ${PATH_TO_GSL_LIB}/libgsl.so
