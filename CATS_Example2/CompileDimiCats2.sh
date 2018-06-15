#!/bin/bash
rm -rf obj/
rm -rf ./bin/

PATH_ROOT="/home/dmihaylov/root"
PATH_GSL_INCLUDE="/usr/include/gsl"
PATH_GSL_LIB="/usr/lib"
PATH_HOME="/home/dmihaylov/Dudek_Ubuntu/DLM_GitHub/CATS_Example2"
PATH_DLM="/home/dmihaylov/Dudek_Ubuntu/DLM_GitHub"

mkdir ./obj
mkdir ./obj/Debug
mkdir ./bin
mkdir ./bin/Debug

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_HOME/CatsExample2.cpp -o obj/Debug/CatsExample2.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS/CATS.cpp -o obj/Debug/CATS.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS/CATStools.cpp -o obj/Debug/CATStools.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/DLM_CppTools/DLM_CppTools.cpp -o obj/Debug/DLM_CppTools.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/DLM_RootTools/DLM_DrawingTools.cpp -o obj/Debug/DLM_DrawingTools.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/DLM_RootTools/DLM_SubPads.cpp -o obj/Debug/DLM_SubPads.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS_Extentions/DLM_StefanoPotentials.cpp -o obj/Debug/DLM_StefanoPotentials.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS_Extentions/DLM_Potentials.cpp -o obj/Debug/DLM_Potentials.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS_Extentions/DLM_Source.cpp -o obj/Debug/DLM_Source.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS_Extentions/DLM_ResponseMatrix.cpp -o obj/Debug/DLM_ResponseMatrix.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS_Extentions/DLM_WfModel.cpp -o obj/Debug/DLM_WfModel.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS_Extentions/DLM_CkModels.cpp -o obj/Debug/DLM_CkModels.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS_Extentions/DLM_CkDecomposition.cpp -o obj/Debug/DLM_CkDecomposition.o

g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_DLM/CATS_Extentions/DLM_Fitters.cpp -o obj/Debug/DLM_Fitters.o

#to add further of your files to the compile list, find all occurrences of main copy-paste them, changing the name
g++ -Wall -fexceptions -g -I$PATH_ROOT/include -I$PATH_DLM/CATS -I$PATH_DLM/CATS_Extentions -I$PATH_DLM/DLM_CppTools -I$PATH_DLM/DLM_RootTools -I$PATH_GSL_INCLUDE -ICATSpotentials -I$PATH_HOME/ -c $PATH_HOME/main.cpp -o obj/Debug/main.o

g++ -L$PATH_ROOT/lib -L$PATH_GSL_LIB -o bin/Debug/Using_CATS obj/Debug/CatsExample2.o obj/Debug/CATS.o obj/Debug/CATStools.o obj/Debug/DLM_CppTools.o obj/Debug/DLM_DrawingTools.o obj/Debug/DLM_SubPads.o obj/Debug/DLM_StefanoPotentials.o obj/Debug/DLM_Potentials.o obj/Debug/DLM_Source.o obj/Debug/DLM_ResponseMatrix.o obj/Debug/DLM_WfModel.o obj/Debug/DLM_CkModels.o obj/Debug/DLM_CkDecomposition.o obj/Debug/DLM_Fitters.o obj/Debug/main.o $PATH_ROOT/lib/libASImage.so $PATH_ROOT/lib/libASImageGui.so $PATH_ROOT/lib/libCint.so $PATH_ROOT/lib/libCintex.so $PATH_ROOT/lib/libcomplexDict.so $PATH_ROOT/lib/libCore.so $PATH_ROOT/lib/libdequeDict.so $PATH_ROOT/lib/libEG.so $PATH_ROOT/lib/libFitPanel.so $PATH_ROOT/lib/libFoam.so $PATH_ROOT/lib/libFumili.so $PATH_ROOT/lib/libGed.so $PATH_ROOT/lib/libGenVector.so $PATH_ROOT/lib/libGeom.so $PATH_ROOT/lib/libGeomBuilder.so $PATH_ROOT/lib/libGeomPainter.so $PATH_ROOT/lib/libGpad.so $PATH_ROOT/lib/libGraf.so $PATH_ROOT/lib/libGraf3d.so $PATH_ROOT/lib/libGui.so $PATH_ROOT/lib/libGuiBld.so $PATH_ROOT/lib/libGuiHtml.so $PATH_ROOT/lib/libGX11.so $PATH_ROOT/lib/libGX11TTF.so $PATH_ROOT/lib/libHist.so $PATH_ROOT/lib/libHistPainter.so $PATH_ROOT/lib/libHtml.so $PATH_ROOT/lib/liblistDict.so $PATH_ROOT/lib/libmap2Dict.so $PATH_ROOT/lib/libmapDict.so $PATH_ROOT/lib/libMathCore.so $PATH_ROOT/lib/libmathtext.a $PATH_ROOT/lib/libMatrix.so $PATH_ROOT/lib/libMemStat.so $PATH_ROOT/lib/libMinuit.so $PATH_ROOT/lib/libMLP.so $PATH_ROOT/lib/libmultimap2Dict.so $PATH_ROOT/lib/libmultimapDict.so $PATH_ROOT/lib/libmultisetDict.so $PATH_ROOT/lib/libNet.so $PATH_ROOT/lib/libNew.so $PATH_ROOT/lib/libPhysics.so $PATH_ROOT/lib/libPostscript.so $PATH_ROOT/lib/libProof.so $PATH_ROOT/lib/libProofBench.so $PATH_ROOT/lib/libProofDraw.so $PATH_ROOT/lib/libProofPlayer.so $PATH_ROOT/lib/libQuadp.so $PATH_ROOT/lib/libRecorder.so $PATH_ROOT/lib/libReflex.so $PATH_ROOT/lib/libReflexDict.so $PATH_ROOT/lib/libRint.so $PATH_ROOT/lib/libRIO.so $PATH_ROOT/lib/libRootAuth.so $PATH_ROOT/lib/libSessionViewer.so $PATH_ROOT/lib/libsetDict.so $PATH_ROOT/lib/libSmatrix.so $PATH_ROOT/lib/libSpectrum.so $PATH_ROOT/lib/libSpectrumPainter.so $PATH_ROOT/lib/libSPlot.so $PATH_ROOT/lib/libSQLIO.so $PATH_ROOT/lib/libThread.so $PATH_ROOT/lib/libTree.so $PATH_ROOT/lib/libTreePlayer.so $PATH_ROOT/lib/libTreeViewer.so $PATH_ROOT/lib/libvalarrayDict.so $PATH_ROOT/lib/libvectorDict.so $PATH_ROOT/lib/libVMC.so $PATH_ROOT/lib/libX3d.so $PATH_ROOT/lib/libXMLIO.so $PATH_GSL_LIB/libgsl.a $PATH_GSL_LIB/libgsl.so $PATH_GSL_LIB/libgslcblas.a $PATH_GSL_LIB/libgslcblas.so

