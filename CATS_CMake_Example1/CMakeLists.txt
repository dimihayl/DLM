CMAKE_MINIMUM_REQUIRED(VERSION 1.0)

project(COMPILATOR)
# SET PATHS #
SET(ROOT_INCLUDE "/home/dmihaylov/root/include")#where are all ROOT related .h files
SET(ROOT_LIB "/home/dmihaylov/root/lib")#where are the ROOT related .a and .so files
SET(GSL_INCLUDE "/usr/include/gsl")#where are all GSL related .h files
SET(GSL_LIB "/usr/lib")#where are the GSL .a and .so files
SET(CATS_INCLUDE "/home/dmihaylov/Temp/CATS_TEST1/include")#where are all CATS related .h files
SET(CATS_LIB "/home/dmihaylov/Temp/CATS_TEST1/lib")#where are the CATS related .a files

add_executable(EXECUTABLE main.cpp
DimiFit.cpp
#ADD BELOW ANY OTHER .CPP FILES THAT YOU WOULD LIKE TO COMPILE (FULL PATH)
)
# DO NOT CHANGE THE REST #

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 0)
SET(VERSION_PATCH 0)
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CATS_INCLUDE})
include_directories(${GSL_INCLUDE})
include_directories(${ROOT_INCLUDE})

target_link_libraries(EXECUTABLE -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev -L${ROOT_LIB} libAfterImage.a libASImage.so libASImageGui.so libCint.so libCintex.so libcomplexDict.so libCore.so libdequeDict.so libEG.so libFitPanel.so libFoam.so libFumili.so libGed.so libGenetic.so libGenVector.so libGeom.so libGeomBuilder.so libGeomPainter.so libGpad.so libGraf.so libGraf3d.so libGui.so libGuiBld.so libGuiHtml.so libGX11.so libGX11TTF.so libHist.so libHistPainter.so libHtml.so liblistDict.so liblzma.a libmap2Dict.so libmapDict.so libMathCore.so libmathtext.a libMatrix.so libMemStat.so libMinuit.so libMLP.so libmultimap2Dict.so libmultimapDict.so libmultisetDict.so libNet.so libNew.so libpcre.a libPhysics.so libPostscript.so libProof.so libProofBench.so libProofDraw.so libProofPlayer.so libQuadp.so libRecorder.so libReflex.so libReflexDict.so libRint.so libRIO.so libRootAuth.so libSessionViewer.so libsetDict.so libSmatrix.so libSpectrum.so libSpectrumPainter.so libSPlot.so libSQLIO.so libSrvAuth.so libThread.so libTMVA.so libTree.so libTreePlayer.so libTreeViewer.so libvalarrayDict.so libvectorDict.so libVMC.so libX3d.so libXMLIO.so -L${GSL_LIB} libgsl.a libgsl.so)

install (TARGETS EXECUTABLE DESTINATION ${PROJECT_SOURCE_DIR}/bin)

