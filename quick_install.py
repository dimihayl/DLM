#!/usr/bin/env python
import sys,os,subprocess,time,multiprocessing,string


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#def create_cmake(PATH_DLM,PATH_CATS,PATH_GSL_INC,PATH_GSL_LIB,PATH_ROOT):



def is_gsl_include_okay(PATH_GSL_INC):
    gsl_dawson=False
    gsl_coulomb=False
    gsl_bessel=False
    gsl_legendre=False
    gsl_gamma=False
    OS_OUT = os.popen('ls -1 -U '+PATH_GSL_INC)
    OS_READ = OS_OUT.read()
    OS_SPLIT = OS_READ.splitlines()
    for strng in OS_SPLIT:
        if strng=='gsl_sf_dawson.h':
            gsl_dawson=True
        if strng=='gsl_sf_coulomb.h':
            gsl_coulomb=True
        if strng=='gsl_sf_bessel.h':
            gsl_bessel=True
        if strng=='gsl_sf_legendre.h':
            gsl_legendre=True
        if strng=='gsl_sf_gamma.h':
            gsl_gamma=True
    return (gsl_dawson and gsl_coulomb and gsl_bessel and gsl_legendre and gsl_gamma)

def is_gsl_lib_okay(PATH_GSL_LIB):
    OS_OUT = os.popen('ls -1 -U '+PATH_GSL_LIB)
    OS_READ = OS_OUT.read()
    OS_SPLIT = OS_READ.splitlines()
    for strng in OS_SPLIT:
        if strng=='libgsl.a' or strng=='libgsl.so':
            return True
    return False

def get_gsl_lib(PATH_GSL_LIB):
    libname = ''
    if '/libgsl.a' in PATH_GSL_LIB:
        libname = '/libgsl.a'
    elif '/libgsl.so' in PATH_GSL_LIB:
        libname = '/libgsl.so'
    if libname=='':
        return 'libgsl:'
    PATH_GSL_LIB = PATH_GSL_LIB.replace(libname,'')
    PATH_GSL_LIB = PATH_GSL_LIB.replace('\n','')
    return PATH_GSL_LIB




def quick_install(type):
    print(bcolors.BOLD+'CATS installer'+bcolors.ENDC)
    print(' To work properly, this files has to be executed from the DLM repository!')
    print(' Type "abort" to terminate the script.')
    print(' Continue? (y/n) ', end = '')
    answer = ''
    while True:
        answer = input()
        if answer=='abort':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return
        elif answer.lower()=="y" or answer.lower()=="yes":
            answer='y'
        elif answer.lower()=="n" or answer.lower()=="no":
            answer='n'
        else:
            answer=''

        if answer=='y':
            break
        elif answer=='n':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return
        else:
            print(' yes/no (y/n) answers only: ', end = '')

    install_mode = 'basic'
    if len(type)==2:
        install_mode = type[1]
    elif len(type)>2:
        print(bcolors.FAIL+'Too many input arguments, only one required!'+bcolors.ENDC)
        return

    #basic installs CATS + its non-root extensions
    #root includes any extensions that include ROOT
    #dev the full repository is installed
    if install_mode!='basic' and install_mode!='root' and install_mode!='dev':
        print(bcolors.FAIL+'Supported installation types: basic, root or dev'+bcolors.ENDC)
        return

    print(' Installation mode: '+install_mode)

    install_lvl=0
    if install_mode=='root':
        install_lvl=1
    elif install_mode=='dev':
        install_lvl=2

    #root installation
    root_loaded = True;
    PATH_ROOT = "root:"
    if install_lvl>=1:
        OS_OUT = os.popen('whereis root')
        OS_READ = OS_OUT.read()
        OS_SPLIT = OS_READ.split(" ")
        PATH_ROOT = OS_SPLIT[-1]
        #if "root:" in PATH_ROOT:
            #print(bcolors.FAIL+' ROOT installation not found!'+bcolors.ENDC)
            #print(' Please "source" the "thisroot.sh" of your ROOT installation before proceeding!')
        while "root:" in PATH_ROOT:
            #print(bcolors.FAIL+' ROOT installation not found!'+bcolors.ENDC)
            OS_OUT = os.popen('whereis root')
            OS_READ = OS_OUT.read()
            OS_SPLIT = OS_READ.split(": ")
            PATH_ROOT = OS_SPLIT[-1]
            if "root:" in PATH_ROOT:
                root_loaded = False
                print(bcolors.FAIL+' ROOT installation not found!'+bcolors.ENDC)
                print(' Please provide the path to the "thisroot.sh" of you ROOT installation:')
                print(' PATH_ROOT = ', end = '')
                PATH_ROOT = input()
                if PATH_ROOT=='abort':
                    print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
                    return
                #OS_CMD = '. '+PATH_ROOT+'/thisroot.sh'
                #print(OS_CMD)
                #os.system(OS_CMD)
                if not os.path.exists(PATH_ROOT+'/thisroot.sh'):
                    PATH_ROOT = "root:"

        #remove the name of the executable from the path
        if root_loaded:
            PATH_ROOT = PATH_ROOT[0:-5]
            if PATH_ROOT[-1]=="/":
                PATH_ROOT = PATH_ROOT[0:-1]
            print(bcolors.OKGREEN+' ROOT location found: '+bcolors.ENDC+PATH_ROOT)
        else:
            print(bcolors.WARNING+' ROOT found, but is not loaded in the environment'+bcolors.ENDC)
            print('  Solution (1): "source $PATH_ROOT/thisroot.sh" before using CATS')
            print('  Solution (2): Add the above command in your ~/.bashrc')
    #ROOT
        #remove the /bin
        PATH_ROOT = PATH_ROOT[0:-4]

    #path to the includes of GSL
    PATH_GSL_INC = ""
    PATH_GSL_INC = "gsl:"
    OS_OUT = os.popen('whereis gsl')
    #OS_OUT = "gsl:"
    OS_READ = OS_OUT.read()
    OS_SPLIT = OS_READ.split(" ")
    for PATH_GSL_INC in OS_SPLIT:
        if "gsl:" in PATH_GSL_INC:
            continue
        if is_gsl_include_okay(PATH_GSL_INC):
            print(bcolors.OKGREEN+' GSL location found: '+bcolors.ENDC+PATH_GSL_INC)
            break
    suggest = True
    while "gsl:" in PATH_GSL_INC:
        print(bcolors.FAIL+' GNU-GSL not found!'+bcolors.ENDC)
        if suggest:
            print(' Suggestion: intall using "sudo apt-get install libgsl-dev"')
            print(' Type "abort" to terminate the installation.')
            print(' Alternatively, type the path to the GSL "include" folder, containing "gsl_*.h" files')
            suggest = False
        print(' ', end = '')
        PATH_GSL_INC = input()
        if PATH_GSL_INC=='abort':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return

        if is_gsl_include_okay(PATH_GSL_INC):
            print(bcolors.OKGREEN+' GSL location valid!'+bcolors.ENDC)
        else:
            PATH_GSL_INC = "gsl:"

    #path to the lib of GSL
    PATH_GSL_LIB = ""
    PATH_GSL_LIB = "libgsl:"
    OS_OUT = os.popen('whereis libgsl')
    #OS_OUT = "gsl:"
    OS_READ = OS_OUT.read()
    OS_SPLIT = OS_READ.split(" ")
    for PATH_GSL_LIB in OS_SPLIT:
        if "libgsl:" in PATH_GSL_LIB:
            continue
        PATH_GSL_LIB = get_gsl_lib(PATH_GSL_LIB)
        if PATH_GSL_LIB!='libgsl:':
            print(bcolors.OKGREEN+' LIBGSL location found: '+bcolors.ENDC+PATH_GSL_LIB)
            break
    suggest = True
    while "libgsl:" in PATH_GSL_LIB:
        print(bcolors.FAIL+' LIBGSL not found!'+bcolors.ENDC)
        if suggest:
            print(' Suggestion: intall using "sudo apt-get install libgsl-dev"')
            print(' Type "abort" to terminate the installation.')
            print(' Alternatively, type the path to the GSL "lib" folder, containing "libgsl.so" and/or "libgsl.a" files')
            suggest = False
        print(' ', end = '')
        PATH_GSL_LIB = input()
        if PATH_GSL_LIB=='abort':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return

        if is_gsl_lib_okay(PATH_GSL_LIB):
            print(bcolors.OKGREEN+' LIBGSL location valid!'+bcolors.ENDC)
        else:
            PATH_GSL_LIB = "libgsl:"

    omp_okay = True
    if install_lvl>=2:
        omp_okay = (subprocess.check_call("apt show libomp-dev", shell=True, stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT))==0
        #>/dev/null 2>/dev/null

    if not omp_okay:
        print(bcolors.FAIL+'OpenMP needed, but NOT found!'+bcolors.ENDC)
        print(' Suggestion: intall using "sudo apt-get install libomp-dev"')
        return

    #path to the DLM repository
    OS_OUT = os.popen('pwd')
    PATH_DLM = OS_OUT.read()
    PATH_DLM = PATH_DLM.replace('\n','')
    #path to the installation folder for CATS
    PATH_CATS = PATH_DLM+'/install'
    PATH_CATS_CMAKE = PATH_CATS+'/CMake'

    time.sleep(1)
    print('--------------------')
    print(bcolors.OKGREEN+'The automatic setup was successful!'+bcolors.ENDC)
    print(bcolors.WARNING+'Please verify the settings ones more!'+bcolors.ENDC)
    print('--------------------')
    time.sleep(1)
    print(' DLM repository: '+PATH_DLM)
    print(' CATS to be installed in: '+PATH_CATS)
    print(' Using the following software:')
    print('  ROOT: '+PATH_ROOT)
    print('  GSL (include): '+PATH_GSL_INC)
    print('  GSL (lib): '+PATH_GSL_LIB)
    print(' Are these setting OKAY? (y/n) ', end = '')
    answer = ''
    while True:
        answer = input()
        if answer=='abort':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return
        elif answer.lower()=="y" or answer.lower()=="yes":
            answer='y'
        elif answer.lower()=="n" or answer.lower()=="no":
            answer='n'
        else:
            answer=''

        if answer=='y':
            break
        elif answer=='n':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return
        else:
            print(' yes/no (y/n) answers only: ', end = '')

    print('--------------------')
    print(bcolors.BOLD+'Setup the installation folder'+bcolors.ENDC)

    IncDir = []
    IncDir.append('CATS')
    IncDir.append('CATS_Extentions')
    IncDir.append('DLM_CppTools')
    IncDir.append('DLM_MathTools')
    if install_lvl>=1:
        IncDir.append('DLM_RootTools')
        IncDir.append('DLM_FemtoTools')
    if install_lvl>=2:
        IncDir.append('CECA')
    IncFile={}
    for dir in IncDir:
        if dir=='CATS':
            IncFile[dir] = ['CATS', 'CATStools']
        if dir=='CATS_Extentions':
            IncFile[dir] = ['DLM_Ck', 'DLM_ResponseMatrix', 'DLM_CkDecomp', 'DLM_Source', 'DLM_StefanoPotentials', 'DLM_Potentials', 'DLM_CkModels', 'DLM_WfModel']
            if install_lvl>=1:
                IncFile[dir].extend(['DLM_CkDecomposition', 'DLM_Fitters'])
        if dir=='DLM_CppTools':
            IncFile[dir] = ['DLM_CppTools']
            if install_lvl>=2:
                IncFile[dir].extend(['DLM_OmpTools'])
        if dir=='DLM_MathTools':
            IncFile[dir] = ['DLM_Histo','DLM_Integration','DLM_Random','DLM_Bessel','DLM_MathFunctions','DLM_Unfold','DLM_Fit','DLM_RootFinder']
        if dir=='DLM_RootTools':
            if install_lvl>=1:
                IncFile[dir] = ['DLM_CRAB_PM','DLM_SubPads','DLM_DrawingTools','DLM_HistoAnalysis','DLM_RootWrapper','DLM_MultiFit','DLM_DecayMatrix','DLM_RootFit']
        if dir=='DLM_FemtoTools':
            if install_lvl>=1:
                IncFile[dir] = ['CommonAnaFunctions']
        if dir=='CECA':
            if install_lvl>=2:
                IncFile[dir] = ['CECA','TREPNI']


    # checking whether folder exists or not
    # deletes if empty, creates it from scratch
    if os.path.exists(PATH_CATS):
        # checking whether the folder is empty or not
        if len(os.listdir(PATH_CATS)) != 0:
            print(bcolors.WARNING+' The installation folder is not empty, it will be DELETED!'+bcolors.ENDC)
            print('Proceed? (y/n) ', end = '')
            answer = ''
            while True:
                answer = input()
                if answer=='abort':
                    print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
                    return
                elif answer.lower()=="y" or answer.lower()=="yes":
                    answer='y'
                elif answer.lower()=="n" or answer.lower()=="no":
                    answer='n'
                else:
                    answer=''

                if answer=='y':
                    break
                elif answer=='n':
                    print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
                    return
                else:
                    print(' yes/no (y/n) answers only: ', end = '')
        os.system('rm -rf '+PATH_CATS)
    os.mkdir(PATH_CATS)
    print(' Folder '+PATH_CATS+' created')
    os.mkdir(PATH_CATS_CMAKE)
    print(' Folder '+PATH_CATS_CMAKE+' created')

    cmakelists = open(PATH_CATS_CMAKE+'/CMakeLists.txt', 'w')
    #cmakelists.Write('')

    cmakelists.write('CMAKE_MINIMUM_REQUIRED(VERSION 2.8)\n')
    cmakelists.write('if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")\n')
    cmakelists.write(' set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")\n')
    cmakelists.write(' set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")\n')
    cmakelists.write('endif()\n')
    cmakelists.write('\n')
    cmakelists.write('SET(PROJECT_NAME "CATS3")\n')
    cmakelists.write('\n')
    cmakelists.write('project(${PROJECT_NAME})\n')
    cmakelists.write('# SET PATHS #\n')
    cmakelists.write('SET(DLM_REPO "'+PATH_DLM+'")\n')
    cmakelists.write('SET(CATS_INSTALL "'+PATH_CATS+'")\n')
    cmakelists.write('SET(GSL_INCLUDE "'+PATH_GSL_INC+'")\n')
    cmakelists.write('SET(GSL_LIB "'+PATH_GSL_LIB+'")\n')
    if install_lvl>=1:
        cmakelists.write('SET(ROOT_PATH "'+PATH_ROOT+'")\n')

    cmakelists.write('\n')
    cmakelists.write('SET(VERSION_MAJOR 3)\n')
    cmakelists.write('SET(VERSION "${VERSION_MAJOR}")\n')
    cmakelists.write('SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CATS_INSTALL})\n')
    cmakelists.write('\n')

    cmakelists.write('include_directories(${GSL_INCLUDE})\n')
    for strng in IncDir:
        cmakelists.write('include_directories(${DLM_REPO}/'+strng+')\n')
    if install_lvl>=1:
        cmakelists.write('include_directories(${ROOT_PATH}/include)\n')

    cmakelists.write('add_library(CATS SHARED   \n')
    for dir in IncDir:
        for file in IncFile[dir]:
            cmakelists.write('              ${DLM_REPO}/'+dir+'/'+file+'.cpp\n')
    cmakelists.write('              )\n')

    cmakelists.write('\n')
    if install_lvl>=1:
        cmakelists.write('execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --cflags" OUTPUT_VARIABLE CFLAGS)\n')
        cmakelists.write('execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --libs" OUTPUT_VARIABLE LIBS)\n')
        cmakelists.write('string(REGEX REPLACE "\n$" "" CFLAGS "${CFLAGS}")\n')
        cmakelists.write('string(REGEX REPLACE "\n$" "" LIBS "${LIBS}")\n')

    if install_lvl==0:
        cmakelists.write('set(CFLAGS " -O2 -std=c++11")\n')
    elif install_lvl==1:
        cmakelists.write('set(CFLAGS " -O2 -std=c++11 ${CFLAGS}")\n')
    elif install_lvl==2:
        cmakelists.write('set(CFLAGS " -O2 -lgomp -pthread -fopenmp -std=c++11 ${CFLAGS}")\n')

    if install_lvl>=2:
        cmakelists.write('set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")\n')
        cmakelists.write('set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")\n')
        cmakelists.write('set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")\n')


    cmakelists.write('set_target_properties( CATS PROPERTIES COMPILE_FLAGS ${CFLAGS} ARCHIVE_OUTPUT_DIRECTORY "${CATS_INSTALL}/lib" LIBRARY_OUTPUT_DIRECTORY "${CATS_INSTALL}/lib" RUNTIME_OUTPUT_DIRECTORY "${CATS_INSTALL}/bin")\n')
    if install_lvl==0:
        cmakelists.write('target_link_libraries(  CATS -L${GSL_LIB} -lgsl -lgslcblas)\n')
    else:
        cmakelists.write('target_link_libraries(  CATS -L${GSL_LIB} -lgsl -lgslcblas ${LIBS})\n')

    for dir in IncDir:
        for file in IncFile[dir]:
            cmakelists.write('install (FILES ${DLM_REPO}/'+dir+'/'+file+'.h DESTINATION ${CATS_INSTALL}/include)\n')
    cmakelists.write('install (PROGRAMS ${DLM_REPO}/bin/cats-config DESTINATION ${CATS_INSTALL}/bin)\n')

    cmakelists.write('\n')

    cmakelists.write('file(WRITE ${CATS_INSTALL}/bin/CMakeDLM.txt "${PROJECT_NAME}\\n")\n')
    cmakelists.write('file(APPEND ${CATS_INSTALL}/bin/CMakeDLM.txt "${CATS_INSTALL}\\n")\n')
    cmakelists.write('file(APPEND ${CATS_INSTALL}/bin/CMakeDLM.txt "${ROOT_PATH}\\n")\n')
    cmakelists.write('file(APPEND ${CATS_INSTALL}/bin/CMakeDLM.txt "${GSL_INCLUDE}\\n")\n')
    cmakelists.write('file(APPEND ${CATS_INSTALL}/bin/CMakeDLM.txt "${GSL_LIB}\\n")\n')
    cmakelists.write('\n')

    cmakelists.close()

    print()
    print(bcolors.OKGREEN+'READY TO INSTALL!'+bcolors.ENDC)

    #return

    print(' Proceed? (y/n) ', end = '')
    answer = ''
    while True:
        answer = input()
        if answer=='abort':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return
        elif answer.lower()=="y" or answer.lower()=="yes":
            answer='y'
        elif answer.lower()=="n" or answer.lower()=="no":
            answer='n'
        else:
            answer=''

        if answer=='y':
            break
        elif answer=='n':
            print(bcolors.FAIL+'Installation canceled'+bcolors.ENDC)
            return
        else:
            print(' yes/no (y/n) answers only: ', end = '')

    #os.system('cd install/CMake/')
    os.chdir('./install/CMake/')
    if os.system('cmake .'):
        print(bcolors.FAIL+'"cmake" failed'+bcolors.ENDC)
        return
    if os.system('make -j'+str(multiprocessing.cpu_count())):
        print(bcolors.FAIL+'"make" failed'+bcolors.ENDC)
        return
    if os.system('make install'):
        print(bcolors.FAIL+'"make install" failed'+bcolors.ENDC)
        return
    print(bcolors.OKGREEN+'CATS was successfully installed!'+bcolors.ENDC)


quick_install(sys.argv)
