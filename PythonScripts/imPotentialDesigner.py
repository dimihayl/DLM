import os,sys
import re
import ROOT
import math
import optuna
import subprocess
import time
import shutil

#!! this executable needs to be compiled. It should read the files produced by this file, and spit out a ROOT file
# containing the result for the scattering length given specific settings
script_name = "/Users/sartozza/Software/cats/CatsSource/bin/imPotentialDesigner"

FmToNu = 5.067731237e-3

# studyR = optuna.create_study(direction='minimize', sampler=optuna.samplers.RandomSampler())
studyT = optuna.create_study(direction='minimize', sampler=optuna.samplers.TPESampler())

# these are the ranges of the potential parameters, which we would like to study
par1_range = [0, 0]
par2_range = [0, 0]
par3_range = [0, 0]
par4_range = [0, 0]
par5_range = [0, 0]
par6_range = [0, 0]

# Target function to be optimized by Optuna
def target_function(trial):
    # Sample the parameters within the specified ranges
    v_par1 = trial.suggest_float('v_par1', *par1_range)
    v_par2 = trial.suggest_float('v_par2', *par2_range)
    v_par3 = trial.suggest_float('v_par3', *par3_range)
    v_par4 = trial.suggest_float('v_par4', *par4_range)
    v_par5 = trial.suggest_float('v_par5', *par5_range)
    v_par6 = trial.suggest_float('v_par6', *par6_range)
    return 0

def run_script(args):
    subprocess.Popen([script_name] + [args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def parse_float_range(range_str):
    return [float(value) for value in range_str.split()]

def get_yes_no_input(prompt):
    while True:
        user_input = input(prompt).strip().lower()
        if user_input in ["yes", "y"]:
            return True
        elif user_input in ["no", "n"]:
            return False
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")

def main():
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        return

    BaseFileName = sys.argv[1]
    filename_arg = BaseFileName+'.txt'
    # ROOT_FILE_NAME = sys.argv[1]+'.root'

    fRe_goal = 0
    fRe_err = 0
    fIm_goal = 0
    fIm_err = 0
    # d0Re_goal = 0
    # d0Re_err = 0
    pot = ""
    M1 = 0
    M2 = 0
    PW = 0
    q1q2 = 0
    kMin = 0
    kMax = 100
    kMaxFit = 100
    kBin = 50
    CPU = 1
    FIT_TIMEOUT = 1#in minutes
    TIMEOUT = 20#in minutes
    WakeUp = 2.5#in seconds

    if not os.path.exists(filename_arg):
        sys.exit(f"The file '{filename_arg}' does not exist.")

    # Open the file for reading
    with open(filename_arg, 'r') as file:
        for line in file:
            single_line = line.strip()
            line_split = single_line.split()
            # print(line_split)
            if line_split[0].lower()=='fre_goal':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of fRe_goal')
                try:
                    fRe_goal = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: fRe_goal must be a number.")
            if line_split[0].lower()=='fre_err':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of fRe_err')
                try:
                    fRe_err = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: fRe_err must be a number.")
            if line_split[0].lower()=='fim_goal':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of fIm_goal')
                try:
                    fIm_goal = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: fIm_goal must be a number.")
            if line_split[0].lower()=='fim_err':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of fIm_err')
                try:
                    fIm_err = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: fIm_err must be a number.")
            # if line_split[0].lower()=='d0re_goal':
            #     if len(line_split)!=2:
            #         sys.exit('Wrong set up of d0Re_goal')
            #     try:
            #         d0Re_goal = float(line_split[1])  # Try converting the string to a float
            #     except ValueError:
            #         sys.exit("Error: d0Re_goal must be a number.")
            # if line_split[0].lower()=='d0re_err':
            #     if len(line_split)!=2:
            #         sys.exit('Wrong set up of d0Re_err')
            #     try:
            #         d0Re_err = float(line_split[1])  # Try converting the string to a float
            #     except ValueError:
            #         sys.exit("Error: d0Re_err must be a number.")         
            if line_split[0].lower()=='par1':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par1')
                try:
                    par1_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par1 must be a number.")
                if len(line_split)==3:
                    try:
                        par1_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par1 must be a number.")
                else:
                    par1_range[1] = par1_range[0]
            if line_split[0].lower()=='par2':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par2')
                try:
                    par2_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par2 must be a number.")
                if len(line_split)==3:
                    try:
                        par2_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par2 must be a number.")
                else:
                    par2_range[1] = par2_range[0]
            if line_split[0].lower()=='par3':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par3')
                try:
                    par3_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par3 must be a number.")
                if len(line_split)==3:
                    try:
                        par3_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par3 must be a number.")
                else:
                    par3_range[1] = par3_range[0]
            if line_split[0].lower()=='par4':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par4')
                try:
                    par4_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par4 must be a number.")
                if len(line_split)==3:
                    try:
                        par4_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par4 must be a number.")
                else:
                    par4_range[1] = par4_range[0]
            if line_split[0].lower()=='par5':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par5')
                try:
                    par5_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par5 must be a number.")
                if len(line_split)==3:
                    try:
                        par5_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par5 must be a number.")
                else:
                    par5_range[1] = par5_range[0]
            if line_split[0].lower()=='par6':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par6')
                try:
                    par6_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par6 must be a number.")
                if len(line_split)==3:
                    try:
                        par6_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par6 must be a number.")
                else:
                    par6_range[1] = par6_range[0]
            if line_split[0].lower()=='m1':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of M1')
                try:
                    M1 = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: M1 must be a number.")
                if M1<=0:
                    sys.exit("Error: M1 must be non-zero and positive.")
            if line_split[0].lower()=='m2':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of M2')
                try:
                    M2 = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: M2 must be a number.")
                if M2<=0:
                    sys.exit("Error: M2 must be non-zero and positive.")
            if line_split[0].lower()=='kmin':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of kMin')
                try:
                    kMin = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: kMin must be a number.")
                if kMin<0:
                    sys.exit("Error: kMin must be non-negative.")
            if line_split[0].lower()=='kmax':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of kMax')
                try:
                    kMax = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: kMax must be a number.")
                if kMax<=0:
                    sys.exit("Error: kMax must be positive.")
            if line_split[0].lower()=='kbin':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of kBin')
                if not line_split[1].isdigit():
                    sys.exit('kBin is not a digit')
                kBin = int(line_split[1])
                if kBin<=0:
                    sys.exit("Error: kBin must be positive.")
            if line_split[0].lower()=='eps':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of eps')
                try:
                    eps = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: eps must be a number.")
                if eps<1e-10 or eps>1e-5:
                    sys.exit("Error: eps must be >1e-10 and <1e-5 for CATS to properly work.")
            if line_split[0].lower()=='pot':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of pot')
                pot = line_split[1]
                #!! DEFINE THE LIST OF POTENTIALS
                if pot not in ["ComplexGaussian"]:
                    sys.exit("Error: pot should be selected from 'ComplexGaussian'")
            if line_split[0].lower()=='pw':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of pw')
                pw_str = line_split[1]
                if pw_str.lower() not in ["s","p","d","f"]:
                    sys.exit("Error: pw should be selected from 's, p, d, f'")
                if pw_str.lower()=='s':
                    PW = 0
                if pw_str.lower()=='p':
                    PW = 1
                if pw_str.lower()=='d':
                    PW = 2
                if pw_str.lower()=='f':
                    PW = 3
            if line_split[0].lower()=='coulomb':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of Coulomb')
                if not line_split[1].isdigit():
                    sys.exit('Coulomb is not a digit')
                q1q2 = int(line_split[1])
            if line_split[0].lower()=='cpu':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of CPU')
                if not line_split[1].isdigit():
                    sys.exit('CPU is not a digit')
                CPU = int(line_split[1])
                if CPU<=0:
                    sys.exit("Error: CPU must be positive.")
                if CPU>32:
                    sys.exit("Error: CPU must be <= 32.")

            if line_split[0].lower()=='timeout':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of TIMEOUT')
                if not line_split[1].isdigit():
                    sys.exit('TIMEOUT is not a digit')
                TIMEOUT = int(line_split[1])
                if TIMEOUT<=0:
                    sys.exit("Error: TIMEOUT must be positive.")
            if line_split[0].lower()=='fit_timeout':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of FIT_TIMEOUT')
                if not line_split[1].isdigit():
                    sys.exit('FIT_TIMEOUT is not a digit')
                FIT_TIMEOUT = int(line_split[1])
                if FIT_TIMEOUT<=0:
                    sys.exit("Error: FIT_TIMEOUT must be positive.")
            if line_split[0].lower()=='wakeup':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of WakeUp')
                try:
                    WakeUp = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: WakeUp must be a number.")

        if par1_range[0]>par1_range[1]:
            sys.exit('Error: par1 minimum value is larger than the maximum.')
        if par2_range[0]>par2_range[1]:
            sys.exit('Error: par2 minimum value is larger than the maximum.')
        if par3_range[0]>par3_range[1]:
            sys.exit('Error: par3 minimum value is larger than the maximum.')
        if par4_range[0]>par4_range[1]:
            sys.exit('Error: par4 minimum value is larger than the maximum.')
        if par5_range[0]>par5_range[1]:
            sys.exit('Error: par5 minimum value is larger than the maximum.')
        if par6_range[0]>par6_range[1]:
            sys.exit('Error: par6 minimum value is larger than the maximum.')

        #!! YOU CAN MAKE SOME SAFETY FEATURES HERE, OR REMOVE THOSE IF THEY MAKE NO SENSE
        # if par2_range[0] <= 0:
        # sys.exit('par2 is the width of the potential and should be non-zero and positive.')

        if pot not in ["ComplexGaussian"]:
            sys.exit("Error: pot should be selected from 'ComplexGaussian'")
        if M1<=0 or M2<=0:
            sys.exit("Error: Both M1 and M2 should be defined and be strictly positive.")

        if kMin<0:
            sys.exit("Error: kMin cannot be negative!")
        if kMin>100:
            sys.exit("Error: kMin should be small, certainly not above 100 MeV!")
        if kMin>=kMax:
            sys.exit("Error: kMin should be smaller than kMax!")
        if kMaxFit>kMax:
            kMaxFit = kMax
        if kMax<0:
            sys.exit("Error: kMax cannot be negative!")
        if kBin<=0:
            sys.exit("Error: kBin should be positive!")

        #!! BASED ON PASSED EXPERINENCE, some minimum error is required, else you will stall
        if fRe_err==0:
            fRe_err = abs(0.01*fRe_goal)
            if fRe_err < 0.01:
                fRe_err = 0.01
        if fIm_err==0:
            fIm_err = abs(0.01*fIm_goal)
            if fIm_err < 0.01:
                fIm_err = 0.01
        # if d0Re_err==0:
        #     d0Re_err = abs(0.01*d0Re_goal)
        #     if d0Re_err < 0.01:
        #         d0Re_err = 0.01

        fRe_err_current = fRe_err*10
        fIm_err_current = fIm_err*10
        # d0Re_err_current = d0Re_err*10
        
        UniqueID = 0
        BestEstimator = 1e64
        Best_fRe = 1e64
        Best_fIm = 1e64
        # Best_d0Re = 1e64
        Best_par1 = 0
        Best_par2 = 0
        Best_par3 = 0
        Best_par4 = 0
        Best_par5 = 0
        Best_par6 = 0
        Progress_fRe = 0
        Progress_fIm = 0
        Progress_d0Re = 0        
        TotExeTime = 0
        StartTime = time.time()/60.
        # we re-iterate until we reach our convergence criteria
        while (fRe_err_current>fRe_err or fIm_err_current>fIm_err) and TotExeTime<TIMEOUT:
            TrialsCPU = []
            FilesIn = []
            FilesOut = []
            CurrentPar1 = []
            CurrentPar2 = []
            CurrentPar3 = []
            CurrentPar4 = []
            CurrentPar5 = []
            CurrentPar6 = []
            for iCPU in range(0,CPU):
                UniqueID += 1
                fileBaseJob = BaseFileName+"_"+str(UniqueID)
                #!! This is the input for the cpp script
                fileNameIn = fileBaseJob+".fit"
                #!! This is the output from the cpp script
                fileNameOut = fileBaseJob+".root"
                FilesIn.append(fileNameIn)
                FilesOut.append(fileNameOut)
                # EstimatorsCPU.append(0)

                trial = studyT.ask()
                TrialsCPU.append(trial)
                target_function(trial)

                try:
                    shutil.copy(filename_arg, fileNameIn)
                except FileNotFoundError:
                    print(f"Source file '{filename_arg}' not found.")
                except Exception as e:
                    sys.exit(f"An error occurred: {e}")

                # Check if the file exists
                if os.path.exists(fileNameOut):
                    try:
                        # Delete the file
                        os.remove(fileNameOut)
                    except OSError:
                        pass  # Ignore errors silently if they occur during file deletion

                with open(fileNameIn, 'a') as fjob:
                    fjob.write( "UniqueID   "+str(UniqueID)+"\n")
                    fjob.write( "v_par1     "+str(trial.params['v_par1'])+"\n")
                    fjob.write( "v_par2     "+str(trial.params['v_par2'])+"\n")
                    fjob.write( "v_par3     "+str(trial.params['v_par3'])+"\n")
                    fjob.write( "v_par4     "+str(trial.params['v_par4'])+"\n")
                    fjob.write( "v_par5     "+str(trial.params['v_par5'])+"\n")
                    fjob.write( "v_par6     "+str(trial.params['v_par6'])+"\n")

                    CurrentPar1.append(trial.params['v_par1'])
                    CurrentPar2.append(trial.params['v_par2'])
                    CurrentPar3.append(trial.params['v_par3'])
                    CurrentPar4.append(trial.params['v_par4'])
                    CurrentPar5.append(trial.params['v_par5'])
                    CurrentPar6.append(trial.params['v_par6'])

                run_script(fileBaseJob)
            # for iCPU in range(0,CPU)

            # wait for all jobs to finish
            StartWaitTime = time.time()/60.
            NumRunningJobs = CPU
            while NumRunningJobs>0 and (time.time()/60.-StartWaitTime)<FIT_TIMEOUT:
                NumRunningJobs = CPU
                for iCPU in range(0,CPU):
                    if os.path.exists(FilesOut[iCPU]):
                        try:
                            NumRunningJobs -= 1
                            time.sleep(WakeUp*0.1)
                        except OSError:
                            time.sleep(WakeUp)
            time.sleep(0.5)
            if NumRunningJobs>0:
                print('WARNING: Some jobs have failed')
            # after all jobs are done, we collect the output
            for iCPU in range(0,CPU):
                if os.path.exists(FilesOut[iCPU]):
                    try:
                        #!! the root file provided by the executable should have a histogram, where bin 1 is the fRe and bin 2 is the fIm
                        root_file = ROOT.TFile.Open(FilesOut[iCPU])
                        hScattPars = root_file.Get("hScattPars")
                        hPotentialReal = root_file.Get("hPotentialReal")
                        hPotentialImag = root_file.Get("hPotentialImag")

                        Estimator = 0.0
                        Current_fRe = hScattPars.GetBinContent(1)
                        Current_fIm = hScattPars.GetBinContent(2)
                        # Current_d0Re = hScattPars.GetBinContent(3)
                        Estimator += ( (Current_fRe-fRe_goal)/fRe_err )**2
                        Estimator += ( (Current_fIm-fIm_goal)/fIm_err )**2
                        # Estimator += ( (Current_d0Re-d0Re_goal)/d0Re_err )**2

                        if BestEstimator>Estimator:

                            BestEstimator = Estimator
                            Best_par1 = CurrentPar1[iCPU]
                            Best_par2 = CurrentPar2[iCPU]
                            Best_par3 = CurrentPar3[iCPU]
                            Best_par4 = CurrentPar4[iCPU]
                            Best_par5 = CurrentPar5[iCPU]
                            Best_par6 = CurrentPar6[iCPU]                            

                            Best_fRe = Current_fRe
                            Best_fIm = Current_fIm
                            # Best_d0Re = Current_d0Re
                            
                            fRe_err_current = abs(Best_fRe-fRe_goal)
                            fIm_err_current = abs(Best_fIm-fIm_goal)
                            # d0Re_err_current = abs(Best_d0Re-d0Re_goal)
                            # print("fRe_err_current = ",fRe_err_current,"\n")
                            # print("fIm_err_current = ",fIm_err_current,"\n")

                            if fRe_err_current==0:
                                Progress_fRe = 1
                            else:
                                Progress_fRe = fRe_err/fRe_err_current
                            if fIm_err_current==0: 
                                Progress_fIm = 1
                            else:
                                Progress_fIm = fIm_err/fIm_err_current
                                
                            # if d0Re_err_current==0:
                            #     Progress_d0Re = 1
                            # else:
                            #     Progress_d0Re = d0Re_err/d0Re_err_current

                            # The message you want to update
                            if Progress_fRe>1:
                                Progress_fRe = 1
                            if Progress_fIm>1:
                                Progress_fIm = 1
                            # if Progress_d0Re>1:
                            #     Progress_d0Re = 1
                            print(f"Progress: fRe({Progress_fRe*100:.1f} %) = {Best_fRe:.3f} fm     fIm({Progress_fIm*100:.1f} %) = {Best_fIm:.3f} fm", end="\r")

                            best_file = ROOT.TFile.Open(BaseFileName+"_BestSolution.root","recreate")
                            hScattPars.Write()
                            hPotentialReal.Write()
                            hPotentialImag.Write()
                            root_file.cd()
                            shutil.copy(FilesIn[iCPU], BaseFileName+"_BestSolution.txt")
                            #!! SAVE FULL INFO, INCLUDING BEST PARAMETERS (needs to be updated)
                            with open(BaseFileName+"_BestSolution.txt", 'a') as fupdate:
                                fupdate.write( "fRe_best     "+str(Best_fRe)+"\n")
                                fupdate.write( "fIm_best     "+str(Best_fIm)+"\n")
                                # fupdate.write( "d0Re_best     "+str(Best_d0Re)+"\n")

                        studyT.tell(TrialsCPU[iCPU],Estimator)
                        root_file.Close()
                        # clean up
                        os.remove(FilesOut[iCPU])
                        os.remove(FilesIn[iCPU])
                    except OSError:
                        pass
            # for iCPU in range(0,CPU)
            TotExeTime = time.time()/60.-StartTime
        # while f_err_current>f_err or d_err_current>d_err

        print("\nSUMMARY:")
        if fRe_err_current<=fRe_err and fIm_err_current<=fIm_err:
            print(" - The desired convergence was reached!")
        else:
            print(" - The script terminated due to timeout!")
        print(f" - Obtained fRe = {Best_fRe:.3f} fm")
        print(f" - Obtained fIm = {Best_fIm:.3f} fm")
        # print(f" - Obtained d0Re = {Best_d0Re:.3f} fm")
        print(" - The potential was "+pot)
        print(f"   par1 = {Best_par1:.2f}")
        print(f"   par2 = {Best_par2:.5f}")
        print(f"   par3 = {Best_par3:.2f}")
        print(f"   par4 = {Best_par4:.5f}")
        print(f"   par5 = {Best_par5:.2f}")
        print(f"   par6 = {Best_par6:.5f}")        
        print(" - Detailed output in the "+BaseFileName+"_BestSolution"+".txt/root files.")


if __name__ == "__main__":
    main()
