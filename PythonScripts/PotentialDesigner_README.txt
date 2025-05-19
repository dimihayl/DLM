Running this script will allow you to create a potential, that has some specific phase shifts or scattering parameters. Below two examples, how to run the script for A) known phase shifts or B) known scattering parameters. This python script will execute multiple C++ CATS scripts in parallel, to do the fitting. The whole process is controlled by a settings file, described below. Important: there a limited potential types that are currently possible to use, if you need more, either implement within DLM_FemtoTools/CommonAnaFunctions.cpp as well as the python script in this folder, or contact the developers. The currently allowed potential types are:
------------
Gauss:
par1, par2 are the amplitude (MeV) and width (fm)
------------
DoubleGauss:
par1, par2 are the amplitude (MeV) and width (fm) of the first Gauss
par3, par4 are the amplitude (MeV) and width (fm) of the second Gauss
------------
Yukawa:
A/r exp(-alpha r)
par1 = A and is dimensionless
par2 = alpha in MeV
------------
YukawaDLM:
apporximately the same as Yukawa,
but modifies at very low distances to avoid a singularity for numerical stability
------------
UsmaniCore:
an Usmani potential (tuned to pLambda), where we can change the repulsive core:
par1 is the amplitude (MeV)
par2 is the range of the core (in fm)
par3 is the slope of the core (in fm)



A) --------------------------
To design a new potential, go to your work folder and run:
python3 <PATH>/PotentialDesigner.py <filename>

The <filename> is a settings file, an example of a such is provided below, where
we read the phase shifts from a TH1F saved within a ROOT file.
The number on the first line is the desired chi2 to reach
--------------------------
PhaseShifts	/home/input_path/random_root_file.root	h_TH1F_object_name	1
M1	938.272
M2	1115.683
kMin	0
kMax	90
kBin	45
par1	2134	2140
par2	0.1	0.6
par3	0.05	0.35
par4	0
pw	s
pot 	UsmaniCore
WakeUp  1
TIMEOUT 15
--------------------------


B) --------------------------
To design a new potential, go to your work folder and run:
python3 <PATH>/PotentialDesigner.py <filename>

The <filename> is a settings file, an example of a such is provided below:
--------------------------
f_goal	2.91
f_err	0.005
d_goal  2.71
d_err	0.2
M1	938.272
M2	1115.683
kMin	0
kMax	90
kBin	45
par1	2134	2140
par2	0.1	0.6
par3	0.05	0.35
par4	0
pw	s
pot 	UsmaniCore
WakeUp  1
TIMEOUT 15
--------------------------
Detailed description of these parameters:
f_goal +/- f_err is the scattering length you want to obtain
d_goal +/- d_err is the effective range you want to obtain
M1/M2 mass of the first/second particle
kMin/kMax is the fit range within the effective range expansion. Typically 0-90 is fine
kBin how many bins we have within the effective range expansion. Typically works well with a bin width of 2 MeV
par1/2/3/4 are the 4 parameters of the potential that this scripts allows you to fit. What they really are, depends on the type of potential
pw which pw is fitted (s,l,d,f)
pot is the potential type
WakeUp is how often (in sec) the script checks if the C++ code provided output. Typically 1 second is good.
TIMEOUT in minutes, is the max time to let the script run.
nPars optional parameter, the number of expansion parameters in the effective range expansion, relevant to B). Default 3.
CPU optional parameter, the number of threads to use by the system. By default it uses all available threads - 1.
