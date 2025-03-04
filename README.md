iCATS is an extended functionality based on the CATS framework aiming at delivering the two-particle correlation function for any local complex potential and for any emitting source profile.
iCATS exploit all the functionalities of CATS, hence please refer to the main code for any additional details.
In order to use this novel functionality, please use the following setter function on your CATS object (Kitty):

Kitty.SetUsingiCATS(true)

By default this function is set to false.
Examples of how to implement a complex potential can be found in the CATS_Extensions/DLM_Potentials.cpp(h) for a Gaussian Complex Yukawa-type potential (https://github.com/dimihayl/DLM/blob/devImCats/CATS_Extentions/DLM_Potentials.cpp#L1900).
