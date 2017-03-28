# Chemical Equilibrium with Applications for MATLAB

This project aims to modify NASA's [Chemical Equilibrium with Applications](https://www.grc.nasa.gov/WWW/CEAWeb/) to be able to use it's functionalities with MATLAB in a much easier and powerful way.

## MEX
The **_MEX_** folder holds the modified cea2.f file with it's compiled mex64 file and the CEA class to set the input parameters and run the program properly. This folder also holds an example file called *initialize.m*.

### <span style="color:red;"> IMPORTANT </span>
*The mex version on CEA only supports the **rocket** application at this time.*

## DEV
The **_DEV_** folder holds doccumentation for the developement of the API and the tracking of changes of the MEX files.

## OLD
The **_OLD_** folder holds the following files which will be deleted in the comming weeks

>In the **_OLD/NasaSource_** folder there is all the original files downloaded straight from NASA's CEA website. It was supposed to be better commented, however I never got around to doing it.

>The **_OLD/CEA_Wrapper_** folder holds a quickly made matlab wrapper of CEA which consists of a function that writes the input file, calls the compiled CEA program, and translates the output file into a struct matrix which is returned to the user. Will be deleted since the MEX method is faster and much more accurate for the output data.

Created as a project in collaboration with École Polytechnique de Montréal