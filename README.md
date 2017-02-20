<h1>Chemical Equilibrium with Applications for MATLAB</h1>

<p>This project aims to modify NASA's <a href="https://www.grc.nasa.gov/WWW/CEAWeb/">Chemical Equilibrium with Applications</a> to be able to use it's functionalities with MATLAB in a much easier and powerful way.</p>

<p>In the <emph>NasaSource</emph> folder there is all the original files downloaded straight from NASA's CEA website. However the FORTRAN source codes are commented for better comprehension of the way CEA works.</p>

<p>The <emph>CEA_Wrapper</emph> folder holds a quickly made matlab wrapper of CEA which consists of a function that writes the input file, calls the compiled CEA program, and translates the output file into a struct matrix which is returned to the user.</p>

<p>The <emph>MEX</emph> folder holds the modified cea2.f file and a subfolder which will hold the mex64 files and their dependencies once done.</p>

<p>Created as a project in collaboration with École Polytechnique de Montréal</p>