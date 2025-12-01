# Mass Spring System solver

## Usage
After cloning first run, run in the command line:

git submodule update --init

Then build the project with cmake:

mkdir build
cd build
cmake ..
make
cd ..


You can then run and plot the mass-spring system for a given time and nr of steps automatically with the shell script "runmassspring.sh":

~/team01$ ./runmassspring.sh tend_relativetopi steps method

tend_relativetopi and steps are doubles, method is either "explicit" or "improved", for example:


~/team01$ ./runmassspring.sh 4 100 explicit
