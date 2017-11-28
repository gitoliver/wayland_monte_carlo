# wayland_monte_carlo
Instructions:
Install gems and gmml. As of Nov 16, 2017 you need the dev version of both i.e. clone -b gems-dev and clone-b gmml-dev

As part of installing gems you will have set GEMSHOME

Also do
export LD_LIBRARY_PATH=$GEMSHOME/gmml/bin/

open main.cpp and edit this line:
#include "/home/asdf/gems/gmml/includes/MolecularModeling/assembly.hpp"
To reflect your system path. I hope to figure out how to not have this as a requirement.

Now you can compile:
g++ -std=c++0x -I$GEMSHOME/gmml/includes/* -L$GEMSHOME/gmml//bin/ main.cpp -lgmml -o main.exe
