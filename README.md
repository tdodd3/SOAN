# SOAN
Subsets Of Adjacent Nodes - A program for rapidly computing suboptimal paths in a protein dynamic network.

# How to use SOAN from the command line
SOAN is written in Python and is designed to be called from the command line. Here is a example 

python soan.py -i adjM.dat -s 19 -t 235 -l 1 -k 1000 -o paths.txt 

At any point, the list of command line parameters and default values can be viewed by,

python soan.py -h

usage: soan.py [-h] [-i INPUT] [-s SOURCE] [-t TARGET] [-l LEVEL] [-k NPATHS]
               [-o OUTPUT] [-m MAXLEVEL]

optional arguments:

-h, --help   show this help message and exit

-i INPUT     path to adjacency matrix input file, default='network.dat'

-s SOURCE    source node (can be multiple sources), default=1

-t TARGET    target node, default=1

-l LEVEL     level of neighbors to expand to, default=1

-k NPATHS    number of paths to find, default=1000

-o OUTPUT    output file name, default='paths.txt'

-m MAXLEVEL  Max level to expand to (default=5), this should only be changed
               if the calculation fails
