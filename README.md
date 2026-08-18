# Sequence Dotplot

Python script for plotting DNA sequence dotplots, similar to https://en.vectorbuilder.com/tool/sequence-dot-plot.html. 

## Installation 

Clone the repository locally and run from your computer

''' bash 
git clone 
'''

The script requires numpy and matplotlib. 

'''bash 
pip install numpy
pip install matpotlib 
'''

## Usage

* **dotplot.py** is the main python script - run this to produce your plot
* **sequence1.txt** and **sequence2.txt** are text files where you should input the sequence you would like to compare
* **options.txt** are user configurable options
  * Window size controls how many nucleotides need to match for a dot to be plotted
  * You can optionally add you own gene annotations in the format label:start,end
 
The files currently contain an example input for the Tasmanian Devil MHC locus which produces the following plot:
<img width="1689" height="918" alt="image" src="https://github.com/user-attachments/assets/28ccb1bc-4e98-46e5-a014-027977b46f08" />
