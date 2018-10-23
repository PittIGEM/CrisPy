# CrisPy
Software tool to analyze Sanger sequencing readouts of CRISPR experiments

The aim of our software tool (CrisPy) was to provide a fast and accurate way for characterizing CRISPR on-target and off-target edits through Sanger sequencing. Software tools have attempted to numerically characterize DNA traces in the past. Mark Crowe’s SeqDoc (Perl, 2005) provided functionality for normalizing and comparing two traces, while Timothy Lu’s Sequalizer (Matlab, 2017) used these comparisons to estimate nucleotide mutation frequencies. Our system builds upon this previous work by utilizing modified versions of these algorithms in a free and concise python module. Moreover, this module creates new functionality by predicting and ranking likely off-target sequences (using Gibbs free energy) and finding their mutation frequencies. We tested this tool using a plasmid containing purposeful off-target sequences and consequently gained insight into our recording circuit’s behavior. It is our hope that synthetic biologists can use this tool to quickly test for off-target mutation frequencies that exist in their system.



To set up:
1. Download the biopython library https://biopython.org/wiki/Download
2. Download this repo and run CrisPyApp.py to use the GUI, or directly use classes defined in CrisPy.py
3. Upload your reference and test trace (.ab1 files) to compare the two
4. Specify your target sequence and base pairs that will be mutated (if nothing entered, will default to whatever is predefined in CrisPyApp) 
5. Analyze with results printed to screen. The target sequence will be displayed first, along with its mutation frequencies; next off-targets will be displayed in order of the normalized match score (higher is a closer match).




Email Evan Becker (ewb12@pitt.edu) for questions



Testing Results:

![alt text](https://raw.githubusercontent.com/MiscM/CrisPy/master/CrisPy.png)

