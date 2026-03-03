# find_and_extract_location_in_assembly
Python script which extracts fasta files from genome assemblies using minimap2 alignment results

This script takes two main parameters corresponding to the input assembly folder in which assemblies are in fasta format and can have the .fa, .fasta or .fna extensions and can be gzip compressed or not, and a reference folder in which are the genome portions (called reference in the script help) in fasta format which the user wants to locate in the assemblies. The naming convention for assemblies is Assembly_name.hapX.extension and for the reference Reference_name.extension.

The script will align each reference on each assembly and merge alignments using a distance threshold (default 10kb). It will find the longuest alignment location and extract it in fasta format. This fasta files will be names Assemblyname.hapX_Referencename.fa. 

The typical use case for the script is to extract the immunoglobuline loci from a set of assemblies. 

The script has several other parameters preseented hereunder. 

<pre>
  usage: process_references_assemblies.py [-h] -a ASSEMBLIES [-t THREADS] -r REFERENCES [-p PAF] [-d DIST] [-o OUTPUT] [-x REMOVEPAF] [-f FLANQUING]

Align references to assemblies and extract best hits.

options:
  -h, --help            show this help message and exit
  -a ASSEMBLIES, --assemblies ASSEMBLIES
                        Directory containing genome assemblies
  -t THREADS, --threads THREADS
                        Number of threads for minimap2
  -r REFERENCES, --references REFERENCES
                        Directory containing reference sequences
  -p PAF, --paf PAF     Directory in which paf files will be written
  -d DIST, --dist DIST  Max distance to merge alignments, default 10Kb
  -o OUTPUT, --output OUTPUT
                        Output summary file
  -x REMOVEPAF, --removepaf REMOVEPAF
                        Removes paf intermediate paf files
  -f FLANQUING, --flanquing FLANQUING
                        Flanking region size to add

</pre>
