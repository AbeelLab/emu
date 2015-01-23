# Emu
Emu is a tool that normalizes alternate representations of large insertions and deletions (indels) and substitutions across multiple genome samples. This is done through a process called canonicalization. To get running, Emu requires a file containing a list of VCF files and the appropriate reference genome used to generate those VCF files.

Quick-run recipe:

  `java -jar emu.jar extract -o <output directory> -i <file containing list of VCF files> -r <FASTA-formatted reference genome>`

  `java -jar emu.jar canonicalize -o <output directory> -i ../AllExtractedLSVs.reference_extended.vcf`

  `java -jar emu.jar vcfify -o <output directory> -c ../CanonicalizedLSVs.vcf -s ../AllExtractedSNVs.vcf -v <file containing list of VCF files>`

# Downloads and documentation

Versioned production releases: https://github.com/AbeelLab/emu/releases

Nightly development builds: http://abeellab.org/jenkins/emu/

Source-code: https://github.com/AbeelLab/emu.git

Submit bugs and features requests: https://github.com/AbeelLab/emu/issues

Documentation and test run: see sections below:

Citation: 
we have recently submitted a manuscript of Emu to the ISMB 2015 Conference Proceedings. For temporary citation, cite our published extended abstract (www.biomedcentral.com/1471-2105/16/S2/A8/abstract).

# Getting Started
## Installing Emu
Download the (latest) binary release (https://github.com/AbeelLab/emu/releases) and extract it in your directory of choice. To verify proper installation, execute 'javar -jar emu.jar'. You should see the following CLI:

```
EEEEEEEEEEEEEEEEEEEEEEMMMMMMMM               MMMMMMMMUUUUUUUU     UUUUUUUU
E::::::::::::::::::::EM:::::::M             M:::::::MU::::::U     U::::::U
E::::::::::::::::::::EM::::::::M           M::::::::MU::::::U     U::::::U
EE::::::EEEEEEEEE::::EM:::::::::M         M:::::::::MUU:::::U     U:::::UU
  E:::::E       EEEEEEM::::::::::M       M::::::::::M U:::::U     U:::::U
  E:::::E             M:::::::::::M     M:::::::::::M U:::::D     D:::::U
  E::::::EEEEEEEEEE   M:::::::M::::M   M::::M:::::::M U:::::D     D:::::U
  E:::::::::::::::E   M::::::M M::::M M::::M M::::::M U:::::D     D:::::U
  E:::::::::::::::E   M::::::M  M::::M::::M  M::::::M U:::::D     D:::::U
  E::::::EEEEEEEEEE   M::::::M   M:::::::M   M::::::M U:::::D     D:::::U
  E:::::E             M::::::M    M:::::M    M::::::M U:::::D     D:::::U
  E:::::E       EEEEEEM::::::M     MMMMM     M::::::M U::::::U   U::::::U
EE::::::EEEEEEEE:::::EM::::::M               M::::::M U:::::::UUU:::::::U
E::::::::::::::::::::EM::::::M               M::::::M  UU:::::::::::::UU
E::::::::::::::::::::EM::::::M               M::::::M    UU:::::::::UU
EEEEEEEEEEEEEEEEEEEEEEMMMMMMMM               MMMMMMMM      UUUUUUUUU

Emu is an algorithm that canonicalizes alternate representations of large sequence variants across microbial genomes.

    COMMANDS             DESCRIPTION
    extract              Extracts all PASS-only LSVs and SNPs.
    canonicalize         Identifies and canonicalizes alternate representations of LSVs across a set of microbial genomes.
    lsv-matrix           Creates a LSV-to-sample matrix file
    vcfify               Creates VCF files containing normalized LSVs and extracted SNPs.
    merge-identicals     Merges identical uncanonicalized LSVs across a set of microbial genomes.
    metrics              Creates data file to summarize emu's canonicalization performance.
    lsv-matrix           Creates an LSV-to-sample matrix.
```

##Running Emu
###Extracting large variants
To extract large variants use the `extract` command. 

During extraction, flanking sequences are respectively added to the reference and alternative sequences of identified large variants. SNPs are redirected a separate output file called, `AllExtractedSNPs.vcf`.

Required arguments:
```
-o or --output
  Output directory
  
-i or --input
  File containing list of VCF file paths, one per line. NOTE: each pathname must end in the string '.vcf'.

-r or --reference-genome
  Input file that contains a list of FASTA-formatted reference genomes used in the VCF files. 
  Emu can parse multi-contig reference genomes.
```	        
Optional arguments:
```
-s or --LSV-size
  Minimum size for a large variant to be considered 'large' (default value is 2nt).
  
-F or --flanking-size
  Input flanking sequence size (nt) that will be attached to LSVs (default is 100). 
  This is generally used to analyze insertions in context of the reference genome.
  
-L or --LSV-output-prefix
  Prefix name for extacted LSV output file (default is 'AllExtractedLSVs.reference_extended.vcf').
```
###Normalizing large variants
To normalize large variants, use the `canonicalize` argument. 

Canonicalization is performed in two sets: variants that are complete or have at least X nt (see `--min-comparison-length argument`) and variants that are incomplete and have less than X nt. The X value is 20 nt by default. The latter set undergoes a much stricter criterion during canonicalization (e.g. large variants must have the same coordinates, type, reference, and alternative sequences to be canonicalized). The output of both sets are merged into one VCF file called 'CanonicalizedLSVs.vcf'

For information regarding canonicalization of large variants with nested variation, see the optional argument, `--max_leniency`.

For information regarding log files, see the optional argument, `--document-logs`.

Required arguments:
```
-o or --output
  Output directory.
  
-i or --input
  Input file extracted LSVs only. By default, this file is called 'AllExtractedLSVs.reference_extended.vcf'.
  
-l or --max_leniency
  Input maximum percent difference (default is 0.01). This can be tuned down to avoid overzealous 
  canonicalization of variants with nested variation (e.g. insertion and substitutions with SNPs).
```

Optional arguments:
```
-c or --min-comparison-length <value>
  Canonicalization is performed in two sets, variants that are complete (default is 20nt).
  
-p or --output-prefix
  Prefix name for canonicalized LSV output file (default is 'CanonicalizedLSVs').
  
-d or --document-logs
  Turn on/off documentation of every canonicalized event (this parameter is turned on by default).
  Each canonicalized event will contain a log file with information of raw variants that were merged,
  if any. Every canonicalized even is given an ID which is the name of the log files. They can be 
  found in the 'normalization_log' directory created inside the given output directory.

-s or --max-sub-size
  Largest substitution size allowed for normalization (default is 14999nt).
  
-x or --cross-contig
  Argument that allows  canonicalization of large variants across contigs. By default, this
  parameter is turned off.
```
###VCFify
Finally, the command `vcfify` will create individual standard VCF files that contain a strain's normalized large variants and initial SNPs. This command requires the following arguments:
```
-o or --output
  Output directory
  
-c or --canonicalized
  Input file that contains normalized variants (generally CanonicalizedLSVs.vcf).
  
-s or --small
  Input file that contains SNPs and indels. This file is an output of (AllExtractedSmall.emu)
  
-v or --vcf-list
  Input file that contains list of initial VCF files (same file given in the extract command).
  
-p or --output-prefix
  Ouput prefix for each VCF file (Default is 'canonicalized').
```
###Emu metrics
We provide pair of commands to summarize normalization in a given data set. The steps are outlined below:

First, use the `merge-identicals` command which merges raw large variants that were consistently predicted across a given data set and redirects the output to a separate VCF file. This command needs only two arguments:
```
-o or --output
  Output directory
  
-u or --uncanonicalized
  Input file containing raw large variant predictions (generally AllExtractedLSVS.reference_extended.vcf).
```

Finally, use the `metrics` command which compares the merged large variant (raw) predictions file with the normalized VCF file and outputs a text file containing general metrics summarizing canonicalization in a given data set. In other words, summarizes the number of unique calls before and after the application of Emu on a given data set. This commands requires three arguments:
```
-o or --output
  Input output directory.
  
-u or --uncanonicalized
  Input VCF file that contains uncanonicalized variants.
  
-c or --canonicalized
  Input VCF file that contains canonicalized variants
```
##Test run
To test run Emu, we provide two separate VCF files containing extracted large variants and SNPs (the output files of Emu's `extract` command) from three previously published clinical strains of Mycobacterium tuberculosis (http://www.nature.com/ng/journal/v45/n10/full/ng.2735.html). 

###Instructions:
Download the (latest) Emu binary release and extract it in a directory of your choice https://github.com/AbeelLab/emu/releases

Extract the two provided files to a directory of your choice (download 'Test-Files.zip' at https://github.com/AbeelLab/emu/releases). 

Create a txt file containing three arbitrary sample names, one per line. In real usage, this file should instead contain the full pathname that points to the VCF files from a desired data set. But for now, use arbitrary names like `Sample1.vcf, Sample2.vcf, and Sample3.vcf`. NOTE: each name must end in the string '.vcf'.

Execute the following command (note: -Xmx value was arbitrarily chosen):
`java -Xmx2g -jar <fill emu directory>.jar <fill output directory> canonicalize -o <fill output directory> -i AllExtractedLSVs.reference_extended.vcf`

Create separate VCF files for the arbitrary sample names:
`java -Xmx2g -jar <fill emu directory>.jar <fill output directory> vcfify -o <fill appropriate directory> -c <fill appropriate directory>/CanonicalizedLSVs.vcf -s <fill appropriate directory>/AllExtractedSNVs.vcf -v <fill appropriate directory>/vcfList.testfiles.txt`

