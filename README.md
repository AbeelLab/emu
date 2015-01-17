# Emu
Emu is a tool that canonicalizes non-identical representations of large insertions and deletions (indels) and substitutions across multiple genome samples. It requires file containing a list of VCF files and the appropriate reference genome used to generate the VCF files.

Quick-run recipe:

  `java -jar emu.jar extract -o <output directory> -i <file containing list of VCF files> -r <FASTA-formatted reference genome>`

  `java -jar emu.jar canonicalize -o <output directory> -i ../AllExtractedLSVs.reference_extended.vcf`

  `java -jar emu.jar vcfify -o <output directory> -c ../CanonicalizedLSVs.vcf -s ../AllExtractedSNVs.vcf -v <file containing list of VCF files>`

# Downloads and documentation

Versioned production releases: https://github.com/AbeelLab/emu/releases

Nightly development builds: http://abeellab.org/jenkins/emu/

Source-code: https://github.com/AbeelLab/emu.git

Submit bugs and features requests: https://github.com/AbeelLab/emu/issues

Documentation: below

# Getting Started
## Installing Emu
Download the (latest) binary release (https://github.com/AbeelLab/emu/blob/master/README.md) and extract it in your directory of choice. To verify proper installation, execute 'javar -jar emu.jar'. You should see the following CLI:

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
Please contact Alex Salazar at salazar@broadinstitute.org for any questions.

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
  File containing list of VCF file paths, one per line.

-r <value> | --reference-genome <value>
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
###Canonicalizing large variants
To canonicalize large variants, use the `canonicalize` argument. 

Canonicalization is performed in two sets: variants that are complete or have at least X nt (see `--min-comparison-length argument`) and variants that are incomplete and have less than X nt. The X value is 20 nt by default. The latter set undergoes a much stricter criterion during canonicalization (e.g. large variants must have the same coordinates, type, reference, and alternative sequences to be canonicalized).

For information regarding canonicalization of large variants with nested variation, see the `--max_leniency` argument. 
For information regarding log files, see the `--document-logs` argument.

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
