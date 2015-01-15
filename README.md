# emu
Emu is a tool that canonicalizes non-identical representations of large insertions and deletions (indels) and substitutions across multiple genome samples. All that is needed to get started is a file containing a list of VCF files and the appropriate reference genome used to generate the VCF files.

Quick-run recipe:

  `java -jar emu.jar extract -o <output directory> -i <file containing list of VCF files> -r <FASTA-formatted reference genome>`

  `java -jar emu.jar canonicalize -o <output directory> -i ../AllExtractedLSVs.reference_extended.vcf`

  `java -jar emu.jar vcfify -o <output directory> -c ../CanonicalizedLSVs.vcf -s ../AllExtractedSNVs.vcf -v <file containing list of VCF files>`

# Commands:

1). extract

java -jar emu.jar extract -o <output directory> -i <file containing list of VCF files>

Emu will extract all PASS-only large sequence variants (LSVs) and small variants such as SNPs and indels.
To run this command, Emu requires an ouput directory and a file containing a list of VCF files. There are two outputfiles: AllExtractedSmall.emu and AllExtractedLSPs.emu

2). normalize

java -jar emu.jar normalize -o <output directory> -i ../AllExtractedLSPs.emu -t <boolean value>

Emu will perform cross-genome comparisons of the extracted LSVs and differentiate between alternate representations and unique LSVs.
To run this command, Emu requires an output directory and the "AllExtractedLSPs.emu" file which is outputtted by the "extract" command.
Optionally, a boolean value in the "-t" field will make Emu create a subdirectory inside the output directory called "normalization_log" which will contain a log file for each LSV processed through Emu. The log file will contain a summary of which LSVs (if any) were canonicalized and which weren't. 
The output file is "Normalized.emu" and contains every unique LSV from the samples inputted. 

3). vcfify

java -jar emu.jar vcfify -o <output directory> -n ../Normalized.emu -s ../AllExtractedSmall.emu -v <file containing list of VCF files>

Emu will combine the "Normalize.emu" and "AllExtractedSmall.emu" files and reformat the data to VCF format.
To run this command, Emu requires an output directory, the "Normalized.emu" file, the "AllExtractedSmall.emu" file, and the file containing the list of VCF file used in "extract" command.
The output file is "<sample name>_Normalized.vcf" which is the new VCF file after resolving for alternate representations. 

Data preparation for excel visualization:

1). statistics

java -jar emu.jar statistics -o <output directory> -e ../AllExtractedLSPs.emu -n ../Normalized.emu

Emu will perform statistics for a side-by-side comparisons of the number of unique LSVs before and after resolving for alternate representations.
To run this command, Emu requires an ouput directory, the "AllExtractedLSPs.emu" file (output result from the "extract" command), and the "Normalized.emu" file (output result from "normalize" command). The output file is "Emu_Statistics.txt"; results can be visualized in Excel.

2). gwp

java -jar emu.jar gwp -o <output directory> -i ../Normalized.emu

Emu will perform data processing to visualize LSV type, size, and their genome-wide distribution.
To run this command, Emu requires an ouput directory and the "Normalized.emu" file (output result from "normalize" command).
The output files is "Gwp_data.txt" and can be visualized in Excel. 