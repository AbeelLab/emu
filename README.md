# emu
Emu is a tool that canonicalizes non-identical representations of large insertions and deletions (indels) and substitutions across multiple genome samples. All that is needed to get started is a file containing a list of VCF files and the appropriate reference genome used to generate the VCF files.

Quick-run recipe:

  `java -jar emu.jar extract -o <output directory> -i <file containing list of VCF files> -r <FASTA-formatted reference genome>`

  `java -jar emu.jar canonicalize -o <output directory> -i ../AllExtractedLSVs.reference_extended.vcf`

  `java -jar emu.jar vcfify -o <output directory> -c ../CanonicalizedLSVs.vcf -s ../AllExtractedSNVs.vcf -v <file containing list of VCF files>`
