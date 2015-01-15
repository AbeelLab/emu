package abeellab.emu

import atk.util.Tool
import java.io.File


object Emu {

  /**
   * input: string of VCF files directory, directory where to write output
   * passes inputs to pipeline
   *
   */

    def main(args: Array[String]): Unit = {

      if (args.length == 0) {

        println(listInstructions);println(help)
      } else {
        args(0) match {
          case "list" => println(listInstructions)
          case "help" => println(listInstructions);println(help)          
          case "extract" => println("Extracting large sequence variants from VCF files."); ExtractLSVs.main(args.drop(1))
          case "merge-identicals" => println("Merging identical LSVs."); MergeIdenticalLSVs.main(args.drop(1))                               
          case "canonicalize" => println("Canonicalizing large sequence variants.");Canonicalize.main(args.drop(1))                   
          case "lsv-matrix" => println("Creating LSV-to-sample matrix file."); LSVmatrix.main(args.drop(1))
          case "metrics" => println("Computing Emu metrics."); EmuMetrics.main(args.drop(1))
          case "vcfify" => println("VCFifying Emu output files."); VCFify.main(args.drop(1))
          case "fasta2vcf" => println("Converting FASTA files to VCF."); FASTA2VCF.main(args.drop(1))
          case _ => println(listInstructions); println(help)
        }
      }

    }
    
     val listInstructions = ("Usage: java -jar emu.jar [command] [command parameters...]\n")
     
     val help = (
         
        //Thanks to Text ASCII Art Generator (http://patorjk.com/blog/software/)
"""
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
         
"""
         +
"Emu is an algorithm that canonicalizes alternate representations of large sequence variants across microbial genomes.\n"
         +

"Please contact Alex Salazar at salazar@broadinstitute.org for any questions.\n\n"         
         +
"    COMMANDS             DESCRIPTION\n"
         +
"    extract              Extracts all PASS-only LSVs and SNPs.\n"+
"    canonicalize         Identifies and canonicalizes alternate representations of LSVs across a set of microbial genomes.\n"+
"    lsv-matrix           Creates a LSV-to-sample matrix file\n"+
"    vcfify               Creates VCF files containing normalized LSVs and extracted SNPs.\n"+ //+
"    merge-identicals     Merges identical uncanonicalized LSVs across a set of microbial genomes.\n"+
"    metrics              Creates data file to summarize emu's canonicalization performance.\n"+
"    lsv-matrix           Creates an LSV-to-sample matrix."

)
println()

}