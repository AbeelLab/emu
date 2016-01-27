package abeellab.emu

import atk.util.Tool
import java.io.File


object CLI {

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
          case "help" => println(listInstructions);println(help)          
          case "extract" => println("Extracting large sequence variants from VCF files."); ExtractLSVs.main(args.drop(1))
          case "normalize" => println("Normalizing large variants."); Normalize.main(args.drop(1))
          case "vcfify" => println ("Creating VCF files for each sample"); VCFify.main(args.drop(1))
        
          
          case _ => println(help)
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
"    normalize         Identifies and canonicalizes alternate representations of LSVs across a set of microbial genomes.\n"+
"    vcfify               Creates VCF files containing normalized LSVs and extracted SNPs.\n"


)
println()

}