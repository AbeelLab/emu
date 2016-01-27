package abeellab.emu

import atk.util.Tool
import java.io.File


object ManuscriptCLI {

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
         
          case "identify-simulated-indels" => println("Identifying simulated indels."); IdentifyLSVs.main((args.drop(1)))
          case "validate-kmers" => println("Extracting and validating kmer uniquness from sequence files."); 
          							ValidateKmers.main((args.drop(1))) 
          case "modify-reference" => println("Modifying reference genome with simulated indels."); 
          								ModifyReferenceGenome.main((args.drop(1)))
          								
          case "print-distance-metrics" => println("Printing distance metrics for simulated indels.");
          									PrintBreakpointDistanceMetrics.main((args.drop(1)))
          									
          case "filter-subs" => println("Println initiating substitution filtering.");
          							FilterSubstitutions.main((args.drop(1)))
          							
          case "vcf2conserved" => ConservedRegions.main((args.drop(1)))
          case "nearest-variant-distance" => FindNearestDistanceRange.main((args.drop(1)))
          case "fasta2vcf" => fasta2vcf.main((args.drop(1)))
          case "infer-genotype" => InferGenotype.main((args.drop(1)))
          case "genotype-robustness" => GenotypeRobustness.main((args.drop(1)))
          case "check-batch-jobs" => CheckBatchJobs.main((args.drop(1)))
          case "flanking-densities" => FlankingDensities.main((args.drop(1)))
          case "link-densities" => LinkMissedAndDensities.main((args.drop(1)))
          case "check-consistency" => CheckConsistency.main((args.drop(1)))
          case "format-list" => FormatTarget.main((args.drop(1)))
          case _ => println(help)
        }
      }

    }
    
     val listInstructions = ("Usage: java -jar emu.jar [command] [command parameters...]\n")
     
     val help ="help"

}