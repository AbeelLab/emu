package abeellab.emu

import atk.util.Tool
import java.io.File
import atk.compbio.vcf.VCFLine
import java.io.PrintWriter
import be.abeel.bioinformatics.FastaIterator

object VariantOverlapPartioner extends Tool with EmuUtils {

  case class Config(
    val largeVariants: File = null,
    val windowSize: Int = 10,
    val reference: File = null,
    val intDir: File = null,
    val onlyNormalize: Boolean = false,
    val leniency: Double = 0.1,
    val outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar emu.jar partition") {
      //command line arguments (see 'text' for description)
      opt[File]('l', "large-variants") required () action { (x, c) =>
        c.copy(largeVariants = x)
      } text ("Single vcf file containing all large variants (output of 'ExtractLSVs')")
      opt[File]('r', "reference-genome") required () action { (x, c) =>
        c.copy(reference = x)
      } text ("Reference genome that variants were called against.")
      opt[File]('i', "intermediate-directory") required () action { (x, c) =>
        c.copy(intDir = x)
      } text ("Clean directory for intermediate files. Will be used for normalization procedure")
      opt[File]('o', "output-directory") required () action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory for final output files")
      opt[Int]('w', "window-size") action { (x, c) =>
        c.copy(windowSize = x)
      } text ("OPTIONAL. Window size [nt] for insertions to search for potential alternate representations (default is 10).")
    }

    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      assume(config.outputDir.exists() && config.outputDir.isDirectory(),
        "Output directory does not exist or invalid directory!")
      //make sure large variant file exists and is a proper file
      assume(config.largeVariants.exists() && config.largeVariants.isFile(),
        "LSV file does not exist or it is an invalid file!")
      partitionVariants(config)
    }
  }

  //create intermediate output file for singletons (large variants that don't overlap with anything)
  /**Method to partition large variants file into smaller independent bins of overlapping variants*/
  def partition(remainingVariants: List[VCFLine], set_count: Int, singletons: List[VCFLine], config: Config): List[VCFLine] = {
    if (remainingVariants.isEmpty) singletons
    else {
      //get the first variant in the current list of variants
      val reference_variant = remainingVariants.head
      //check all other variants that overlap with the head variant and sort by position
      val overlappingSet = remainingVariants.filter(query => isOverlap(reference_variant, query, config.windowSize))
        .sortBy(_.pos)
      //if no overlapping set is found
      assert(overlappingSet.size != 0)
      if (overlappingSet.size == 1) {
        println("Singleton found with in chromosome " + reference_variant.refGenome + "at position " + reference_variant.pos)        
        partition(remainingVariants diff overlappingSet, set_count, singletons ::: List(overlappingSet.head), config)
      } else {
        println("Overlapping set found in chromosome " + reference_variant.refGenome + " at position " + reference_variant.pos)
        //get starting and ending coordinates for local reference sequence
        val (start, end) = {
          //make tuples with the starting and end coordinate of the variants, sort in ascending by starting 
          val tmp = overlappingSet.map(x => (x.pos, x.pos + x.refLength)).sortBy(_._1)
          //for the end coordinate, take the maximum of the second value in the tuple
          (tmp.head._1 - 10, tmp.map(_._2).max + 10)
        }
        //  println(start,end)
        //get local chromosome
        val chromosome = overlappingSet.map(_.refGenome).distinct
        assert(chromosome.size == 1, "Overlapping set contains different chromsomes:\n" + overlappingSet.mkString("\n"))
        //load reference genome as a FastaIterator
        val reference = new FastaIterator(config.reference)
        val local_reference_sequence = fastaParser(reference, start, end, chromosome.head)
        //create intermediate output file
        val pw = new PrintWriter(config.intDir + "/intermediate.overlappingSet." + set_count + ".vcf")
        pw.println(start + ";" + end + ";" + local_reference_sequence)
        pw.println(overlappingSet.map(_.toString).mkString("\n"))
        pw.close
        partition(remainingVariants diff overlappingSet, set_count + 1, singletons, config)
      }
    }
  }

  /**Method that partitions a single vcf files into independent sets of overlapping variants*/
  def partitionVariants(config: Config) {
    println("Loading large variants file")
    //load all large variants file and convert to vcf lines
    val largeVariants = tLines(config.largeVariants).map(new VCFLine(_))
    //create output file
    val singletons = new PrintWriter(config.outputDir + "/singletons.vcf")
    //Partition variants, obtain singletons, sorty by chromosome and location, redirect to output file
    singletons.println(
        partition(largeVariants, 1, List(), config)
        .sortBy(r => (r.refGenome, r.pos)).mkString("\n"))
    println("Partitioning was successfully completed!")
    singletons.close()

  }

}