package abeellab.emu

import atk.util.Tool
import java.io.File
import atk.compbio.vcf.VCFLine
import java.io.PrintWriter
import be.abeel.bioinformatics.FastaIterator

object Normalize extends Tool with EmuUtils {

  case class Config(
    val largeVariants: File = null,
    val windowSize: Int = 10,
    val reference: File = null,
    val intDir: File = null,
    val onlyNormalize: Boolean = false,
    val leniency: Double = 0.1,
    val outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar emu.jar normalize") {
      //command line arguments (see 'text' for description)
      opt[File]('l', "large-variants") required () action { (x, c) =>
        c.copy(largeVariants = x)
      } text ("Single vcf file containing all large variants (output of 'ExtractLSVs')")
      opt[File]('r', "reference-genome") required () action { (x, c) =>
        c.copy(reference = x)
      } text ("Reference genome that variants were called against.")
      opt[File]('i', "intermediate-directory") required () action { (x, c) =>
        c.copy(intDir = x)
      } text ("Directory for intermediate files")
      opt[File]('o', "output-directory") required () action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory for final output files")
      opt[Int]('w', "window-size") action { (x, c) =>
        c.copy(windowSize = x)
      } text ("OPTIONAL. Window size [nt] for insertions to search for potential alternate representations (default is 10).")
      opt[Double]('l', "percent-leniency") action { (x, c) =>
        c.copy(leniency = x)
      } text ("OPTIONAL. Maximum percent of mismatches allowed when substitutions or insertions are involved during comparison (default is 0.1).")
      opt[Unit]("normalization-only") action { (_, c) =>
        c.copy(onlyNormalize = true)
      } text ("Do not perform variant partioning and only perform the normalization procedure (this is turned off by default).")

    }

    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      assume(config.outputDir.exists() && config.outputDir.isDirectory(),
        "Output directory does not exist or invalid directory!")
      //make sure large variant file exists and is a proper file
      assume(config.largeVariants.exists() && config.largeVariants.isFile(),
        "Extended LSVs file does not exist or it is an invalid file!")
      val normalized_output = new PrintWriter(config.outputDir + "/intermediate.normalized.vcf")
      if (!config.onlyNormalize) partitionVariants(config, normalized_output)
      normalize(config, normalized_output)
      normalized_output.close()
      val pw = new PrintWriter(config.outputDir + "/normalized.vcf")
      pw.println(tLines(config.outputDir + "/intermediate.normalized.vcf")
        .map(new VCFLine(_)).sortBy(r => (r.refGenome, r.pos)).mkString("\n"))
      pw.close
    }
  }

  /**Main method to normalize large variants*/
  def normalize(config: Config, normalized_output: PrintWriter) {
    //get all overlapping variant files
    val intermediate_dir = config.intDir.listFiles().map(_.getAbsolutePath())
      .filter(_.contains("intermediate.overlappingSet"))
    //in parallel, iterate through each overlapping set and normalize
    intermediate_dir.foreach(file => {
      //println("Processing file: " + new java.io.File(file).getName())
      //load lines into local reference and nest case classes
      val (localReference: LocalReference, variants: List[Nest]) = {
        //load variants
        val tmp = tLines(file)
        //get the header (local reference line)
        val header = tmp.head.split(";")
        //convert header to local reference case class
        (new LocalReference(header(0).toInt, header(1).toInt, header(2).toString),
          //convert the rest of the lines into Nest case class
          tmp.tail.map(variant => {
            //for each variant, conver to EmuVCF case class
            val emuVCF = toEmuVCF(variant)
            //calculate the local position to the relative to the local reference
            val (localStart, localEnd) = (emuVCF.variant.pos - header(0).toInt,
              ((header(1).toInt - header(0).toInt)) - (header(1).toInt - (emuVCF.variant.pos + emuVCF.variant.refLength)))
            //convert to Nest case class
            new Nest(emuVCF, localStart, localEnd, 1, List(emuVCF.sampleID))
          }))
      }

      val normalized = normalizeLSVs(localReference, variants, List(), config.leniency)
      println("Overlapping set at " + variants.head.variant.variant.pos + ": " + "normalized " + variants.size +
        " variants to " + normalized.size + " variants (" + file + ")")
      normalized.foreach(x => normalized_output.println(x.variant.variant.toString + "#&#" + x.sampleIDs.mkString(";")))
    })
    println("Normalization was sucessfully completed!")
  }

  /**Method that partitions a single vcf files into independent sets of overlapping variants*/
  def partitionVariants(config: Config, normalized_output: PrintWriter) {
    println("Loading large variants file")
    //load all large variants file and convert to vcf lines
    val largeVariants = tLines(config.largeVariants).map(new VCFLine(_))
    //create intermediate output file for singletons (large variants that don't overlap with anything)
    /**Method to partition large variants file into smaller independent bins of overlapping variants*/
    def partition(remainingVariants: List[VCFLine], set_count: Int): String = {
      if (remainingVariants.isEmpty) "Exit"
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
          normalized_output.println(overlappingSet.head.toString)
          partition(remainingVariants diff overlappingSet, set_count)
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
          val pw = new PrintWriter(config.outputDir + "/intermediate.overlappingSet." + set_count + ".vcf")
          pw.println(start + ";" + end + ";" + local_reference_sequence)
          pw.println(overlappingSet.map(_.toString).mkString("\n"))
          pw.close
          partition(remainingVariants diff overlappingSet, set_count + 1)
        }
      }
    }
    //Partition variants
    partition(largeVariants, 1)
    println("Partitioning was successfully completed!")

  }
}