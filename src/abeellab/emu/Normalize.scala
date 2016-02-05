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
      opt[File]('i', "intermediate-directory") required () action { (x, c) =>
        c.copy(intDir = x)
      } text ("Directory for intermediate files")
      opt[File]('o', "output-directory") required () action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory for final output files")
      opt[Double]('l', "percent-leniency") action { (x, c) =>
        c.copy(leniency = x)
      } text ("OPTIONAL. Maximum percent of mismatches allowed when substitutions or insertions are involved during comparison (default is 0.1).")
    }

    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      assume(config.outputDir.exists() && config.outputDir.isDirectory(),
        "Output directory does not exist or invalid directory!")
      //make sure large variant file exists and is a proper file
      assume(config.largeVariants.exists() && config.largeVariants.isFile(),
        "Extended LSVs file does not exist or it is an invalid file!")
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
      val (localReference: LocalReference, variants: List[Nest]) = toLocalReferenceAndNest(file)
      //runs the normalization methods
      val normalized = normalizeLSVs(localReference, variants, List(), config.leniency)
      println("Overlapping set at " + variants.head.variant.variant.pos + ": " + "normalized " + variants.size +
        " variants to " + normalized.size + " variants (" + file + ")")
      normalized.foreach(x => normalized_output.println(x.variant.variant.toString + "#&#" + x.sampleIDs.mkString(";")))
    })
    println("Normalization was sucessfully completed!")
  }
}