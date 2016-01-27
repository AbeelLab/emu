package abeellab.emu

import atk.util.Tool
import java.io.File
import atk.compbio.vcf.VCFLine
import java.io.PrintWriter

object ExtractLSVs extends Tool with EmuUtils {

  case class Config(
    val vcfList: File = null,
    val min_size: Int = 2,
    val outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar emu.jar extract") {
      //command line arguments (see 'text' for description)
      opt[File]('v', "vcf-list") required () action { (x, c) =>
        c.copy(vcfList = x)
      } text ("Input text file containing paths to vcf file")

      opt[File]('s', "min-lsv-size") action { (x, c) =>
        c.copy(vcfList = x)
      } text ("OPTIONAL. Minimum size (nt) for a variant to be considered an large sequence variant (LSV). Default is 2 nt")

      opt[File]('o', "output-directory") required () action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Working output directory")
    }

    parser.parse(args, Config()).map { config =>
      assume(config.outputDir.exists() && config.outputDir.isDirectory(), "Output directory does not exist or invalid directory!")
      assume(config.vcfList.exists() && config.vcfList.isFile(), "Text file does not exist or invalid file!")
      extract_lsvs(config)
    }
  }

  /** Extract all LSVs from given list of vcf file path and outputs into a single, merged vcf file */
  def extract_lsvs(config: Config) {
    println("Extracting large variants")
    //get text file containing a list of all vcf paths
    val vcf_list = tLines(config.vcfList).map(toVCFlist(_))
    val size = vcf_list.size
    //create temporary file that will contain all extracted LSVs
    val large = new PrintWriter(config.outputDir + "/.intermediate.allLargeVariants.vcf")
    //create temporary file that will contain all extracted LSVs
    val small = new PrintWriter(config.outputDir + "/.intermediate.allSmallVariants.vcf")
    //for each file in the list, extract all LSVs and output it the file created above
    vcf_list.foreach(file => {
      println(("Processing: " + (vcf_list.indexOf(file) + 1) + "/" + size + " files"))
      //load vcf file, extract large variants, and add sample identifier to end of vcf line
      tLines(file.path).map(new VCFLine(_))
        .foreach(x => {
          if (isProperVariant(x)) {
            if (isLSV(x)) large.println(x.toString + "#&#" + file.sampleID)
            else small.println(x.toString + "#&#" + file.sampleID)
          } else small.println(x.toString + "#&#" + file.sampleID)
        })

    })
    large.close(); small.close()
/**
    println("Sorting large variants")

    val pwx = new PrintWriter(config.outputDir + "/.intermediate.allLargeVariants.vcf")
    //open temp file, sort by position, and output final file
    pwx.println(tLines(config.outputDir + "/.intermediate.tmp.allLargeVariants.vcf")
      .map(new VCFLine(_)).sortBy(_.pos).map(_.toString()).mkString("\n"))
    pwx.close

    //remove intermediate file
    val temp = new File(config.outputDir + "/.intermediate.tmp.allLargeVariants.vcf")
    temp.delete()
*/
    println("Successfully completed.")

  }

}