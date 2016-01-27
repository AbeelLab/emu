package abeellab.emu

import java.io.File
import atk.util.Tool
import java.io.PrintWriter
import atk.compbio.vcf.VCFLine

object VCFify extends Tool with EmuUtils {

  case class Config(
    val normalized: File = null,
    val snvs: File = null,
    val vcfList: File = null,
    val outputDir: File = null,
    val suffix: String = "normalized")

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar emu.jar vcfify") {
      //CLI arguments
      opt[File]('n', "normalized-vcf") required () action { (x, c) =>         
        c.copy(normalized = x) } text ("Input file that contains normalized variants.")
      opt[File]('s', "snv-vcf") required () action { (x, c) => 
        c.copy(snvs = x) } text ("Input file that contains single nucleotide variants")
      opt[File]('v', "vcf-list") required () action { (x, c) => 
        c.copy(vcfList = x) } text ("Input file where column 1 is the sample identifier and column 2 is the absolute path.")
      opt[File]('o', "output") required () action { (x, c) => 
        c.copy(outputDir = x) } text ("Output directory")
      opt[String]('p', "output-prefix") optional () action { (x, c) =>
        c.copy(suffix = x)
      } text ("Output suffix string for each VCF file (Default is 'normalized').")
    }

    parser.parse(args, Config()).map { config =>
      //verify that the specified output directory is a directory
      if (config.outputDir.exists())
        assume(config.outputDir.isDirectory(), "The output directory is not a directory: " + config.outputDir)
      else {
        //if it does not exist, create it
        config.outputDir.mkdirs()
        assume(config.outputDir.exists(), "Could not create output directory :" + config.outputDir)
      }
      
      //verify that each input files exists and is a proper file
      assume(config.normalized.exists(), "The input file does not exists!"+ config.normalized)
      assume(config.normalized.isFile(), "The input file is not a file: " + config.normalized)
      assume(config.snvs.exists(), "The input file does not exists!"+ config.snvs)
      assume(config.snvs.isFile(), "The input file is not a file: " + config.snvs)
      assume(config.vcfList.exists(), "The input file does not exists!"+ config.vcfList)
      assume(config.vcfList.isFile(), "The input file is not a file: " + config.vcfList)
      
      vcfify(config)
    }
  }
  
  /**Main method to create VCF files for each specified sample*/
  def vcfify(config: Config){
    println("Loading input files")
    //load file containing list of vcf files with their sample identifiers
    val vcfList = tLines(config.vcfList).map(toVCFlist(_))
    //load normalized large variants and single nucleotide variants files
    val (normalized, snvs) = (tLines(config.normalized), tLines(config.snvs))
    //iterate through each sample, and create a vcf file
    vcfList.foreach(sample => {
      println("Creating VCF file for " + sample.sampleID + "(" + 
          (vcfList.indexWhere(_.sampleID == sample.sampleID)+1) + "/" + vcfList.size + ")")
      //create output file
      val pw = new PrintWriter(config.outputDir + "/" + sample.sampleID + "." + config.suffix + ".vcf")
      //get all variants from both large normalized variants and small variants belonging to the current sample
      val (large, small) = (normalized.filter(_.contains(sample.sampleID)), snvs.filter(_.contains(sample.sampleID)))
      //output to file
      pw.println(
          //merge the large and small variants and sort by their position
          (large:::small).map(new VCFLine(_)).sortBy(r => (r.refGenome, r.pos))
          //replace the sample identifier at the end of the vcf line and output to file
          .map(x => {val index = x.toString.indexOf("#&#"); x.toString.substring(0, index)}).mkString("\n")
          )
       pw.close
    })
    println("Sucessfully completed.")
  }

}