package abeellab.emu

import java.io.File
import atk.util.Tool
import java.io.PrintWriter

object VCFify extends Tool {

  case class Config(val large: File = null, val output: File = null, small: File = null,
    val vcf_list: File = null, val _prefix: String = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar emu.jar vcfify") {
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output directory")
      opt[File]('c', "canonicalized") required () action { (x, c) => c.copy(large = x) } text ("Input file that contains normalized variants." +
        "This file is an output of (Normalized.emu)")
      opt[File]('s', "small") required () action { (x, c) => c.copy(small = x) } text ("Input file that contains SNPs and indels. " +
        "This file is an output of (AllExtractedSmall.emu)")
      opt[File]('v', "vcf-list") required () action { (x, c) => c.copy(vcf_list = x) } text ("Input file that contains list of VCF files. ")
      opt[String]('p', "output-prefix") optional () action { (x, c) =>
        c.copy(_prefix = x)
      } text ("Ouput prefix for each VCF file (Default is 'canonicalized').")
    }

    parser.parse(args, Config()).map { config =>
      if (config.output.exists())
        assume(config.output.isDirectory(), "The output directory is not a directory: " + config.output)
      else {
        config.output.mkdirs()
        assume(config.output.exists(), "Could not create output directory :" + config.output)
      }
      assume(config.small.exists(), "The input file does not exists!")
      assume(config.small.isFile(), "The input file is not a file: " + config.small)
      assume(config.large.exists(), "The input file does not exists!")
      assume(config.large.isFile(), "The input file is not a file: " + config.large)
      assume(config.vcf_list.exists(), "The input file does not exists!")
      assume(config.vcf_list.isFile(), "The input file is not a file: " + config.vcf_list)

      val prefix: String ={
        if (config._prefix == null) ".canonicalized"
        else 
          "." + config._prefix.toString
      }
          

      Emu2VCF(config.output.toString() + "/", config.large, config.small, config.vcf_list, prefix)
    }
  }

  def Emu2VCF(outputDir: String, LSVs: java.io.File, small: java.io.File, _sample_list: java.io.File,
      prefix: String) {
    def getSampleName(name: String): String =
      //assumes that the name is if preceeded by the directory full path and has an extension name (e.g. ".vcf")
      { val temp = name.substring(name.lastIndexOf("/") + 1); temp.substring(0, temp.lastIndexOf(".")) }
    println("Getting sample names.")
    val sample_list = tLines(_sample_list)
    println("Getting variants.")
    val large_variants = tLines(LSVs)
    val small_variants = tLines(small)
    println("Converting Emu to VCF format:")
    sample_list.foreach(sample => {
      val sample_names = getSampleName(sample)
      println(sample_names)
      val large_variants_subset = large_variants.filter(lsv => lsv.contains(sample_names))
      val small_variants_subset = small_variants.filter(lsv => lsv.contains(sample_names))
      val concatted_variants = (large_variants_subset ::: small_variants_subset)
        .sortBy(lsv => lsv.split("\t").toList(1).toInt)
      val pw = new PrintWriter(outputDir + sample_names + prefix + ".vcf")
      pw.println("#CHROM" + "\t" + "POS" + "\t" + "ID" + "\t" + "REF" + "\t" +
        "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" +
        "FORMAT" + "\t" + "SAMPLE")
      pw.println(concatted_variants.mkString("\n"))
      pw.close
    })
    println("VCFify complete!")
  }

}