package abeellab.emu

import atk.util.Tool
import java.io.PrintWriter
import java.io.File

/**
 * Given an input directory of .vcf files, it will iterate through every .vcf file and extract all PASS-ONLY
 * LSVs (deletion, insertion, or substitution). It will then output corresponding files for each .vcf file
 * which contains all LSVs within the respective .vcf file. It will also output a file containing
 * all the LSVs identified from the input directory.
 */

object Canonicalize extends Tool {

  case class Config(val input: File = null, val output: File = null, val _max_leniency: String = null,
    val _minimum_comparison_length: String = null, val _output_prefix: String = null,
    val documentLogs: String = null, val _maxSubSize: String = null, val _cross_contig: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar emu.jar canonicalize") {
      //requires an output directory
      opt[File]('o', "output") required () action { (x, c) =>
        c.copy(output = x)
      } text ("Output directory.")

      //requires an LSV-only VCF file (generally AllExtractedLSVs_reference_extended.vcf)
      opt[File]('i', "input") required () action { (x, c) =>
        c.copy(input = x)
      } text ("Input file extracted LSVs only.")

      //optional a value for the maximum percent difference allowed (e.g., for five percent difference, input "-l 0.05") 
      opt[String]('l', "max_leniency") optional () action { (x, c) =>
        c.copy(_max_leniency = x)
      } text ("[OPTIONAL] Input maximum percent difference (default is 0.01).")

      //optional a value for the minimum context length for LSVs with incomplete alternatie sequence 
      opt[String]('c', "min-comparison-length") optional () action { (x, c) =>
        c.copy(_minimum_comparison_length = x)
      } text ("[OPTIONAL] Input minimum comparison length (default is 20).")

      //optional an prefix name to the output file containing the cananocalized LSVs
      opt[String]('p', "output-prefix") optional () action { (x, c) =>
        c.copy(_output_prefix = x)
      } text ("[OPTIONAL] Prefix name for canonicalized LSV output file (default is 'CanonicalizedLSVs').")

      //optional a string of 'true' or 'false'
      opt[String]('d', "document-logs") optional () action { (x, c) =>
        c.copy(documentLogs = x)
      } text ("[OPTIONAL] Turn on/off documentation of normalization process (default is 'true').")

      //optional int of maximum substitution size allowed for normalization
      opt[String]('s', "max-sub-size") optional () action { (x, c) =>
        c.copy(documentLogs = x)
      } text ("[OPTIONAL] Largest substitution size allowed for normalization (default is 14999).")

      //optional parameter to allow cross-contig canonicalization
      opt[String]('x', "cross-contig") optional () action { (x, c) =>
        c.copy(_cross_contig = x)
      } text ("[OPTIONAL] Boolean whether to allow cross contig canonicalization of LSVs (default is false).")

    }

    parser.parse(args, Config()).map { config =>
      //checks that the file and directory exists. If the directory does not exists, it creates it.
      if (config.output.exists())
        assume(config.output.isDirectory(), "The output directory is not a directory: " + config.output)
      else {
        config.output.mkdirs()
        assume(config.output.exists(), "Could not create output directory :" + config.output)
      }

      //checks that the input file exists
      assume(config.input.exists(), "The input does not exist!: " + config.input)
      assume(config.input.isFile(), "The input is not a file: " + config.input)

      val documentation_directory = new File(config.output.toString() + "/normalization_log")
      if (documentation_directory.exists())
        throw new IllegalStateException("A '/normalization_log' directory already exists " +
          "inside the specified output directory! ")
      else {
        documentation_directory.mkdirs()
        assume(documentation_directory.exists(), "Could not create output directory :" +
          documentation_directory.getAbsoluteFile())
      }

      //double checks that the desired output file names does not exist already
      assume(!(new File(config.output.getAbsolutePath() + "/" + config._output_prefix + ".vcf")).exists,
        "The desired extracted LSV output file name already exists!")

      //if no output prefix filename or max percent leniency is given, it uses default values
      val output_prefix = if (config._output_prefix == null) "CanonicalizedLSVs" else config._output_prefix

      val max_leniency: Double =
        if (config._max_leniency == null) 0.01
        else {
          assume(config._max_leniency.toDouble.isValidInt,
            "Max difference leniency for comparison parameter is an invalid number.")
          config._max_leniency.toDouble
        }

      val minimum_comparison_length: Int =
        if (config._minimum_comparison_length == null) 20
        else {
          assume(config._max_leniency.toInt.isValidInt,
            "Minimum comparison length for incomplete LSVs is an invalid number")
          config._minimum_comparison_length.toInt
        }

      val boolean4LogFile = {
        if (config.documentLogs == null) true
        else {
          config.documentLogs match {
            case "false" => false
            case "true" => true
            case _ => throw new IllegalArgumentException("Unrecognized boolean field for 'document-logs' parameter." +
              "If used, enter 'true' or 'false'")
          }
        }
      }

      val maxSubSize = {
        if (config._maxSubSize == null) 14999
        else {
          assume(config._max_leniency.toInt.isValidInt,
            "Maxmimum substitution size for normalization parameter is an invalid number.")
          config._maxSubSize.toInt
        }
      }

      val cross_contig = {
        if (config._cross_contig == null) false
        else {
          config._cross_contig match {
            case "false" => false
            case "true" => true
            case _ => throw new IllegalArgumentException("Unrecognized boolean field for 'cross-contig' parameter." +
              "If parameter is in used, enter 'true' or 'false'")
          }
        }
      }

      LSV_canonicalizer.canonicalize(config.input, max_leniency, minimum_comparison_length,
        maxSubSize, output_prefix, boolean4LogFile, config.output.toString() + "/", cross_contig)
    }

  }
}