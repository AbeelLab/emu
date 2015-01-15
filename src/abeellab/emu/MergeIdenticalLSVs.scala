package abeellab.emu

import atk.util.Tool
import java.io.File
import java.io.PrintWriter
import be.abeel.util.HashMap2D

object MergeIdenticalLSVs extends Tool {

  case class Config(val output: File = null, val uncanonicalized: File = null)
  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar emu.jar merge-identicals") {
      opt[File]('o', "output") required () action { (x, c) =>
        c.copy(output = x)
      } text ("Working directory")
      opt[File]('u', "uncanonicalized") required () action { (x, c) =>
        c.copy(uncanonicalized = x)
      } text ("Input data file " +
        "of extracted LSVs with no unique ID (generally AllExtractedLSVS.reference_extended.vcf).")

    }

    parser.parse(args, Config()).map { config =>
      if (config.output.exists())
        assume(config.output.isDirectory(), "The output directory is not a directory: " + config.output)
      else {
        config.output.mkdirs()
        assume(config.output.exists(), "Could not create output directory: " + config.output)
      }

      //checks that the input file exists
      assume(config.uncanonicalized.exists(), "The uncanonicalized file does not exist!: " +
        config.uncanonicalized)
      assume(config.uncanonicalized.isFile(), "The uncanonicalized file is not a file!: " +
        config.uncanonicalized)

      exportUN(config.output.toString + "/", config.uncanonicalized)
    }
  }

  def exportUN(outputDir: String, file: java.io.File) = {
    
    println("Starting merge.")
    val ID = new HashMap2D[String, String, Int]
    
    def isIdentical(reference_variant: normalizationTools4CompleteLSVs,
      LSV: normalizationTools4CompleteLSVs): Boolean = {
      reference_variant.position_true == LSV.position_true &&
        reference_variant.reference_seq_true == LSV.reference_seq_true &&
        reference_variant.alternative_seq_true == LSV.alternative_seq_true
    }

    def includeID(_variant: String): String = {
      val variant = new VCFvariant(_variant)
      val id = variant.variant_type + "." + variant.reference_seq.size + "." + variant.position +
        "." + variant.difference
      if (ID.containsKey("LSV_", id)) ID.put("LSV_", id, ID.get("LSV_", id) + 1) else ID.put("LSV_", id, 1)
      variant.convertString.replace(variant.info, variant.info + ("ID:LSV_" + ID.get("LSV_", id) + "." + id + ";"))
    }

    val pw = new PrintWriter(outputDir + "UncanonicalizedLSVs.vcf")
    
    def find_unique(lsvs: List[normalizationTools4CompleteLSVs]): String = {
      if (lsvs.isEmpty) "done"
      else {
        val reference_variant = lsvs.head
        val identical_variants = lsvs.filter(lsv => isIdentical(reference_variant, lsv))
          .sortBy(lsv => lsv.position_true)
        val left_most_variant = identical_variants(0)
        val additionalInfo_sample: String => String = representative_name =>
          identical_variants.map(lsv => lsv.emuSampleName).
            filterNot(name => name == representative_name).mkString(",")            
        val chromosome_names = identical_variants.map(lsv => lsv.chromosome).distinct
        val IDincluded = includeID(left_most_variant.convertString_true.replace(left_most_variant.info,
          left_most_variant.info + "," + additionalInfo_sample(left_most_variant.emuSampleName) + ";"))
        if (chromosome_names.size == 1) pw.println(IDincluded)
        else {
          val chromosome_names_modified = chromosome_names.mkString(".")
          pw.println(IDincluded.replace(left_most_variant.chromosome, chromosome_names_modified))
        }
        find_unique(lsvs diff identical_variants)
      }
    }
    val uncanonicalizedLSVs = tLines(file).map(lsv => new normalizationTools4CompleteLSVs(lsv,0,0,0))
    println(uncanonicalizedLSVs.size)
    find_unique(uncanonicalizedLSVs)
    pw.close
    
  }

}