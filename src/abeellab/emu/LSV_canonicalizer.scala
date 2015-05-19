package abeellab.emu

import CanonicalizationMethods._
import java.io.PrintWriter
import atk.util.Tool

object LSV_canonicalizer extends Tool {

  private def determineCrossContig(parameter: Boolean, chrm1: String, chrm2: String): Boolean = {
    if (parameter) true else chrm1 == chrm2
  }

  private def isOverlap4below(reference_variant: normalizationTools4CompleteLSVs,
    LSV: normalizationTools4CompleteLSVs): Boolean = {
    reference_variant.position_true == LSV.position_true &&
      reference_variant.variant_type == LSV.variant_type
  }

  private def isIdentical(reference_variant: normalizationTools4CompleteLSVs, LSV: normalizationTools4CompleteLSVs): Boolean = {
    reference_variant.position_true == LSV.position_true &&
      reference_variant.reference_seq_true == LSV.reference_seq_true &&
      reference_variant.alternative_seq_true == LSV.alternative_seq_true
  }

  //assumes that there variants are insertions
  private def logicallyDifferent(reference_variant: normalizationTools4CompleteLSVs,
    LSV: normalizationTools4CompleteLSVs, leniency: Double): Boolean = {
    if (reference_variant.variant_type == "Insertion" &&
      LSV.variant_type == "Insertion" &&
      !reference_variant.isIncomplete &&
      !LSV.isIncomplete &&
      reference_variant.position_true == LSV.position_true) {
      val smallest = List(reference_variant.alternative_seq.size, LSV.alternative_seq.size).min
      (reference_variant.alternative_seq_size_true - LSV.alternative_seq_size_true).abs > smallest * leniency
    } else false
  }

  def canonicalizeBelowThreshold(LSVs: List[normalizationTools4CompleteLSVs],
    canonicalized_output_file: PrintWriter, boolean4LogFile: Boolean, outputDir: String,
    maximum_Substitution_size_allowed: Int, boolean_cross_contig: Boolean) {

    def recursive(remaining_variants: List[normalizationTools4CompleteLSVs]): String = {
      if (remaining_variants.isEmpty) "done"
      else {
        println("Remaining number of LSVs: " + remaining_variants.size)
        //head of the list if the reference variant
        val reference_variant = remaining_variants.head
        println("Current LSV: " + reference_variant.isIncomputable + "\t" + reference_variant.variant_type + "\t" + reference_variant.difference_true + "\t" +
          reference_variant.position_true)
        if (reference_variant.isIncomputable) {
          println("LSV is incomputable.")
          val merge = new VCFvariant(includeID(reference_variant.convertString_true))
          val idmod = merge.info + "FAILURE:SIZEOVERLOAD;"
          canonicalized_output_file.println(merge.convertString.replace(merge.info, idmod))
          recursive(remaining_variants.tail)
        } else {

          //filter list of variants that overlap in respects to the reference_variant
          val overlapping_variants = {
            remaining_variants.filter(LSV => {
              determineCrossContig(boolean_cross_contig, reference_variant.chromosome, LSV.chromosome) &&
                isOverlap4below(reference_variant, LSV)
            })
          }
          val alternate_representations = {
            overlapping_variants.filter(LSVs => {
              if (isIdentical(reference_variant, LSVs)) true
              else {
                val extendedSequences = reference_variant.readyComparison(LSVs)
                if (extendedSequences.isEmpty) false
                else computeDifferences(extendedSequences) <= computeLeniency(extendedSequences, 0)
              }
            })
          }

          val modified_list_alt_reps = {
            if (alternate_representations.size == 0) reference_variant :: alternate_representations
            else alternate_representations
          }

          val canonicalized_LSV = merge(modified_list_alt_reps)
          if (boolean4LogFile)
            log_under(canonicalized_LSV, reference_variant, overlapping_variants, outputDir)
          canonicalized_output_file.println(canonicalized_LSV.mkString("\n"))
          recursive(remaining_variants.tail diff alternate_representations)
        }
      }
    }
    recursive(LSVs)
  }

  def canonicalizeAboveThreshold(LSVs: List[normalizationTools4CompleteLSVs], max_leniency: Double,
    canonicalized_output_file: PrintWriter, boolean4LogFile: Boolean, outputDir: String,
    maximum_Substitution_size_allowed: Int, boolean_cross_contig: Boolean) {

    def recursive(remaining_variants: List[normalizationTools4CompleteLSVs]): String = {
      if (remaining_variants.isEmpty) "done"
      else {
        println("Remaining: " + remaining_variants.size)
        //head of the list if the reference variant
        val reference_variant = remaining_variants.head
        //standard output job status
        println(reference_variant.isIncomputable + "\t" + reference_variant.variant_type + "\t" + reference_variant.difference_true + "\t" +
          reference_variant.position_true)
        //if the LSV's reference or alternative sequence is above the maximum size, do not attempt to canonicalize
        if (reference_variant.isIncomputable) {
          println("LSV is incomputable.")
          val merge = new VCFvariant(includeID(reference_variant.convertString_true))
          val idmod = merge.info + "FAILURE:SIZEOVERLOAD;"
          canonicalized_output_file.println(merge.convertString.replace(merge.info, idmod))
          recursive(remaining_variants.tail)
        } else {
          //find overlaps of the reference LSV
          val overlapping_variants = {
            remaining_variants.filter(LSV => {
              if (LSV.isIncomputable) false
              //cant canonicalize intra-sample LSVs              
              else {
                determineCrossContig(boolean_cross_contig, reference_variant.chromosome, LSV.chromosome) &&
                  reference_variant.isOverlap(LSV)
              }
            })
          }
          //find true alternate representations of the overlapping LSVs
          val alternate_representations = {
            overlapping_variants.filter(LSVs => {
              if (logicallyDifferent(reference_variant, LSVs, max_leniency)) false
              else if (isIdentical(reference_variant, LSVs)) true
              else {
                val extendedSequences = reference_variant.readyComparison(LSVs)
                if (extendedSequences.isEmpty) false
                else {
                  computeDifferences(extendedSequences) <= computeLeniency(extendedSequences, max_leniency)
                }
              }
            })
          }
          //if there are no alternate reprsentations, add the reference variant to this list.
          //otherwise, add all alternate representation to this list
          val modified_list_alt_reps = if (alternate_representations.size == 0) reference_variant :: alternate_representations else alternate_representations
          val canonicalized_LSVs = merge(modified_list_alt_reps)
          if (boolean4LogFile) 
            log_over(canonicalized_LSVs, reference_variant, modified_list_alt_reps, overlapping_variants, outputDir)
          canonicalized_output_file.println(canonicalized_LSVs.mkString("\n"))
          recursive((remaining_variants diff modified_list_alt_reps).sortBy(lsv => lsv.position))
        }
      }
    }
    recursive(LSVs)
  }

  def canonicalize(_extractedLSVs: java.io.File, max_leniency: Double, minimum_comparison_length: Int,
    maximum_Substitution_size_allowed: Int, output_prefix: String,
    boolean4LogFile: Boolean, outputDir: String, boolean_cross_contig: Boolean) {

    println("Input file: " + _extractedLSVs.getAbsolutePath())
    println("Maximum percent difference allowed: " + max_leniency)
    println("Minimum comparison length: " + minimum_comparison_length + "nt")
    println("Document normalization process: " + boolean4LogFile)
    println("Maximum substitution size: " + maximum_Substitution_size_allowed + "nt")
    println("Output directory: " + outputDir)
    println("Output prefix for file: " + output_prefix)

    def isAboveThreshold(variant: normalizationTools4CompleteLSVs): Boolean = {
      if (variant.isIncomplete) variant.alternative_seq_size_true >= minimum_comparison_length
      else true
    }

    //get extracted LSV file, map each one to be an instance of VCFvariant class, and sort them by genome coordinate
    println("Mapping and sorting LSVs to appropriate class.")
    val extractedLSVs = tLines(_extractedLSVs)
      .map(lsv => new normalizationTools4CompleteLSVs(lsv, max_leniency, minimum_comparison_length,
        maximum_Substitution_size_allowed))
      .sortBy(lsv => lsv.position)
    println("Creating output file.")
    val canonicalized_output_file = new PrintWriter(outputDir + output_prefix + ".vcf")
    //separate LSVs that are above and below the
    println("Seperating LSVs above and below threshold")
    val aboveThreshold = extractedLSVs.filter(lsv => isAboveThreshold(lsv))
    val belowThreshold = (extractedLSVs diff aboveThreshold)
    println("Beginning canonicalization...")
    canonicalizeAboveThreshold(aboveThreshold, max_leniency, canonicalized_output_file, boolean4LogFile,
      outputDir, maximum_Substitution_size_allowed, boolean_cross_contig)
    canonicalizeBelowThreshold(belowThreshold, canonicalized_output_file, boolean4LogFile,
      outputDir, maximum_Substitution_size_allowed, boolean_cross_contig)
    canonicalized_output_file.close
    println("Canonicalization is complete!")

  }

}