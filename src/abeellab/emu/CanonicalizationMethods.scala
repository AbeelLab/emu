package abeellab.emu

import atk.util.Tool
import java.io.PrintWriter
import java.io.File
import be.abeel.util.HashMap2D

object CanonicalizationMethods extends Tool {

  private val cache = new HashMap2D[String, String, Int]

  private val ID = new HashMap2D[String, String, Int]

  private def leniency_evaluator2(str1: String, str2: String, leniency: Double): Double = {
    val smallest_size = List(str1.size, str2.size).min
    smallest_size match {
      case x if (x < 10) => 0
      case x if (x > 9 && x < 101) => (leniency * 100)
      case x if (x > 101) => (smallest_size * leniency)
    }
  }

  private def leniency_evaluator4(str1: String, str2: String, str3: String, str4: String, leniency: Double): Double = {
    val smallest_size = List(str1.size + str2.size, str3.size + str4.size).min
    smallest_size match {
      case x if (x < 10) => 0
      case x if (x > 9 && x < 101) => leniency * 100
      case x if (x > 101) => (smallest_size * leniency)
    }
  }

  def includeID(_variant: String): String = {
    val variant = new VCFvariant(_variant)
    val id = variant.variant_type + "." + variant.reference_seq.size + "." + variant.position +
      "." + variant.difference

    if (ID.containsKey("LSV_", id)) ID.put("LSV_", id, ID.get("LSV_", id) + 1) else ID.put("LSV_", id, 1)

    variant.convertString.replace(variant.info, variant.info + ("ID:LSV_" + ID.get("LSV_", id) + "." + id + ";"))
  }

  def log_over(merge: String, reference_variant: normalizationTools4CompleteLSVs,
    alt_reps: List[normalizationTools4CompleteLSVs], overlaps: List[normalizationTools4CompleteLSVs],
    outputDir: String) {
    val normalized = new VCFvariant(merge)
    val name = normalized.info.substring(normalized.info.indexOf("ID:") + 3, normalized.info.lastIndexOf(";"))
    val sup = new PrintWriter(outputDir + "normalization_log/" + name + ".txt")
    sup.println("Merged variant:\n" + merge)
    sup.println("Reference variant:\n" + reference_variant.convertString_true)
    sup.println("Above alternative sequence threshold:\ttrue")
    sup.println("Overlapping variants count:\t" + overlaps.size)
    sup.println("Alternate representations count:\t" + alt_reps.size)
    sup.println("Alternate representations:")
    alt_reps.foreach(f => sup.println(f.convertString_true))
    sup.println("Non-alternate representations:")
    (overlaps diff alt_reps).foreach(f => sup.println(f.convertString_true))
    sup.close()
  }

  def log_under(merge: String, reference_variant: normalizationTools4CompleteLSVs,
    identical_variants: List[normalizationTools4CompleteLSVs], outputDir: String) {

    val normalized = new VCFvariant(merge)
    val dir = new java.io.File(outputDir + "normalization_log/").listFiles.toList.map { f => f.getName() }
    val name = normalized.info.substring(normalized.info.indexOf("ID:") + 3, normalized.info.lastIndexOf(";"))
    val sup = new PrintWriter(outputDir + "normalization_log/" + name)
    sup.println("Merged variant:\t" + merge)
    sup.println("Reference variant:\t" + reference_variant.convertString_true)
    sup.println("Above alternative sequence threshold:\tfalse")
    sup.println("Alternate representations count:\t" + identical_variants.size)
    sup.println("Alternate representations:\t" + identical_variants.size)
    sup.println("Alternate representations variants:")
    identical_variants.foreach(f => sup.println(f.convertString_true))
    sup.close()
  }

  private def mergeSameLSVtype(alternate_representations: List[normalizationTools4CompleteLSVs], LSV_type: String): String = {
    val additionalInfo_sample: String => String = representative_name =>{
      val temp = alternate_representations.map(lsv => lsv.emuSampleName)
                 .filterNot(name => name == representative_name).mkString(",")
       if(temp.size == 0) temp else ","+temp
    }
    val representative_LSV: normalizationTools4CompleteLSVs = {
      if (LSV_type == "Deletion") {
        alternate_representations.sortBy(lsv => lsv.reference_seq_true.size).reverse(0)
      } else if (LSV_type == "Insertion") {
        alternate_representations.sortBy(lsv => lsv.alternative_seq_size_true).reverse(0)
      } else {
        alternate_representations.sortBy(lsv => lsv.difference_true).reverse(0)
      }
    }
      val chromosome_names = alternate_representations.map(lsv => lsv.chromosome).distinct
      val IDincluded = includeID(representative_LSV.convertString_true.replace(representative_LSV.info,
    		  		representative_LSV.info + "," + additionalInfo_sample(representative_LSV.emuSampleName)+";"))
      if(chromosome_names.size == 1) IDincluded
      else{
        val chromosome_names_modified = chromosome_names.mkString(".")
        IDincluded.replace(representative_LSV.chromosome, chromosome_names_modified)
      }
  }

  private def getMostParsimoniousLSV(alternate_representations: List[normalizationTools4CompleteLSVs], LSV_type: String): String = {
    val additionalInfo_sample: String => String = representative_name =>{
      val temp = alternate_representations.map(lsv => lsv.emuSampleName)
                 .filterNot(name => name == representative_name).mkString(",")
       if(temp.size == 0) temp else ","+temp
    }
    val representative_LSV: normalizationTools4CompleteLSVs = {
      if (LSV_type == "Deletion.Substitution") {
        alternate_representations.filter(lsv => lsv.variant_type == "Deletion").sortBy(lsv => lsv.reference_seq_true.size).reverse(0)
      } else {
        assert(LSV_type == "Insertion.Substitution", "Expected Insertion and Substitution mix type in merging")
        alternate_representations.filter(lsv => lsv.variant_type == "Insertion").sortBy(lsv => lsv.alternative_seq_size_true).reverse(0)
      }
    }
    val chromosome_names = alternate_representations.map(lsv => lsv.chromosome).distinct
    val IDincluded = includeID(representative_LSV.convertString_true.replace(representative_LSV.info,
      representative_LSV.info + additionalInfo_sample(representative_LSV.emuSampleName) + ";"))
    if(chromosome_names.size == 1) IDincluded
      else{
        val chromosome_names_modified = chromosome_names.mkString(".")
        IDincluded.replace(representative_LSV.chromosome, chromosome_names_modified)
      }  
  }

  def merge(alternate_representations: List[normalizationTools4CompleteLSVs]): String = {
    val types = alternate_representations.map(lsv => lsv.variant_type).distinct
    if (types.size == 1) {
      types(0) match {
        case "Deletion" => mergeSameLSVtype(alternate_representations, "Deletion")
        case "Insertion" => mergeSameLSVtype(alternate_representations, "Insertion")
        case "Substitution" => mergeSameLSVtype(alternate_representations, "Substitution")
      }
    } else {
      assert(!(types.contains("Insertion") && types.contains("Deletion")), "Insertion and deletion present in merging")
      if (types.contains("Deletion") && types.contains("Substitution")) getMostParsimoniousLSV(alternate_representations, "Deletion.Substitution")
      else {
        assert(types.contains("Insertion") && types.contains("Substitution"), "Expected Insertion and Substitution mix type in merging")
        getMostParsimoniousLSV(alternate_representations, "Insertion.Substitution")
      }
    }
  }

  def computeLeniency(sequencesList: List[String], leniency: Double): Double = {
    assert(sequencesList.size == 2 || sequencesList.size == 4, "List size of sequences two compare is not 2 or 4")
    sequencesList.size match {
      case 2 => leniency_evaluator2(sequencesList(0), sequencesList(1), leniency)
      case 4 => leniency_evaluator4(sequencesList(0), sequencesList(1), sequencesList(2), sequencesList(3), leniency)
    }
  }

  //add levenshtein distance and compute leniency independently
  /**
   * \Thanks to wikipidia! (http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Scala)
   */
  private def levenshtein(str1: String, str2: String): Int = {
    if (str1.equals(str2)) 0
    else if (cache.containsKey(str1, str2)) cache.get(str1, str2)
    else {
      val lenStr1 = str1.length
      val lenStr2 = str2.length

      val d: Array[Array[Int]] = Array.ofDim(lenStr1 + 1, lenStr2 + 1)

      for (i <- 0 to lenStr1) d(i)(0) = i
      for (j <- 0 to lenStr2) d(0)(j) = j

      for (i <- 1 to lenStr1; j <- 1 to lenStr2) {
        val cost = if (str1(i - 1) == str2(j - 1)) 0 else 1

        d(i)(j) = min(
          d(i - 1)(j) + 1, // deletion
          d(i)(j - 1) + 1, // insertions
          d(i - 1)(j - 1) + cost // substitution
          )

      }
      val distance = d(lenStr1)(lenStr2)
      cache.put(str1, str2, distance)
      distance
    }
  }

  private def min(nums: Int*): Int = nums.min

  def computeDifferences(sequencesList: List[String]): Int = {
    assert(sequencesList.size == 2 || sequencesList.size == 4, "List size of sequences two compare is not 2 or 4")
    sequencesList.size match {
      case 2 => levenshtein(sequencesList(0), sequencesList(1))
      case 4 => levenshtein(sequencesList(0), sequencesList(1)) + levenshtein(sequencesList(2), sequencesList(3))
    }
  }

}