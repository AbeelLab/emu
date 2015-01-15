package abeellab.emu

//10nt sequence marker for LSVs with incomplete alternative sequence
trait generalInfo {
  val incomplete_marker: String = "NNNNNNNNNN"
}

//class for analyzing vcf formatted variants
class VCFvariant(_variant: String) extends generalInfo { 
  //general information
  val variant = _variant.split("\t").toList
  val chromosome: String = variant(0).toString
  val position: Int = variant(1).toInt
  val id: AnyVal = variant(2)
  val reference_seq: String = variant(3).toString
  val alternative_seq: String = variant(4).toString
  //method that returns a boolean whether the current variant has incomplete alternative sequence
  def isIncomplete: Boolean = alternative_seq.contains(incomplete_marker)
  //computes the appropriate alternative sequence size in case the LSV has incomplete alternative sequence
  val alternative_seq_size: Int = if (isIncomplete) alternative_seq.size - incomplete_marker.size else alternative_seq.size
  //computes the absolute value of the size differences of an LSV's reference and alternative sequence size
  val difference: Int = (reference_seq.size - alternative_seq.size).abs  
  //note: removed val size = difference.abs
  //general info
  val quality: AnyVal = variant(5)
  val filter: String = variant(6).toString
  val info: String = variant(7).toString
  val format: String = variant(8).toString
  val sample: String = variant(9).toString
  //due to inconsistency in vcf format, this field will check and store any additional information in the vcf file for each variant
  val Optionalfield = {
    if(variant.size > 10) {
      variant.mkString("\t").substring(_variant.mkString.indexOf("\t"+variant(10)))
    }
    else ""
  }
  
  //note removed val emuSampleName: String and val flankingSizes
  
  //method to determine the appropriate LSV type
  def isLargeDeletion: Boolean = (reference_seq.size > 1 && alternative_seq.size == 1 && difference > 1)
  def isLargeInsertion: Boolean = (alternative_seq.size > 1 && reference_seq.size == 1 && difference > 1)
  def isLargeSubstitution: Boolean = (reference_seq.size > 1 && alternative_seq.size > 1)
  def isProperVariant: Boolean = reference_seq.matches(("[ATGCN]+")) && alternative_seq.matches(("[ATGCN]+"))
  //field for appropriate variant type
  val variant_type: String = {
    if (isLargeDeletion) "Deletion"
    else if (isLargeInsertion) "Insertion"
    else if (isLargeSubstitution) "Substitution"
    else "SNV"
  }
  val difference_true: Int = {
    if(variant_type == "Deletion" || variant_type == "Insertion") 
      (reference_seq.size - alternative_seq_size).abs
      else reference_seq.size
  }
  //convert variant to original string
  val convertString: String = _variant
}

class normalizationTools4CompleteLSVs(_variant: String, max_leniency: Double, min_comparison_size: Int,
  maxSubSize: Int) extends generalInfo {

  val variant = _variant.split("\t").toList
  val chromosome: String = variant(0).toString
  val position: Int = variant(1).toInt
  val id: AnyVal = variant(2)
  val reference_seq: String = variant(3).toString
  val alternative_seq: String = variant(4).toString
  def isIncomplete: Boolean = alternative_seq.contains(incomplete_marker)
  val alternative_seq_size: Int = if (isIncomplete) alternative_seq.size - incomplete_marker.size else alternative_seq.size
  val difference: Int = (reference_seq.size - alternative_seq.size).abs
  val size: Int = difference.abs
  val quality: AnyVal = variant(5)
  val filter: String = variant(6).toString
  val info: String = variant(7).toString

  val flankingSizes: Int = {
    val temp = info.substring(info.indexOf("FLANKING_SEQ_SIZE=") + 18)
    temp.substring(0, temp.indexOf(";")).toInt
  }
  val format: String = variant(8).toString
  val sample: String = variant(9).toString
  val Optionalfield = {
    if(variant.size > 10) {
      variant.mkString("\t").substring(_variant.mkString.indexOf("\t"+variant(10)))
    }
    else ""
  }
  val emuSampleName: String = if (!info.contains("EMU_SAMPLE_NAME=")) null else info.substring(info.indexOf("EMU_SAMPLE_NAME=") + 16)

  def trimSeqs(seq: String): String = {
    val temp = seq.substring(flankingSizes)
    temp.substring(0, temp.size - flankingSizes)
  }
  val position_true: Int = position + flankingSizes
  val reference_seq_true: String = trimSeqs(reference_seq)
  val alternative_seq_true: String = trimSeqs(alternative_seq)
  val alternative_seq_size_true: Int = {
    alternative_seq_size - (flankingSizes * 2) -
    	(
    	if(alternative_seq_true.contains(incomplete_marker)) incomplete_marker.size
    	else 0
    	)
  }
  //val difference_true: Int = (reference_seq_true.size - alternative_seq_size_true).abs

  def isLargeDeletion: Boolean = (reference_seq_true.size > 1 && alternative_seq_true.size == 1 && difference > 1)

  def isLargeInsertion: Boolean = (alternative_seq_true.size > 1 && reference_seq_true.size == 1 && difference > 1)

  def isLargeSubstitution: Boolean = (reference_seq_true.size > 1 && alternative_seq_true.size > 1)

  val variant_type: String = {
    if (isLargeDeletion) "Deletion"
    else if (isLargeInsertion) "Insertion"
    else if (isLargeSubstitution) "Substitution"
    else "SNV"
  }
   val difference_true: Int = {
    if(variant_type == "Deletion" || variant_type == "Insertion") 
      (reference_seq_true.size - alternative_seq_size_true).abs
      else reference_seq_true.size
  }

  val convertString: String = variant.mkString("\t")
  val convertString_true: String = (chromosome + "\t" + position_true + "\t" + id + "\t" + reference_seq_true+ "\t" +
        alternative_seq_true + "\t" + quality + "\t" + filter + "\t" + info + "\t" + format + "\t" + sample + Optionalfield)        
    
  val isIncomputable: Boolean =
    (variant_type == "Substitution") && (reference_seq.size + alternative_seq_size > maxSubSize)

  //find the overlaping coordinate range of two LSVs
  private def find_overlaps(that: normalizationTools4CompleteLSVs): List[Int] = {
    lazy val location_range1 = Range(position, reference_seq.size + position)
    lazy val location_range2 = Range(that.position, that.reference_seq.size + that.position)
    //the overlapping index range of their reference sequence
    val overlapping_index = location_range1.filter(location_range2.contains(_))
    //return List(-1) if their alt seq don't overlap. 
    if (overlapping_index.length == 0) List(-1)
    //normalize location range relative to the variant with the most-left genome coordinate
    else {
      val earliest = if (position <= that.position) position else that.position
      List(overlapping_index(0) - earliest,
        overlapping_index(overlapping_index.length - 1) - earliest)
    }
  }

  //method to perform reference sequence extension
  private def extend(that: normalizationTools4CompleteLSVs): List[List[String]] = {
    val range = find_overlaps(that)
    if(range.size == 1 && range(0) == -1) List()
    else{
    //method to compute the extension sequence for a given an LSV based on overlapping coordinate range
    def findExtensionSeq(lsv: normalizationTools4CompleteLSVs, isForthat: Boolean, overlap: List[Int]): List[String] = {
      val range_coord = if (isForthat) List(overlap(0) - overlap(0), overlap(1) - overlap(0)) else overlap
      val leftflank = lsv.reference_seq.substring(0, range_coord(0))
      val rightflank = if (range_coord(1) + 1 == lsv.reference_seq.size) "" else lsv.reference_seq.substring(range_coord(1) + 1)
      List(leftflank, rightflank)
    }
    //method to compute the extension sequence for a given an LSV based on overlapping coordinate range
    def findExtensionSeq2(isForthat: Boolean, overlap: List[Int]): List[String] = {
      val range_coord = if (isForthat) List(overlap(0) - overlap(0), overlap(1) - overlap(0)) else overlap
      val leftflank = reference_seq.substring(0, range_coord(0))
      val rightflank = if (range_coord(1) + 1 == reference_seq.size) "" else reference_seq.substring(range_coord(1) + 1)
      List(leftflank, rightflank)
    }
    //compute the core sequence that the two LSVs share
    val core_seq = {
      if ((0 <= range(0) && range(0) <= reference_seq.size) && (0 <= range(1) && range(1) <= reference_seq.size)) {
        reference_seq.subSequence(range(0), range(1) + 1)
      } else that.reference_seq.subSequence(range(0), range(1) + 1)
    }
    //obtain respective extension sequence
    val ext_ref_for_that = findExtensionSeq2(false, range)
    val ext_ref_for_variant = findExtensionSeq(that, true, range)
    List(ext_ref_for_that, ext_ref_for_variant)
  }
  }

  //extends the reference sequence of complete LSVs and performs a Levenshtein distance test
  private def referenceExtensionComparison4Complete(that: normalizationTools4CompleteLSVs): List[String] = {
    val extension = extend(that)
    if(extension.isEmpty) return List()
    //create a list of containing the new alternative sequence after adding the respective reference extension to each LSV
    val ext_alt = List(if (extension(0).length > 1) extension(0)(0) + that.alternative_seq + extension(0)(1)
    else if (extension(0).length == 1) extension(0)(0) + that.alternative_seq else that.alternative_seq,
      if (extension(1).length > 1) extension(1)(0) + alternative_seq + extension(1)(1)
      else if (extension(1).length == 1) extension(1)(0) + alternative_seq else alternative_seq)
    //ouput: List(0) and List(1) are compatabile
    //levenshtein(ext_alt(0), ext_alt(1)) <= compute_leniency(ext_alt(0), ext_alt(1))
    ext_alt
  }

  //method to adjust the available sequences during alternative sequence comparison for LSVs with incomplete alternative sequence
  private def tailorAltSeqs(_sequence1: String, _sequence2: String, incomplete_number: Int): List[String] = {
    //incomplete_number indicates how many incomplete LSVs are present
    if (incomplete_number == 1) {
      //split incomplete seq @ incomplete sequence marker
      val split_incomplete_seq = _sequence2.split(incomplete_marker).toList
      //tailor LSV with complete sequence to the LSV with incomplete sequence
      val complete_left_seq = _sequence1.substring(0, split_incomplete_seq(0).size)
      val complete_right_seq = _sequence1.substring(_sequence1.size - split_incomplete_seq(1).size)
      //output: list(0) and list(1) are complementary. list(2) and list(3) are complementary
      List(split_incomplete_seq(0), complete_left_seq, split_incomplete_seq(1), complete_right_seq)
    } else {
      //split incomplete seq @ incomplete sequence marker
      val sequence1 = _sequence1.split(incomplete_marker).toList
      val sequence2 = _sequence2.split(incomplete_marker).toList
      //compute the minimum size for the available sequence at the left/right of the incomplete marker in the incomplete alternative sequence
      val min_size_left_flanking_sequence = if (sequence1(0).size < sequence2(0).size) sequence1(0).size else sequence2(0).size
      val min_size_right_flanking_sequence = if (sequence1(1).size < sequence2(1).size) sequence1(1).size else sequence2(1).size
      //output: list(0) and list(1) are complementary. list(2) and list(3) are complementary
      List(sequence1(0).substring(0, min_size_left_flanking_sequence), sequence2(0).substring(0, min_size_left_flanking_sequence),
        sequence1(1).substring(sequence1(1).size - min_size_right_flanking_sequence),
        sequence2(1).substring(sequence2(1).size - min_size_right_flanking_sequence))
    }
  }

  //extends reference sequence in incomplete LSvs
  private def referenceExtensionComparison4Incomplete(that: normalizationTools4CompleteLSVs, incomplete_count: Int): List[String] = {
    val extension = extend(that)
    if(extension.isEmpty) return List()
    val ext_alt = List(if (extension(0).length > 1) extension(0)(0) + that.alternative_seq + extension(0)(1)
    else if (extension(0).length == 1) extension(0)(0) + that.alternative_seq else that.alternative_seq,
      if (extension(1).length > 1) extension(1)(0) + alternative_seq + extension(1)(1)
      else if (extension(1).length == 1) extension(1)(0) + alternative_seq else alternative_seq)

    if (incomplete_count == 1) {
      val complete = if (ext_alt(0).contains(incomplete_marker)) ext_alt(0) else ext_alt(1)
      val incomplete = if (ext_alt(0).contains(incomplete_marker)) ext_alt(1) else ext_alt(0)
      val tailored_sequences = tailorAltSeqs(complete, incomplete, incomplete_count)
      tailored_sequences
    } else {
      val tailored_sequences = tailorAltSeqs(ext_alt(0), ext_alt(1), incomplete_count)
      tailored_sequences
    }
  }

  //used x1 <= y2 && x2 <= y1.
def isOverlap(that: normalizationTools4CompleteLSVs): Boolean = {
    if (that.variant_type == "Insertion" && variant_type == "Insertion") {
      (position_true <= (that.position_true + that.alternative_seq_size_true - 1)) && 
      (that.position_true <= (position_true + alternative_seq_size_true - 1))
    } else (position_true <= that.reference_seq_true.size + that.position_true - 1) && 
            (that.position_true <= position_true + reference_seq_true.size - 1)
  }

  def check4ImprobableMerging(that: normalizationTools4CompleteLSVs): Boolean = {
    val incomplete = if (isIncomplete) alternative_seq_size else that.alternative_seq_size
    val complete = if (!isIncomplete) that.alternative_seq_size else alternative_seq_size
    incomplete < complete

  }

  //root method for extending the reference sequence to the alternative sequence depending on the types and whether they are complete or incomplete
  def readyComparison(that: normalizationTools4CompleteLSVs): List[String] = {
    val types_list = List(variant_type, that.variant_type)
    val incomplete_boolean_either = isIncomplete || that.isIncomplete
    val incomplete_boolean_both = isIncomplete && that.isIncomplete
    if (types_list.contains("Insertion") && types_list.contains("Deletion")) List()
    else if (incomplete_boolean_both) referenceExtensionComparison4Incomplete(that, 2)
    else if (incomplete_boolean_either) {
      if (check4ImprobableMerging(that)) referenceExtensionComparison4Incomplete(that, 1) else List()
    } else referenceExtensionComparison4Complete(that)
  }

}