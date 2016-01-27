package abeellab.emu

import atk.compbio.vcf.VCFLine
import be.abeel.bioinformatics.FastaIterator

/**A bundle of methods and case classes that are used throughout Emu*/
trait EmuUtils {
  val variant_type_warning = "Warning: non-long variant found. " +
    "Perhaps a single nucleotide variant sneaked in (this does not affect normalization)."
  /**
   * Case class for concatenated vcf file; uses a regular vcf line except with sample identifier value
   * and boolean value on whether the alternative sequence is complete (e.g. no "N" characters)
   */
  case class EmuVCF(variant: VCFLine, complete: Boolean, sampleID: String)
  /** Method that converts a line to EmuVCF*/
  def toEmuVCF(x: String): EmuVCF = {
    val id_index = x.indexOf("#&#")
    val id = x.substring(id_index + 3)
    val real_vcf = new VCFLine(x.substring(0, id_index))
    new EmuVCF(real_vcf, !real_vcf.alt.contains("N"), id)
  }
  /**case class for list of vcf files. This is a tab-delimited file where column 1 is the sample identifier and column 2 is the file path*/
  case class VCFlist(sampleID: String, path: String)
  /**Method to convert lines to VCFList case class*/
  def toVCFlist(x: String): VCFlist = {
    val tmp = x.split("\t")
    new VCFlist(tmp(0).toString, tmp(1).toString)
  }
  /**case class for Emu normalization data structure*/
  case class Nest(variant: EmuVCF, localStart: Int, localEnd: Int, count: Int, sampleIDs: List[String])
  /**case class for local reference genome*/
  case class LocalReference(start: Int, end: Int, sequence: String)
  /**Check if the reference variant is larger than the query*/
  def isLarger(ref: Nest, que: Nest): Boolean = {
    ref.localStart < que.localStart || ref.localEnd > que.localEnd
  }
  /** Check if the variant is a not a single nucleotide variant */
  def isLSV(x: VCFLine): Boolean = x.variation.strType.contains("Long")
  /** Check whether the variant has a proper reference and alternative vcf column entry*/
  def isProperVariant(x: VCFLine): Boolean = x.ref.matches(("[A|T|G|C|N]+")) && x.alt.matches(("[A|T|G|C|N]+"))
  /** Check if two large variants overlap in terms of their reference positions or 10 nt window if insertion */
  def isOverlap(reference: VCFLine, query: VCFLine, windowSize: Int): Boolean = {
    if (reference.refGenome != query.refGenome) false
    else {
      //obtain appropriate coordinate start and end for deletion, insertion, and substitutions
      def referenceRange: (Int, Int) = {
        reference.variation.strType match {
          //if insertion, add specified window length
          case "LongInsertion" => (reference.pos - windowSize, reference.pos + windowSize)
          //if deletion or substitution, use the variant position and size
          case _ => (reference.pos, reference.pos + reference.refLength)
        }
      }
      //get start and end reference coordinates for the reference and query variant
      val (reference_start, reference_end) = referenceRange
      val (query_start, query_end) = (query.pos, query.pos + query.refLength)
      //check if reference and query overlaps
      reference_start <= query_end && reference_end >= query_start
    }
  }
  /**Method to parse out sequence from a FASTA file*/
  def fastaParser(reference: FastaIterator, start: Int, end: Int, chrm: String): String = {
    //if there are no more chromosomes, exit
    if (!reference.hasNext()) "WARNING: No chromosome found for current set of large variants"
    else {
      //load next chromosome
      val chromosome = reference.next
      //check if the current chromosome is the one in which the large variants were identified in
      if (chromosome.getDescription().contains(chrm)) {
        //extract sequence using specified coordinates and exit
        val size = chromosome.getSequence().size
        //println(size)
        chromosome.getSequence().substring(start - 1, end)
        //else, continue searching for appropriate chromosome
      } else {
        fastaParser(reference, start, end, chrm)
      }
    }
  }

  /**Method to obtain the representative of a set of alternative representations*/
  def getRepresentative(altReps: List[Nest], sampleIDs: List[String]): Nest = {
    //method to find the representative
    def findRepresentative(set: List[Nest]): Nest = {
      //get set of breakpoint locations and their respective frequencies
      val breakpoint_frequency = getBreakpointFrequency(set)
      //make set of all the distinct breakpoint frequencies
      val distinct_frequencies = breakpoint_frequency.map(_._2).toList.distinct
      //when all the variants in the set have different breakpoint locations
      if (distinct_frequencies.size == 1 && distinct_frequencies.head == 1) {
        //get the left most variant and choose that as the representative
        val left_most = set.sortBy(_.variant.variant.pos).head
        new Nest(left_most.variant, left_most.localStart, left_most.localEnd,
          altReps.size, sampleIDs)
      } else {
        //get the maximum frequency
        val max_frequency = distinct_frequencies.sortBy(_.toInt).last
        //find the variant representative
        val representative = {
          //get breakpoint representative based on it's frequency
          val breakpoint_representative = breakpoint_frequency.filter(_._2 == max_frequency).head._1
          set.filter(_.variant.variant.pos == breakpoint_representative).head
        }
        new Nest(representative.variant, representative.localStart, representative.localEnd,
          altReps.size, sampleIDs)
      }
    }
    //method to obtain the frequency of each unique breakpoint in the set
    def getBreakpointFrequency(set: List[Nest]): Map[Int, Int] = {
      set.map(_.variant.variant.pos).groupBy(_.toInt).map(x => (x._1, x._2.size))
    }
    //find the boolean field of whether all alterantive reps are complete/incomplete
    val complete_type = altReps.map(_.variant.complete).distinct
    val subs_filtered = altReps.filter(_.variant.variant.variation.strType != "LongSubstitution")
    //if the alteranative reps are all complete
    if (complete_type.size == 1 && complete_type.head == true) {
      //if there is at least one variant that is not a substitution
      if (subs_filtered.size > 0) {
        findRepresentative(subs_filtered)
        //else, there are nothing but substitutions
      } else {
        findRepresentative(altReps)
      }
      //else, there are incomplete variants
    } else {
      //if there are variants other than substitutions
      if (subs_filtered.size > 0) {
        //get set of complete variants
        val all_complete = subs_filtered.filter(x => x.variant.complete)
        //if there is at least one variant that is complete
        if (all_complete.size > 0) findRepresentative(all_complete)
        else findRepresentative(subs_filtered)
        //else, only substitutions are present
      } else {
        //get set of complete variants
        val all_complete = altReps.filter(x => x.variant.complete)
        //if there is at least one variant that is complete
        if (all_complete.size > 0) findRepresentative(all_complete)
        else findRepresentative(altReps)
      }
    }
  }

  /**Method to check if two variants are identicals*/
  def isIdentical(reference: Nest, query: Nest): Boolean = {
    //check if the position, reference, and alternative sequence are identical
    reference.variant.variant.pos == query.variant.variant.pos &&
      reference.variant.variant.ref == query.variant.variant.ref &&
      reference.variant.variant.alt == query.variant.variant.alt
  }

  /**Method to check the true type of a substitution*/
  def getSubstitutionType(x: Nest): String = {
    //assumption 4: if the alt sequence is larger than the ref sequence OR if the alt sequence
    //is incomplete, then the substitution is an insertion. else, its a deletion
    if (x.variant.variant.altLength > x.variant.variant.refLength ||
      !x.variant.complete) "LongInsertion"
    else "LongDeletion"
  }

  /**Method to adjust incomplete insertions*/
  def adjustIncomplete(x: String, y: String, incomplete: Int): List[(String, String)] = {
    //method to get right and left flanking sequences of an incomplete sequence
    def getFlanks(z: String): (String, String) = {
      //substring up to "N" character marker
      val left_flank = z.substring(0, z.indexOf("N"))
      val right_flank = z.substring(z.indexOf("N") + 10)
      (left_flank, right_flank)
    }
    //method to adjust size when flank differ in length
    def adjustFlankSize(a: (String, String), b: (String, String)): List[(String, String)] = {
      //calculate the minimum size of each flank
      val min_lef_flank_size = List(a._1.size, b._1.size).min
      val min_right_flank_size = List(a._2.size, b._2.size).min
      //adjust each flank based on the minimum size
      List((a._1.substring(0, min_lef_flank_size), a._2.substring(a._2.size - min_right_flank_size)),
        (b._1.substring(0, min_lef_flank_size), b._2.substring(b._2.size - min_right_flank_size)))
    }
    //if both sequences are incomplete
    if (incomplete == 1) {
      val (x_flanks, y_flanks) = (getFlanks(x), getFlanks(y))
      adjustFlankSize(x_flanks, y_flanks)
      //if only one is incomplete
    } else {
      val (complete, incomplete_flanks) = ((x, x), getFlanks(y))
      adjustFlankSize(complete, incomplete_flanks)
    }
  }

  /**Method to find whether two lsvs are alternate representations*/
  def isAltRep(reference: Nest, query: Nest, localRef: LocalReference, leniency: Double): Boolean = {
    //method to replace local reference sequence with variant sequence
    def getVariantAffect(alt: String, start: Int, end: Int): String = {
      val leftFlank = localRef.sequence.substring(0, start)
      val rightFlank = localRef.sequence.substring(end)
      leftFlank + alt + rightFlank
    }
    //Method to adjust the affected sequence of if the reference variant is a substitution
    def adjustAffectedSequence(refSeq: String, queSeq: String): (String, String) = {
      //for substitution, we want to avoid other variation around the flanking regions.
      //therefore, only check the sequence that a deletion or insertion is affecting
      val start_difference = query.localStart - reference.localStart
      val end_difference = (query.localEnd - reference.localEnd).abs
      //Method to adust sequence
      def adjust(x: String): String = {
        val tmp = x.substring(start_difference)
        tmp.substring(0, tmp.size - end_difference) // -1 
      }
      //adjust the sequenced based on the start and end coordinate differences
      val (adjustedRef, adjustedQue) = (adjust(refSeq), adjust(queSeq))
      (adjustedRef, adjustedQue)
    }
    //get the effect of the variant in respects to the local reference sequence
    val reference_effect = getVariantAffect(reference.variant.variant.alt, reference.localStart, reference.localEnd)
    val query_effect = getVariantAffect(query.variant.variant.alt, query.localStart, query.localEnd)

    //Method to handle insertions
    def handleInsertionVariants(): Boolean = {
      val complete_type = List(reference.variant.complete, query.variant.complete).distinct
      //if insertions are both complete
      if (complete_type.size == 1 && complete_type.head == true) {
        levenshtein(reference_effect, query_effect) <= calculateMaxMismatches(reference_effect.size, leniency)
      } //if insertions are both incomplete
      else if (complete_type.size == 1 && complete_type.head == false) {
        //adjust the incomplete insertions
        val adjusted_insertions = adjustIncomplete(reference_effect, query_effect, complete_type.size)
        //compare each flank independently via levenshtein then add the scores
        levenshtein(adjusted_insertions.head._1, adjusted_insertions.last._1) +
          levenshtein(adjusted_insertions.head._2, adjusted_insertions.last._2) <= 
            calculateMaxMismatches(reference_effect.size, leniency)
      } else {
        //obtain the variant that is complete and incomplete
        val complete = if (reference_effect.contains("N")) query_effect else reference_effect
        val incomplete = if (complete == reference_effect) query_effect else reference_effect
        //if the incomplete variant is larger than the complete, automatically return false
        if (incomplete.size - 10 > complete.size) false
        else {
          val adjusted_insertions = adjustIncomplete(complete, incomplete, complete_type.size)
         levenshtein(adjusted_insertions.head._1, adjusted_insertions.last._1) +
            levenshtein(adjusted_insertions.head._2, adjusted_insertions.last._2) <= 
              calculateMaxMismatches(reference_effect.size, leniency)
        }
      }
    }
    //get the variant type of the reference
    reference.variant.variant.variation.strType match {
      //Use the appropriate method for alt representation detection based on variant type
      case "LongInsertion" => {
        //check the variant type of the query
        query.variant.variant.variation.strType match {
          case "LongInsertion" => {
            handleInsertionVariants()
          }
          case "LongDeletion" => false
          case "LongSubstitution" => {
            //if the substitution is a deletion-based variant, then it cannot be normalized with an insertion
            if (getSubstitutionType(query) == "LongDeletion") false
            else handleInsertionVariants()
          }
          case _ => println(variant_type_warning); false
        }
      }
      case "LongDeletion" => {
        //check the variant type of the query
        query.variant.variant.variation.strType match {
          //assumption 1: deletion and insertion are never alternate representations of one another
          case "LongInsertion" => false
          //comparisons between deletions are always exact matches
          case "LongDeletion" => reference_effect == query_effect
          case "LongSubstitution" => {
            //assumption 2: if the substitution is incomplete, then it means it represents an insertion
            //and if the substitution is complete, it represents a deletion. (see assumption 1)
            if (getSubstitutionType(query) == "LongInsertion") false
            else {
              //assumption 3: if the substitution is larger than the query, then it contains nearby variation
              //also, if the query is smaller, it truly is a multi-variant region
              if (isLarger(query, reference) && query_effect.size >= reference_effect.size) {
                //adjust the left flanking sequence accordingly
                //val (ref, que) = adjustAffectedSequence(reference_effect, query_effect)
                reference_effect.size == query_effect.size &&
                  levenshtein(reference_effect, query_effect) <= calculateMaxMismatches(reference_effect.size, leniency)
                //check whether the adjusted sequences are identical
                //ref == que
                //else, check whether the effects are identical
              } else reference_effect == query_effect
            }
          }
          case _ => println("variant_type_warning"); false
        }
      }
      case "LongSubstitution" => {
        //check the variant type of the query
        query.variant.variant.variation.strType match {
          case "LongInsertion" => handleInsertionVariants()
          case "LongDeletion" => {
            //see assumption 3
            if (isLarger(reference, query)) {
              //adjust the left flanking sequence accordingly
              // val (ref, que) = adjustAffectedSequence(reference_effect, query_effect)
              //check whether the adjusted sequences are identical
              //ref == que
              reference_effect.size == query_effect.size &&
                levenshtein(reference_effect, query_effect) <= calculateMaxMismatches(reference_effect.size, leniency)
              //else, check whether the effects are identical
            } else {
              reference_effect == query_effect
            }
          }
          case "LongSubstitution" => {
            //infer substitution type (e.g. deletion, insertion)
            val sub_types = List(getSubstitutionType(reference), getSubstitutionType(query))
            //if the substitutions size is 2, then it is a deletion and an insertion which cannot be normalized
            if (sub_types.distinct.size != 1) false
            else {
              //if the substitutions are deletions, compare via levenshtein distance
              if (sub_types.head == "LongDeletion") {
               levenshtein(reference_effect, query_effect) <= calculateMaxMismatches(reference_effect.size, leniency)
              }
              else {
                // use method to handle insertion-based variants
                handleInsertionVariants()
              }
            }
          }
          case _ => println(variant_type_warning); false

        }
      }
      case _ => println(variant_type_warning); false
    }
  }

  /**Method to normalize large variants*/
  def normalizeLSVs(localRef: LocalReference, remainingVariants: List[Nest], 
      normalized: List[Nest], leniency: Double): List[Nest] = {
    if (remainingVariants.isEmpty) normalized
    else {
      //get the head variant, treat this as the reference variant
      val ref_variant = remainingVariants.head
      //find all variants that are alternate representations
      val alternateReps = remainingVariants.filter(query => isIdentical(ref_variant, query) ||
        isAltRep(ref_variant, query, localRef, leniency))
      //get list of all distinct sample identifiers in the alternate representations
      val sampleIDs = alternateReps.map(_.variant.sampleID).distinct
      //create new normalized Nest case class
      val normalizedVariants = getRepresentative(alternateReps, sampleIDs)
      normalizeLSVs(localRef, remainingVariants diff alternateReps, normalized ::: List(normalizedVariants), leniency)
    }
  }
  /**Method to calculate leniency when substitution and insertions are involved during comparison*/
  def calculateMaxMismatches(size: Int, percent_leniency: Double): Int = (size*percent_leniency).floor.toInt

  /**
   * Levenshtein-Damareu distance algorithm.
   *  Thanks to wikipidia! (http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Scala)
   */
  def levenshtein(str1: String, str2: String): Int = {
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
    distance
  }
  def min(nums: Int*): Int = nums.min
}