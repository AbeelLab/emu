package test

import org.scalatest.FunSuite
import atk.util.Tool
import abeellab.emu.VariantOverlapPartioner
import java.io.File
import abeellab.emu.EmuUtils
import abeellab.emu.Normalize

class NormalizationTest extends FunSuite with EmuUtils {

  val default_leniency = 0.1
  //get all files in intermediate directory
  val dir = new File("test_data/normalization/")
    .listFiles().filter(_.getName().contains("Set"))

  dir.foreach(file => {
    //get file name
    val filename = file.getName()
    //get file id (number of allocated to the end of file name)
    val fileid = filename.substring(filename.indexOf("Set") + 4, filename.indexOf(".vcf")).toInt
    //extract the local reference and variants
    val (localRef: LocalReference, variants: List[Nest]) = toLocalReferenceAndNest(file.getAbsolutePath())
    val refeffect = getVariantAffect(localRef, variants(0).variant.variant.alt, variants(0).localStart, variants(0).localEnd)
    val queryeffect = getVariantAffect(localRef, variants(1).variant.variant.alt, variants(1).localStart, variants(1).localEnd)
    //manually check that the results of the normalization of input variants is correct
    fileid match {
      case 8 => {
        //make sure that the normalization procedure runs to completion without any issues
        test("File 8: Normalization runs to completion") {
          normalizeLSVs(localRef, variants, List(), default_leniency)
        }
        //check that the normalized version is correct/expected
        val normalized_version = normalizeLSVs(localRef, variants, List(), default_leniency)
        test("Checking the normalization results of two substitutions (inferred complete insertions) and two complete insertions") {
          assert(normalized_version(0).variant.variant.ref == variants(3).variant.variant.ref &&
            normalized_version(0).variant.variant.alt == variants(3).variant.variant.alt &&
            normalized_version(0).variant.variant.pos == variants(3).variant.variant.pos &&
            normalized_version(0).sampleIDs.size == 4 &&
            normalized_version(0).sampleIDs.exists(x => x == "vcf1" || x == "vcf2" || x == "vcf3" || x == "vcf4"),
            "Unexpected version of the normalized variant")
        }
      }
      case 9 => {
        //make sure that the normalization procedure runs to completion without any issues
        test("File 9: Normalization runs to completion") {
          normalizeLSVs(localRef, variants, List(), default_leniency)
        }
        //check that the normalized version is correct/expected
        val normalized_version = normalizeLSVs(localRef, variants, List(), default_leniency)
        assert(normalized_version(0).variant.variant.ref == variants(2).variant.variant.ref &&
          normalized_version(0).variant.variant.alt == variants(2).variant.variant.alt &&
          normalized_version(0).variant.variant.pos == variants(2).variant.variant.pos &&
          normalized_version(0).sampleIDs.size == 3 &&
          normalized_version(0).sampleIDs.exists(x => x == "vcf1" || x == "vcf2" || x == "vcf3"),
          "Unexpected version of the normalized variant")
      }

      case 11 => {
        //there are two methods to handle identical variants The first one is by default. The second one is a sanity self check
        test("Testing identical insertions") {
          assert(isIdentical(variants(0), variants(1)),
            "Two identical insertions were not merged with the 'isIdentical' method")
          assert(isAltRep(variants(0), variants(1), localRef, default_leniency),
            "Two identical insertions were not merged with the 'isAltRep' method")
        }
        //make sure that the normalization procedure runs to completion without any issues
        test("File 11: Normalization runs to completion") {
          normalizeLSVs(localRef, variants, List(), default_leniency)
        }
        //get the normalized version of the two variants
        val normalized_version = normalizeLSVs(localRef, variants, List(), default_leniency)
        //check that the normalized version is correct/expected
        test("File 11: Normalized version is correct") {
          assert(normalized_version(0).variant.variant.ref == variants(1).variant.variant.ref &&
            normalized_version(0).variant.variant.alt == variants(1).variant.variant.alt &&
            normalized_version(0).variant.variant.pos == variants(1).variant.variant.pos &&
            normalized_version(0).sampleIDs.size == 2 &&
            normalized_version(0).sampleIDs.exists(x => x == "vcf1" || x == "vcf2"),
            "Unexpected version of the normalized variant")
        }
      }
      case 12 => {
        test("Testing an incomplete substitution and a complete insertion") {
          assert(!isAltRep(variants(0), variants(1), localRef, default_leniency),
            "A different incomplete substitutions and a complete insertion were merged into one")
        }
      }
      case 14 => {
        test("Testing two insertion substitutions with incomplete sequence") {
          assert(!isAltRep(variants(0), variants(1), localRef, default_leniency),
            "Two different substitutions (of insertion type) were merged into one")
        }
      }

      case 15 => {
        test("Testing an insertion and a substitution") {
          assert(isAltRep(variants(0), variants(1), localRef, default_leniency),
            "A substitution and an insertion with the same effect were not merged")
        }
        //make sure that the normalization procedure runs to completion without any issues
        test("File 15: Normalization runs to completion") {
          normalizeLSVs(localRef, variants, List(), default_leniency)
        }
        //get the normalized version of the two variants
        val normalized_version = normalizeLSVs(localRef, variants, List(), default_leniency)
        //check that the normalized version is correct/expected
        test("File 15: Normalized version is correct") {
          assert(normalized_version(0).variant.variant.ref == variants(1).variant.variant.ref &&
            normalized_version(0).variant.variant.alt == variants(1).variant.variant.alt &&
            normalized_version(0).variant.variant.pos == variants(1).variant.variant.pos &&
            normalized_version(0).sampleIDs.size == 2 &&
            normalized_version(0).sampleIDs.exists(x => x == "vcf1" || x == "vcf2"),
            "Unexpected version of the normalized variant")
        }
      }
      case _ => "The file is not in use for testing"
    }
  })

}