package test

import org.scalatest.FunSuite
import atk.util.Tool
import atk.compbio.vcf.VCFLine
import abeellab.emu.ExtractLSVs
import java.io.File

class ExtractionTest extends FunSuite with Tool {
  //arguments for testing the extraction script for Emu
  val arguments = Array("-v", "test_data/test_vcf/allVCFs.txt",
    "-o", "test_data/extraction")
  //first test checks whether the script runs to completion
  test("Process VCF files to completion") {
    ExtractLSVs.main(arguments)
    println("finished")
  }

  //open vcf files
  val lsvs = tLines("test_data/extraction/allLargeVariants.vcf").map(new VCFLine(_))
  val snvs = tLines("test_data/extraction/allSmallVariants.vcf").map(new VCFLine(_))

  //second test checks whether the output files are accurate
  test("Validate number of LSVs and SNVs extracted") {
    assert(lsvs.size == 91, "Unexpected number of LSVs extracted")
    assert(snvs.size == 319, "Unexpected number of SNVs extracted")
    assert(lsvs.size + snvs.size == 410, "Total number of extracted variants do not add up to original amount in VCF files")
  }
  
  //third test manually checks that all lsvs extracted are actually lsvs
  test("Manually validate lsvs by checking ref and alt sizes"){
    val lsvs_ref_alt = lsvs.map(x => (x.ref, x.alt)).filter(y => {
      (y._1.size > 1|| y._2.size > 1)
    })
    assert(lsvs_ref_alt.size == lsvs.size, "Variants in LSVs output file are not all large variants")
  }

}