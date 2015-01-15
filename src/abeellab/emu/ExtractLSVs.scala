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

object ExtractLSVs extends Tool {

  case class Config(val input: File = null, val output: File = null, val _flankingSize: String = null,
    val _referenceGenome: File = null, val _LSV_size: String = null, val _LSV_prefix: String = null, 
    val _SNV_prefix: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar emu.jar extract") {
      //requires an output directory
      opt[File]('o', "output") required () action { (x, c) =>
        c.copy(output = x)
      } text ("Output directory")

      //requires a file containing a list of VCF files (one name per line)
      opt[File]('i', "input") required () action { (x, c) =>
        c.copy(input = x)
      } text ("Input file that contains a list of VCF files to analyze")

      //require a file containing list of reference genomes (one name per line) that were used for the VCF files
      opt[File]('r', "reference-genome") required () action { (x, c) =>
        c.copy(_referenceGenome = x)
      } text ("Input file that contains a list of reference genomes used in the VCF files")
      
      //optional value for the minimum size of a variant to be considered an LSV. Default is any variant greater or equal to 2nt.
      opt[String]('s', "LSV-size") optional () action { (x, c) =>
        c.copy(_LSV_size = x)
      } text ("[OPTIONAL] Input minimum size (nt) for a variant to be considered an LSV (default 2).")
      
      //optional value for the attached flanking sequences to LSVs during extraction
      opt[String]('F', "flanking-size") optional () action { (x, c) =>
        c.copy(_flankingSize = x)
      } text ("[OPTIONAL] Input flanking sequence size (nt) that will be attached to LSVs (default is 100).")

      //optional prefix name to the output file containing extracted LSVs only
      opt[String]('L', "LSV-output-prefix") optional () action { (x, c) =>
        c.copy(_LSV_prefix = x)
      } text ("[OPTIONAL] Prefix name for extacted LSV output file (default is 'AllExtractedLSVs.reference_extended'")

      //optional prefix name to the output file containing extracted SNVs only
      opt[String]('S', "SNV-output-prefix") optional () action { (x, c) =>
        c.copy(_SNV_prefix = x)
      } text ("[OPTIONAL] Prefix name for extacted SNV output file (default is  'AllExtractedSNVs'")
    }

    //checks that the file and directory exists. If the directory does not exists, it creates it.
    parser.parse(args, Config()).map { config =>
      if (config.output.exists())
        assume(config.output.isDirectory(), "The output directory is not a directory: " + config.output)
      else {
        config.output.mkdirs()
        assume(config.output.exists(), "Could not create output directory :" + config.output)
      }
      
      //checks that the input file exists
      assume(config.input.exists(), "The input does not exist!: " + config.input)
      assume(config.input.isFile(), "The input is not a file: " + config.input)

      //double checks that the desired output file names does not exist alreaady
      assume(!(new File(config.output.getAbsolutePath() + "/" + config._LSV_prefix + ".vcf")).exists,
        "The desired extracted LSV output file name already exists!")

      assume(!(new File(config.output.getAbsolutePath() + "/" + config._SNV_prefix + ".vcf")).exists,
        "The desired extracted SNV output file name already exists!")

      val LSV_prefix = if (config._LSV_prefix == null) "AllExtractedLSVs" else config._LSV_prefix

      val SNV_prefix = if (config._SNV_prefix == null) "AllExtractedSNVs" else config._SNV_prefix

      val flankingSize =
        if (config._flankingSize == null) 100
        else {
          assume(config._flankingSize.toInt.isValidInt, "Flanking sequence size is an invalid number")
          config._flankingSize.toInt
        }
      
      val LSV_size = {
        if (config._LSV_size == null) 2
        else {
          assume(config._LSV_size.toInt.isValidInt,
            "Minimum LSV size parameter is an invalid number.")
          config._LSV_size.toInt
        }
      }
      //gets the input file (list of VCFs) and list of reference genomes and transforms them into a list of strings
      val inputFiles = tLines(config.input)

      extractLSVs(inputFiles,config._referenceGenome,config.output.toString()+"/", flankingSize,
          LSV_prefix, SNV_prefix, LSV_size)
    }
  }

  //method to get name of sample.
  private def getSampleName(name: String): String =
    //assumes that the name is if preceeded by the directory full path and has an extension name (e.g. ".vcf")
    { val temp = name.substring(name.lastIndexOf("/") + 1); temp.substring(0, temp.lastIndexOf(".")) }

  def extractLSVs(vcfFiles: List[String], _referenceGenome: java.io.File, outputDir: String,
    flanking_seq_sizes: Int, LSVprefix: String, SNVprefix: String, LSV_size: Int) {
    
    println("Total of " + vcfFiles.size + " VCF file will be analyzed.")
    println("Flanking sequence size is set to " + flanking_seq_sizes + "nt.")
    println("The reference genome used is " + _referenceGenome.getAbsoluteFile())
    println("Mapping scaffolds in reference genome...")
    //get general reference genome data
    val referenceGenome = tLines(_referenceGenome)
    val referenceGenome_size = referenceGenome.size
    assume(referenceGenome(0).contains(">"), "Not a valid fasta formatted reference genome")
    //method to extract appropriate reference genomes for fasta files incase the input fasta file
    //has multiple contigs
    def concatGenome(reference: List[String], indexes: List[Int]): String = {
      if(indexes.size == 1) reference.slice(indexes(0)+1, referenceGenome_size).mkString
      else{
        reference.slice(indexes(0)+1, indexes(1)).mkString
      }
    }
    //identifies the index of the ">" marker in the fasta file which indicates the start of a contig
    val indexes = referenceGenome.zipWithIndex.filter(_._1.contains(">")).map(_._2)              
    //method to recursively create a Map with the contig name as the key and the respective
    //reference genome as a string value.
    def makeMapRecursive(indexes: List[Int], acc: Map[String,String]): Map[String,String]={
      if(indexes.isEmpty) acc
      else{
        makeMapRecursive(indexes.tail, acc + (referenceGenome(indexes.head).substring(1) ->
            concatGenome(referenceGenome, indexes)))
      }
    }
    //store map here
    val mapScaffoldtoGenome = makeMapRecursive(indexes, Map())
      

    //method to modify the VCF information after adding the flanking sequences to the reference and alternative seq.
    def extendLSV(variant: VCFvariant, flankingSeqs: List[String], emu_sample_name: String): String = {
      //add the left and right flanking sequences to the original reference seq
      val new_ref_seq = flankingSeqs(0) + variant.reference_seq + flankingSeqs(1)
      //add the left and right flanking sequences to the original alternative seq
      val new_alt_seq = flankingSeqs(0) + variant.alternative_seq + flankingSeqs(1)
      //compute the new position
      val new_position = variant.position - flankingSeqs(0).size
      //compile the full VCF string for output
      val result = (variant.chromosome + "\t" + new_position + "\t" + variant.id + "\t" + new_ref_seq + "\t" +
        new_alt_seq + "\t" + variant.quality + "\t" + variant.filter + "\t" +
        variant.info+";FLANKING_SEQ_SIZE=100;EMU_SAMPLE_NAME="+ emu_sample_name + "\t" + variant.format + "\t" + 
        variant.sample + variant.Optionalfield)
        result
    }
    
    println("Creating LSV and SNV output files.")
    //creates the output files       
    val pw_all_LSVs = new PrintWriter(outputDir + LSVprefix + ".reference_extended.vcf")
    val pw_all_SNVs = new PrintWriter(outputDir + SNVprefix + ".vcf")

    //gets reference genome and converts into one string; assumes it is in fasta file format
    def getFlankingSeqs(ref: String)(lsv: VCFvariant): List[String] = {
        val coord = lsv.position
        println("Adding flanking sequences to LSV at coordinate " + coord)
        //if left side of the starting coord is too small, take what is available
        if (coord - flanking_seq_sizes < 0) {
          assert(false, "LSV coordinate is to close to the 'starting' genome coordinate. "+
              "Please contact Alex Salazar at salazar@broadinstitute.org regarding this error.")
          val temp = ref.substring(0, coord + flanking_seq_sizes + lsv.reference_seq.size - 1)
          List(temp.substring(0, coord), temp.substring(coord + lsv.reference_seq.size))
        } //if the right side of the starting coord is too small, take what is available
        else if (coord + flanking_seq_sizes + lsv.reference_seq.size > ref.size) {
          assert(false, "LSV coordinate is to close to the 'ending' genome coordinate. "+
              "Please contact Alex Salazar at salazar@broadinstitute.org regarding this error.")
          val temp = ref.substring(coord - flanking_seq_sizes - 1)
          val normalized_coord = coord - (coord - flanking_seq_sizes)
          List(temp.substring(0, normalized_coord), temp.substring(normalized_coord + lsv.reference_seq.size))
        } //take the user input flanking sequence sizes from both ends
        else {         
          val temp = ref.substring(coord - flanking_seq_sizes - 1, coord + flanking_seq_sizes + lsv.reference_seq.size - 1)
          val normalized_coord = coord - (coord - flanking_seq_sizes)
          List(temp.substring(0, normalized_coord), temp.substring(normalized_coord + lsv.reference_seq.size))
        }
      }
    //iterates through each file in the vcfFiles list and splits up LSVs and SNVs while creating individual LSV vcf files
    vcfFiles.foreach(file => {
      val emu_sample_name = getSampleName(file)      
      println("Extracting from " + getSampleName(file) + " VCF file.")
      //val referenceGenome = tLines(_referenceGenomes).filter(line => line.startsWith(">")).flatten.mkString("\n")
      //val referenceGenome_size = referenceGenome.size

      //Not needed? val pw_sample = new PrintWriter(outputDir + LSVprefix + ".reference_extended.vcf")
      //make variant into class VCFvariant
      println("Converting variants to VCFvariant class.")
      val variants: List[VCFvariant] = tLines(file).map(variant => new VCFvariant(variant))                     														                     
      //function to compute the flanking sequence of an LSV based on it's respective reference genome 
      //(see code above) and position      
      println("Separating LSVs and SNVs")
      //determine if a variant is an LSV or SNV.
      //variant must only have the following chars in the ref and alt columns: A,T,G,C,N
      val LSVs = variants.filter(sequence_variant =>
        sequence_variant.isProperVariant && sequence_variant.difference_true >= LSV_size &&
        (sequence_variant.isLargeDeletion || sequence_variant.isLargeInsertion || 
            sequence_variant.isLargeSubstitution))       
      val SNVs = (variants diff LSVs).map(snv => snv.convertString.replace(snv.info, snv.info + 
          ";EMU_SAMPLE_NAME=" + emu_sample_name)) 
      //print all SNVs to respective output file
      SNVs.foreach(snv => pw_all_SNVs.println(snv))
      println("Outputting LSVs")
      //add the flanking sequences to LSVs and output it to its respective file
      LSVs.foreach(lsv => 
        pw_all_LSVs.println(extendLSV(lsv, getFlankingSeqs(mapScaffoldtoGenome(lsv.chromosome))(lsv), emu_sample_name)))

      //Not needed? pw_sample.print(LSVs.map(lsv => lsv.toString).mkString("\n"))
    })
    println("Closing files")
    pw_all_SNVs.close
    pw_all_LSVs.close
    println("Success: Emu extraction complete!")

  }

}