package abeellab.emu


import atk.util.Tool
import java.io.PrintWriter
import java.io.File

object FASTA2VCF extends Tool {

  case class Config(val fastaList: File = null, val output: File = null, val contig: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar emu.jar fasta2VCF") {
      //requires an output directory
      opt[File]('o', "output") required () action { (x, c) =>
        c.copy(output = x)
      } text ("Output directory")

      opt[File]('i', "fastaList") required () action { (x, c) =>
        c.copy(fastaList = x)
      } text ("FASTA-formatted variants")        
      
      opt[String]('c', "contig-name") required () action { (x, c) =>
        c.copy(contig = x)
      } text ("Name of contig/chromosome")  

    }

    parser.parse(args, Config()).map { config =>
      if (config.output.exists())
        assume(config.output.isDirectory(), "The output directory is not a directory: " + config.output)
      else {
        config.output.mkdirs()
        assume(config.output.exists(), "Could not create output directory :" + config.output)
      }

      //checks that the input file exists
      assume(config.fastaList.exists(), "The input does not exist!: " + config.fastaList)
      assume(config.fastaList.isFile(), "The input is not a file: " + config.fastaList)
      
      convertFASTA2VCF(config.fastaList, config.output.toString+"/", config.contig)

    }
  }
  
  def convertFASTA2VCF(fasta_list: java.io.File, outputDir: String, contig: String) {
    //assumes <ProjectName>.<SampleName>.insertions.fa (e.g., EmuBenchmark.SRR0001.insertions.fa
    def getSampleName(name: String): String ={      
    	name.substring(name.indexOf(".")+1, name.lastIndexOf("insertions")-1)
    }
    
    def formatVariant(header: String, alt_seq: String): String ={
      val ref_seq = alt_seq(0)
      val coordinate = {
        val temp = header.substring(header.indexOf("pos")+4, header.size-1)
        if(temp.contains("HET")) temp.substring(0, temp.indexOf("_")) else temp
      }
      List(contig, coordinate, ".", ref_seq, alt_seq, ".", "PASS", "CONVERTED_FROM_FASTA", "GT", "1/1")
      	.mkString("\t")            
    }

    def formatData(data: List[String], pw: PrintWriter) {
      val keys = data.filter(x => x.startsWith(">"))
      val values = data diff keys
      assert(keys.size == values.size, "Size of keys and values do not match!")
      for (index <- 0 to (keys.size - 1)) {
        val vcf_variant = formatVariant(keys(index), values(index))
        pw.println(vcf_variant)
      }
    }

    val sampleName: String => String = name =>{ 
      name.substring(name.indexOf(".")+1, name.lastIndexOf("insertions")-1)
    }
    println("Getting fasta files.")
    val list = tLines(fasta_list)
    list.foreach(file => {
      println("Processing " + sampleName(file))
      val pw_temp = new PrintWriter(outputDir + sampleName(file) + ".emuConverted.vcf")
      formatData(tLines(file), pw_temp)
      pw_temp.close
    })

  }
}
