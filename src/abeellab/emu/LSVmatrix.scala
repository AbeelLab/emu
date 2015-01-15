package abeellab.emu

import java.io.File
import atk.util.Tool
import java.io.PrintWriter

object LSVmatrix extends Tool {

  case class Config(val output: File = null, val canonicalized: File = null, val _samples: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar emu.jar metrics") {
      opt[File]('o', "output") required () action { (x, c) =>
        c.copy(output = x)
      } text ("Input working directory.")

      opt[File]('c', "canonicalized") required () action { (x, c) =>
        c.copy(canonicalized = x)
      } text ("Input file that contains canonicalized LSVs.")

      opt[File]('v', "vcf-list") required () action { (x, c) =>
        c.copy(_samples = x)
      } text ("Input file that contains list of VCF files.")
    }

    parser.parse(args, Config()).map { config =>
      assume(config.output.exists(), "The output directory does not exist!: " + config.output)
      assume(config.output.isDirectory(), "The output directory is not a directory: " + config.output)
      assume(config.canonicalized.exists(), "The input file: " + config.canonicalized + " does not exist!")
      assume(config.canonicalized.isFile(), "The input file is not a file: " + config.canonicalized)

      val samples = tLines(config._samples)
      makeMatrix(config.output.toString + "/", config.canonicalized, samples)
    }
  }
  
  //method to get name of sample.
  private def getSampleName(name: String): String =
    //assumes that the name is if preceeded by the directory full path and has an extension name (e.g. ".vcf")
    { val temp = name.substring(name.lastIndexOf("/") + 1); temp.substring(0, temp.lastIndexOf(".")) }

  private def getID(v: String): String={
      val t = v.substring(v.indexOf("ID:") + 3); t.substring(0, t.indexOf(";"))
    }
  
  def makeMatrix(outputDir: String, lsvs: java.io.File, vcf: List[String]) {        

    val pw = new PrintWriter(outputDir+"lsv_to_sample.matrix")

    println("Preparing LSVs.")
    
    val variants = tLines(lsvs)
      .filter(f => f.contains("ID:")).map(f => f.split("\t").toList(7))    

      println("Total number of variants: " + variants.size)
      
    println("Preparing samples.") 
    val strains = vcf.map(line => getSampleName(line))

    println("Total number of samples: " + strains.size)
    println("Printing LSV unique IDs.")
    pw.print("$$"); variants.foreach(f=>pw.print("\t" + getID(f)))
    
    println("Processing samples:")
    strains.foreach(f => {
      println(f)
      pw.println(); pw.print(f)
      for (variant <- variants) { if (variant.contains(f)) pw.print("\t1") else pw.print("\t0") }
    })

    pw.close

  }

}