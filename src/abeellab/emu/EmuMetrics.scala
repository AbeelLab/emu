package abeellab.emu

import java.io.File
import atk.util.Tool
import java.io.PrintWriter

object EmuMetrics extends Tool {

  case class Config(val output: File = null, val extracted: File = null, val canonicalized: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar emu.jar metrics") {
      opt[File]('o', "output") required () action { (x, c) =>
        c.copy(output = x)
      } text ("Input working directory.")

      opt[File]('u', "uncanonicalized") required () action { (x, c) =>
        c.copy(extracted = x)
      } text ("Input file that contains a uncanonicalized LSVs.")

      opt[File]('c', "canonicalized") required () action { (x, c) =>
        c.copy(canonicalized = x)
      } text ("Input file that contains canonicalized LSVs.")
    }

    parser.parse(args, Config()).map { config =>
      assume(config.output.exists(), "The working directory does not exist!")
      assume(config.output.isDirectory(), "The output directory is not a directory: " + config.output)
      for (i <- List(config.extracted, config.canonicalized)) {
        assume(i.exists(), "The input file: " + i + " does not exist!")
        assume(i.isFile(), "The input file is not a file: " + i)
      }

      //val inputFiles = tLines(config.input)
      emuMetrics(config.output.toString + "/", config.extracted, config.canonicalized)
    }
  }

  def emuMetrics(outputDir: String, _uncanonicalized: java.io.File, _canonicalized: java.io.File) {
    
    def determineBoolean(lsv: VCFvariant, status: String): Boolean ={
      status match{
        case "Incomplete" => lsv.isIncomplete
        case "Complete" => !lsv.isIncomplete
        case _ => true
      } 
    }
   
    def getGeneralTotal(lsv: List[VCFvariant]): String ={
      val total = lsv.size
	    val deletions = lsv.filter(x => x.isLargeDeletion).size
	    val insertions = lsv.filter(x => x.isLargeInsertion).size
	    val substitutions = lsv.filter(x => x.isLargeSubstitution).size
	    List("Total", total, deletions, insertions, substitutions).mkString("\t")
    }
   
	  def getGeneralMetrics(lsv: List[VCFvariant], status: String): String={	     
	    val total = lsv.filter(x => determineBoolean(x, status)).size
	    val deletions = lsv.filter(x => x.isLargeDeletion && determineBoolean(x, status)).size
	    val insertions = lsv.filter(x => x.isLargeInsertion && determineBoolean(x, status)).size
	    val substitutions = lsv.filter(x => x.isLargeSubstitution && determineBoolean(x, status)).size
	    List(status, total, deletions, insertions, substitutions).mkString("\t")
	  }	  
	  
	  val uncanonicalized = tLines(_uncanonicalized).map(lsv => new VCFvariant(lsv))
	  val canonicalized = tLines(_canonicalized).map(lsv => new VCFvariant(lsv))
	  val pw = new PrintWriter(outputDir+"Emu_metrics.txt")
	  pw.println("#This file can be used to visualize Emu's canonicalization performance.\n" +
	      "#Simply copy and paste the data to an excel-like software and graph it using a 2-D graph chart.\n" +
	      "#Complete = LSVs with complete alternative sequence (e.g., no 'N' sequence marker).\n" +
	      "#Incomplete = LSVs with incomplete alternative sequence (e.g., contains 'N' sequence marker).\n#")
	  pw.println("#General Emu Metrics")
	  pw.println("\tAll\tDeletions\tInsertions\tSubstitutions")	  
	  pw.println(getGeneralMetrics(uncanonicalized, "Pre-Emu"))
	  pw.println(getGeneralMetrics(canonicalized, "Post-Emu"))
	  pw.println("\n#Pre-Emu specific metrics")
	  pw.println("\tAll\tDeletions\tInsertions\tSubstitutions")
	  pw.println(getGeneralMetrics(uncanonicalized, "Incomplete"))
	  pw.println(getGeneralMetrics(uncanonicalized, "Complete"))
	  pw.println(getGeneralTotal(uncanonicalized))
	  pw.println("\n#Post-Emu specific metrics")
	  pw.println("\tAll\tDeletions\tInsertions\tSubstitutions")
	  pw.println(getGeneralMetrics(canonicalized, "Incomplete"))
	  pw.println(getGeneralMetrics(canonicalized, "Complete"))
	  pw.println(getGeneralTotal(canonicalized))
	  pw.close
	  println("Emu metrics complete!")
  }
}
  
