package abeellab.emu

trait EmuUtils {
  
  def extractID(str:String):String={
    
    str.split("/").dropRight(1).last
    
  }

}