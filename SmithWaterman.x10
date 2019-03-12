import x10.io.Console;
import x10.array.Array_2;
import x10.io.File;
import x10.io.FileReader;
import x10.io.ReaderIterator;

/**
 * Amino acids do not include gap codon.
 * See https://www.mathworks.com/help/bioinfo/ref/aminolookup.html
 */

class IntParser {
  val ri:ReaderIterator[String];
  var line:String = "";

  public def this(fr:FileReader) {
    ri = fr.lines();
  }

  public def hasNext() {
    return true;
  }

  public def next() {
    if (line.equals("")) {
      line = ri.next().substring(Int.operator_as(1)).trim();
    }

    val i:Int;
    val ind = line.indexOf(' ');
    // Console.OUT.println("Bef: " + ind + " " + line);
    if (ind.equals(Int.operator_as(-1))) {
      i = Int.parse(line);
      line = "";
    } else {
      i = Int.parse(line.substring(Int.operator_as(0), ind));
      line = line.substring(ind).trim();
    }
    // Console.OUT.println("Aft: " + ind + " " + i + " " + line);
    return i;
  }
}

public class SmithWaterman {
  val N:Long; // Length of A
  val M:Long; // Length of B
  val nS:Long; // Number of amino acids
  var H:Array_2[Double]{self!=null};
  var S:Array_2[Int]{self!=null};

  public def this() {
    N = 10;
    M = 10;
    nS = 24;
    H = new Array_2[Double](N+1, M+1);
    S = new Array_2[Int](nS, nS);
  }

  def printH() {
    for (i in 0..N) {
      for (j in 0..M) {
        Console.OUT.printf("%1.4f ", H(i, j));
      }
      Console.OUT.println();
    }
  }

  def printS() {
    for (i in 0..(nS-1)) {
      for (j in 0..(nS-1)) {
        Console.OUT.printf("%d ", S(i, j));
      }
      Console.OUT.println();
    }
  }

  def parseS(fr:FileReader) {
    val ip = new IntParser(fr);
    for (i in 0..(nS-1)) {
      for (j in 0..(nS-1)) {
        S(i, j) = ip.next();
      }
    }
  }

  public static def main(args:Rail[String]):void {
    if (args.size != 1) {
      Console.OUT.println("Usage: SmithWaterman filename");
      return;
    }

    val file = new File(args(0));
    val fr = file.openRead();
    for (i in 1..36) {
      fr.readLine();
    }

    val sw = new SmithWaterman();
    sw.printH();
    sw.parseS(fr);
    sw.printS();

    fr.close();
  }
}
