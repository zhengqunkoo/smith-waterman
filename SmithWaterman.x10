import x10.io.Console;
import x10.array.Array_2;
import x10.io.File;

public class SmithWaterman {
  val N:Long; // Length of A
  val M:Long; // Length of B
  var H:Array_2[Double]{self!=null};

  public def this() {
    N = 10;
    M = 10;
    H = new Array_2[Double](N+1, M+1);

    // TODO not really needed
    for (i in 1..N) {
      H(i, 0) = 0;
    }
    for (j in 0..M) {
      H(0, j) = 0;
    }
  }

  def print() {
    for (i in 0..N) {
      for (j in 0..M) {
        Console.OUT.printf("%1.4f ", H(i, j));
      }
      Console.OUT.println();
    }
  }

  public static def main(args:Rail[String]):void {
    if (args.size != 1) {
      Console.OUT.println("Usage: SmithWaterman filename");
      return;
    }

    val file = new File(args(0));
    val fr = file.openRead();
    for (i in 1..35) {
      fr.readLine();
    }
    for (line in fr.lines()) {
      Console.OUT.println(fr.readLine());
    }

    val sw = new SmithWaterman();
    sw.print();
  }
}
