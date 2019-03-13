import x10.io.Console;
import x10.array.Array_1;
import x10.array.Array_2;
import x10.io.File;
import x10.io.FileReader;
import x10.io.ReaderIterator;
import x10.util.Pair;

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
  var n:Long; // Length of a
  var m:Long; // Length of b
  var a:String;
  var b:String;
  var u:Long; // Gap extension penalty
  var v:Long; // Gap opening penalty
  val nS:Long; // Number of amino acids
  val alphabet:String; // Amino acids
  var w:Array_1[Long]{self!=null};
  var H:Array_2[Cell]{self!=null};
  var S:Array_2[Int]{self!=null};
  var maxH:Cell;

  static struct Cell(score:Long, x:Long, y:Long) {}

  public def this() {
    nS = 24;
    alphabet = "ARNDCQEGHILKMFPSTWYVBZX*";
    S = new Array_2[Int](nS, nS);
    H = new Array_2[Cell](0, 0);
    w = new Array_1[Long](0);
    maxH = Cell(0, 0, 0);
  }

  def maxTwo(i:Long, j:Long) {
    if (i.compareTo(j) > 0) {
      return i;
    } else {
      return j;
    }
  }

  def maxFour(i:Long, j:Long, k:Long, l:Long) {
    var ind:Long = 0;
    val m1:Long = maxTwo(i, j);
    val m2:Long = maxTwo(k, l);
    var m3:Long = maxTwo(m1, m2);

    if (i == m3) {
      ind = 0;
    } else if (j == m3) {
      ind = 1;
    } else if (k == m3) {
      ind = 2;
    } else if (l == m3) {
      ind = 3;
    }

    return new Pair(m3, ind);
  }

  def printH() {
    for (i in 0..n) {
      for (j in 0..m) {
        Console.OUT.printf("(%d %d %d)", H(i, j).score, H(i, j).x, H(i, j).y);
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

  def parseSeq(fr:FileReader, isFirstSeq:Boolean) {
    var allLines:String = "";
    for (line in fr.lines()) {
      allLines += line;
    }
    allLines = allLines.substring(
      Int.operator_as(0),
      allLines.length()-Int.operator_as(1));

    if (isFirstSeq) {
      n = allLines.length();
      a = allLines;
    } else {
      m = allLines.length();
      b = allLines;
    }
  }

  def skipFile(filename:String, lineAfter:Long):FileReader {
    val file = new File(filename);
    val fr = file.openRead();
    for (i in 1..lineAfter) {
      fr.readLine();
    }
    return fr;
  }

  def initH() {
    H = new Array_2[Cell](n+1, m+1);
  }

  def fillH() {
    var maxK:Long = 0;
    var maxKInd:Long = 1;
    var maxLInd:Long = 1;
    var maxL:Long = 0;

    for (i in 1..n) {
      for (j in 1..m) {
        for (k in 1..i) {
          var oldK:Long = maxK;
          maxK = maxTwo(maxK, H(i-k, j).score-w(k-1));
          if (maxK == oldK) {
            maxKInd = k;
          }
        }
        for (l in 1..j) {
          var oldL:Long = maxL;
          maxL = maxTwo(maxL, H(i, j-1).score-w(l-1));
          if (maxL == oldL) {
            maxLInd = l;
          }
        }

        val pair = maxFour(H(i-1, j-1).score + S(
            alphabet.indexOf(a.charAt(Int.operator_as(i-1))),
            alphabet.indexOf(b.charAt(Int.operator_as(j-1)))),
          maxK,
          maxL,
          0);

        var x:Long = 0;
        var y:Long = 0;
        if (pair.second == 0) {
          x = i-1;
          y = j-1;
        } else if (pair.second == 1) {
          x = i-maxKInd;
          y = j;
        } else if (pair.second == 2) {
          x = i;
          y = i-maxLInd;
        }
        // Console.OUT.println("" + pair.first + " " + x + " " + y);

        H(i, j) = new Cell(pair.first, x, y);
        if (pair.first > maxH.score) {
          maxH = new Cell(pair.first, i, j);
        }
      }
    }
  }

  def initW() {
    if (n > m) {
      w = new Array_1[Long](n);
    } else {
      w = new Array_1[Long](m);
    }
  }

  def fillW() {
    for (i in 1..w.rank()) {
      w(i) = u*i+v;
    }
  }

  public static def main(args:Rail[String]):void {
    if (args.size != 5) {
      Console.OUT.println("Usage: SmithWaterman
        fileSeqA
        fileSeqB
        fileSubst
        openPenalty
        extendPenalty");
      return;
    }

    val sw = new SmithWaterman();

    val frA = sw.skipFile(args(0), 23);
    val frB = sw.skipFile(args(1), 23);
    val frS = sw.skipFile(args(2), 36);

    sw.v = Long.parse(args(3));
    sw.u = Long.parse(args(4));

    sw.parseS(frS);
    sw.printS();
    sw.parseSeq(frA, true);
    Console.OUT.println(sw.n);
    Console.OUT.println(sw.a);
    sw.parseSeq(frB, false);
    Console.OUT.println(sw.m);
    Console.OUT.println(sw.b);

    sw.initW();
    sw.fillW();

    sw.initH();
    sw.fillH();
    sw.printH();

    Console.OUT.printf("%d %d %d\n", sw.maxH.score, sw.maxH.x, sw.maxH.y);

    frA.close();
    frB.close();
    frS.close();
  }
}
