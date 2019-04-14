import x10.io.Console;
import x10.array.Array_1;
import x10.array.Array_2;
import x10.io.File;
import x10.io.FileReader;
import x10.io.ReaderIterator;
import x10.util.Pair;
import x10.util.StringBuilder;
import x10.util.Timer;

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
      line = ri.next().substring(1n).trim();
    }

    val i:Int;
    val ind = line.indexOf(' ');
    if (ind.equals(-1n)) {
      i = Int.parse(line);
      line = "";
    } else {
      i = Int.parse(line.substring(0n, ind));
      line = line.substring(ind).trim();
    }
    return i;
  }
}

class CharParser {
  var line:String;

  public def this(s:String) {
    line = s.trim();
  }

  public def hasNext() {
    return line.length() != 0n;
  }

  public def next() {
    val c:Char;
    val ind = line.indexOf(' ');
    if (ind.equals(-1n)) {
      c = line.charAt(0n);
      line = "";
    } else {
      c = line.substring(0n, ind).charAt(0n);
      line = line.substring(ind).trim();
    }
    return c;
  }
}

public class SmithWaterman {
  var n:Long; // Length of a
  var m:Long; // Length of b
  var a:String;
  var b:String;
  var u:Long; // Gap extension penalty
  var v:Long; // Gap opening penalty
  var alphabet:String; // Amino acids
  var w:Array_1[Long]{self!=null};
  var H:Array_2[Cell]{self!=null};
  var S:Array_2[Int]{self!=null};
  var maxH:Cell;

  static struct Cell(score:Long, x:Long, y:Long) {}

  public def this() {
    S = new Array_2[Int](0, 0);
    H = new Array_2[Cell](0, 0);
    w = new Array_1[Long](0);
    maxH = Cell(0, 0, 0);
  }

  def maxTwo(i:Long, j:Long) {
    if (i > j) {
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
        Console.OUT.printf("(%d)", H(i, j).score);
      }
      Console.OUT.println();
    }
  }

  def printS() {
    for (i in 0..(alphabet.length()-1)) {
      for (j in 0..(alphabet.length()-1)) {
        Console.OUT.printf("%d ", S(i, j));
      }
      Console.OUT.println();
    }
  }

  def parseS(fr:FileReader) {
    val ip = new IntParser(fr);
    for (i in 0..(alphabet.length()-1)) {
      for (j in 0..(alphabet.length()-1)) {
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
      0n,
      allLines.length()-1n);

    if (isFirstSeq) {
      n = allLines.length();
      a = allLines;
    } else {
      m = allLines.length();
      b = allLines;
    }
  }

  def skipComments(fr:FileReader) {
    // Read until "text"
    var line:String = fr.readLine();
    while (line.compareTo("text") != 0n) {
      line = fr.readLine();
    };

    // Skip "text"
    line = fr.readLine();

    // Skip comments
    while (true) {
      // Skip leading '@', if any
      if (line.charAt(0n) == '@') {
        line = line.substring(1n);
      }

      // Ignore all comments
      if (line.charAt(0n) == '#') {
        line = fr.readLine();
      } else {
        break;
      }
    }

    return new Pair(fr, line);
  }

  def skipFile(filename:String, isReadAlphabet:Boolean):FileReader {
    val file = new File(filename);
    var fr:FileReader = file.openRead();
    val pair = skipComments(fr);
    fr = pair.first;
    if (isReadAlphabet) {
      val line = pair.second;
      val cp = new CharParser(line);
      val sb = new StringBuilder();
      while (cp.hasNext()) {
        sb.add(cp.next());
      }
      alphabet = sb.toString();
    }
    return fr;
  }

  def initH() {
    H = new Array_2[Cell](n+1, m+1);
  }

  def fillH() {
    for (i in 1..n) {
      for (j in 1..m) {
        var maxK:Long = 0;
        for (k in 1..(i-1)) {
          maxK = maxTwo(maxK, H(k, j).score-w(i-k));
        }
        var maxL:Long = 0;
        for (l in 1..(j-1)) {
          maxL = maxTwo(maxL, H(i, l).score-w(j-l));
        }

        val pair = maxFour(H(i-1, j-1).score + S(
            alphabet.indexOf(a.charAt((i-1) as Int)),
            alphabet.indexOf(b.charAt((j-1) as Int))),
          maxK,
          maxL,
          0);

        var x:Long = 0;
        var y:Long = 0;
        if (pair.second == 0) {
          x = i-1;
          y = j-1;
        } else if (pair.second == 1) {
          x = i-1;
          y = j;
        } else if (pair.second == 2) {
          x = i;
          y = j-1;
        }

        H(i, j) = new Cell(pair.first, x, y);
        if (pair.first > maxH.score) {
          maxH = new Cell(pair.first, i, j);
        }
      }
    }
    Console.OUT.println("----3----");
    Console.OUT.println("----3----");
  }

  def initW() {
    if (n > m) {
      w = new Array_1[Long](n+1);
    } else {
      w = new Array_1[Long](m+1);
    }
  }

  def fillW() {
    for (i in 1..(w.size-1)) {
      w(i) = u*i+v;
    }
  }

  def backtrackH(pair:Pair[StringBuilder, StringBuilder], i:Long, j:Long) {
    val cell = H(i, j);
    if (cell.score == 0) {
      return pair;
    }

    val k = cell.x;
    val l = cell.y;
    var sb1:StringBuilder = new StringBuilder();
    var sb2:StringBuilder = new StringBuilder();

    if (i-k != 1) {
      sb1 = pair.first.add('-');
    } else {
      sb1 = pair.first.add(a.charAt(k as Int));
    }
    if (j-l != 1) {
      sb2 = pair.second.add('-');
    } else {
      sb2 = pair.second.add(b.charAt(l as Int));
    }

    return backtrackH(new Pair[StringBuilder, StringBuilder](sb1, sb2),
      k,
      l);
  }

  def initS() {
    S = new Array_2[Int](alphabet.length(), alphabet.length());
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

    val frA = sw.skipFile(args(0), false);
    val frB = sw.skipFile(args(1), false);
    val frS = sw.skipFile(args(2), true);

    sw.v = Long.parse(args(3));
    sw.u = Long.parse(args(4));

    sw.parseSeq(frA, true);
    Console.OUT.println(sw.n);
    Console.OUT.println(sw.a);
    sw.parseSeq(frB, false);
    Console.OUT.println(sw.m);
    Console.OUT.println(sw.b);


    sw.initS();
    sw.parseS(frS);
    sw.printS();
    sw.initW();
    sw.fillW();

    sw.initH();
    val fillStart = Timer.milliTime();
    sw.fillH();
    val fillStop = Timer.milliTime();
    sw.printH();

    val backtrackStart = Timer.nanoTime();
    val pair = sw.backtrackH(
      Pair[StringBuilder, StringBuilder](
        new StringBuilder(),
        new StringBuilder()),
      sw.maxH.x,
      sw.maxH.y);
    val backtrackStop = Timer.nanoTime();
    Console.OUT.printf("%s\n%s\n", pair.first, pair.second);
    Console.OUT.println("Timing (ns):");
    Console.OUT.printf("fill: %d\n", fillStop - fillStart);
    Console.OUT.printf("backtrack: %d\n", backtrackStop - backtrackStart);
    Console.OUT.printf("maxH: %d\n", sw.maxH.score);
    frA.close();
    frB.close();
    frS.close();
  }
}

