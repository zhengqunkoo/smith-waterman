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

/**
 * Parse ints of a file reader.
 *
 * @param fr: file reader.
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

  // Return next parsed int across newlines.
  public def next() {
    // Read next line if current line empty.
    if (line.equals("")) {
      line = ri.next().substring(1n).trim();
    }

    // Try to find space in line.
    val i:Int;
    val ind = line.indexOf(' ');
    if (ind.equals(-1n)) {
      // No space found means line has last integer.
      // Set line to empty string explicitly.
      i = Int.parse(line);
      line = "";
    } else {
      // If space found, get substring from start of line till space.
      // Get rest of line, then trim off extra spaces.
      i = Int.parse(line.substring(0n, ind));
      line = line.substring(ind).trim();
    }
    return i;
  }
}

/**
 * Parse chars of a string.
 *
 * @param s: string to be parsed.
 */
class CharParser {
  var line:String;

  public def this(s:String) {
    line = s.trim();
  }

  public def hasNext() {
    return line.length() != 0n;
  }

  // Return next parsed char across newlines.
  public def next() {
    val c:Char;
    // Try to find space in line.
    val ind = line.indexOf(' ');
    if (ind.equals(-1n)) {
      // No space found means line has last integer.
      // Set line to empty string explicitly.
      c = line.charAt(0n);
      line = "";
    } else {
      // If space found, get substring from start of line till space.
      // Get rest of line, then trim off extra spaces.
      c = line.substring(0n, ind).charAt(0n);
      line = line.substring(ind).trim();
    }
    return c;
  }
}


/**
 * Smith-Waterman algorithm.
 */
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

  // A cell is an element of a matrix.
  static struct Cell(score:Long, x:Long, y:Long) {}

  public def this() {
    S = new Array_2[Int](0, 0);
    H = new Array_2[Cell](0, 0);
    w = new Array_1[Long](0);
    maxH = Cell(0, 0, 0);
  }

  // Return highest long of two longs, i and j.
  def maxTwo(i:Long, j:Long) {
    if (i > j) {
      return i;
    } else {
      return j;
    }
  }

  // If [i, j, k, l] were in an array,
  // Return highest long, and index of that highest long in that array.
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

  // Print cells of the H matrix.
  def printH() {
    for (i in 0..n) {
      for (j in 0..m) {
        Console.OUT.printf("(%d)", H(i, j).score);
      }
      Console.OUT.println();
    }
  }

  // Print values of the S matrix.
  def printS() {
    for (i in 0..(alphabet.length()-1)) {
      for (j in 0..(alphabet.length()-1)) {
        Console.OUT.printf("%d ", S(i, j));
      }
      Console.OUT.println();
    }
  }

  // Read integers from file reader @param fr into the S matrix.
  // S is of dimension alphabet.length by alphabet.length.
  def parseS(fr:FileReader) {
    val ip = new IntParser(fr);
    for (i in 0..(alphabet.length()-1)) {
      for (j in 0..(alphabet.length()-1)) {
        S(i, j) = ip.next();
      }
    }
  }

  // If @param isFirstSeq, then set @var a to all lines, and @var n to length
  // of a.
  // If @param not isFirstSeq, then set @var b to all lines, and @var m to
  // length of b.
  def parseSeq(fr:FileReader, isFirstSeq:Boolean) {
    // Read all lines in file reader.
    var allLines:String = "";
    for (line in fr.lines()) {
      allLines += line;
    }

    // Ignore newline.
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

  // @return pair:
  //   file reader pointing at second non-comment line
  //   first non-comment line
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
      // Ignore leading '@', if any
      if (line.charAt(0n) == '@') {
        line = line.substring(1n);
      }

      // Ignore all comments.
      // Comments start with '#'.
      if (line.charAt(0n) == '#') {
        line = fr.readLine();
      } else {
        break;
      }
    }

    return new Pair(fr, line);
  }

  // Open file @param filename.
  // If @param isReadAlphabet,
  // get the first non-comment line
  // parse this line into characters, assign to @var alphabet.
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

  // Initialize cells of the H matrix.
  def initH() {
    H = new Array_2[Cell](n+1, m+1);
  }

  // Fill in each cell of H.
  def fillH() {
    for (i in 1..n) {
      for (j in 1..m) {
        // Get maximum score of cells in same column with penalties, @var maxK.
        // Get maximum score of cells in same row with penalties, @var maxL.
        var maxK:Long = 0;
        for (k in 1..(i-1)) {
          maxK = maxTwo(maxK, H(k, j).score-w(i-k));
        }
        var maxL:Long = 0;
        for (l in 1..(j-1)) {
          maxL = maxTwo(maxL, H(i, l).score-w(j-l));
        }

        // Get the cell with maximum score: either
        // 1) Diagonal neighbor, score is diag neighbor's score plus S matrix.
        // 2) Column neighbor, score maxK.
        // 3) Row neighbor, score maxL.
        // 4) current, score 0.
        val pair = maxFour(H(i-1, j-1).score + S(
            alphabet.indexOf(a.charAt((i-1) as Int)),
            alphabet.indexOf(b.charAt((j-1) as Int))),
          maxK,
          maxL,
          0);

        // Store coordinates of the largest scoring neighbor in (x, y).
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

        // Store cell with maximum score, and coordinates of largest scoring
        // neighbor, in H matrix.
        // Update maxH cell if score of this cell exceeds that of maxH.
        H(i, j) = new Cell(pair.first, x, y);
        if (pair.first > maxH.score) {
          maxH = new Cell(pair.first, i, j);
        }
      }
    }
    Console.OUT.println("----3----");
    Console.OUT.println("----3----");
  }

  // Initialize @var w with the larger of two matrix dimensions.
  def initW() {
    if (n > m) {
      w = new Array_1[Long](n+1);
    } else {
      w = new Array_1[Long](m+1);
    }
  }

  // Fill @var w with a linear value.
  def fillW() {
    for (i in 1..(w.size-1)) {
      w(i) = u*i+v;
    }
  }

  // Backtrack through the H matrix, starting from cell (@param i, j).
  // And storing the character corresponding to each cell in @var sb1, sb2.
  // @param pair: first aligned sequence, second aligned sequence.
  def backtrackH(pair:Pair[StringBuilder, StringBuilder], i:Long, j:Long) {
    val cell = H(i, j);

    // Stop if score of current cell is zero.
    if (cell.score == 0) {
      return pair;
    }

    // Get coordinates of largest scoirng neighbor in new temp variables.
    val k = cell.x;
    val l = cell.y;

    // Create new vars for old string builders with appended character.
    var sb1:StringBuilder = new StringBuilder();
    var sb2:StringBuilder = new StringBuilder();

    // Compare current cell coord to largest scoring neighbor coord.
    // Add gap '-' to string builder if the respective coord does not change.
    // Else, append appropriate character from a or b.
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

    // Recursively call with updated string builders and updated cell coords.
    return backtrackH(new Pair[StringBuilder, StringBuilder](sb1, sb2),
      k,
      l);
  }

  // S is of dimension alphabet.length by alphabet.length.
  def initS() {
    S = new Array_2[Int](alphabet.length(), alphabet.length());
  }

  // Program entry point with @param args like in the usage prompt.
  public static def main(args:Rail[String]):void {

    // Usage prompt.
    if (args.size != 5) {
      Console.OUT.println("Usage: SmithWaterman
        fileSeqA
        fileSeqB
        fileSubst
        openPenalty
        extendPenalty");
      return;
    }

    // Create new SW object.
    val sw = new SmithWaterman();

    // Get file readers pointed at the second non-comment line.
    val frA = sw.skipFile(args(0), false);
    val frB = sw.skipFile(args(1), false);
    val frS = sw.skipFile(args(2), true);

    // Get the gap opening and gap extension penalties.
    sw.v = Long.parse(args(3));
    sw.u = Long.parse(args(4));

    // Parse and assign the sequences to the correct variables.
    sw.parseSeq(frA, true);
    Console.OUT.println(sw.n);
    Console.OUT.println(sw.a);
    sw.parseSeq(frB, false);
    Console.OUT.println(sw.m);
    Console.OUT.println(sw.b);

    // Initialize, fill, and print matrices.
    sw.initS();
    sw.parseS(frS);
    sw.printS();
    sw.initW();
    sw.fillW();

    // Time filling H and backtracking across H.
    sw.initH();
    val fillStart = Timer.nanoTime();
    sw.fillH();
    val fillStop = Timer.nanoTime();
    sw.printH();

    val backtrackStart = Timer.nanoTime();
    val pair = sw.backtrackH(
      Pair[StringBuilder, StringBuilder](
        new StringBuilder(),
        new StringBuilder()),
      sw.maxH.x,
      sw.maxH.y);
    val backtrackStop = Timer.nanoTime();

    // Print alignment results and time values.
    Console.OUT.printf("%s\n%s\n", pair.first, pair.second);
    Console.OUT.println("Timing (ns):");
    Console.OUT.printf("fill: %d\n", fillStop - fillStart);
    Console.OUT.printf("backtrack: %d\n", backtrackStop - backtrackStart);
    Console.OUT.printf("maxH: %d\n", sw.maxH.score);

    // Close all file readers.
    frA.close();
    frB.close();
    frS.close();
  }
}
