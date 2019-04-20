import x10.io.Console;
import x10.array.Array_1;
import x10.array.Array_2;
import x10.io.File;
import x10.io.FileReader;
import x10.io.ReaderIterator;
import x10.util.Pair;
import x10.util.StringBuilder;
import x10.util.Timer;
import x10.lang.Rail;

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

// Base vector class
// Contains methods and fields for vectorized operations
class Vector {
  val VECTOR_SIZE = 8;
  var v: Rail[Long]{self!=null};                // Contains vector values
  public operator this(i:Long) = this.v(i);     // Vector element accessor
  public operator this(i:Long) = (newval: Long) // Vector element assignment accessor
  {
    v(i) = newval;
  }

  // Literal constructor
  public def this(i: Long){
    v = new Rail[Long](VECTOR_SIZE);
    for (var p: Long = 0; p < VECTOR_SIZE; p++)
    {
      v(p) = i;
    }
  }

  // Copy constructor
  public def this(vec: Rail[Long]){
    this.v = new Rail[Long](VECTOR_SIZE);
    for (var p: Long = 0; p < VECTOR_SIZE; p++)
    {
      this.v(p) = vec(p);
    }
  }

  // Subtraction operator
  // Returns newly allocated vector
  public operator this - (vect:Vector) = new Vector(new Rail[Long](
    VECTOR_SIZE, (i:Long) => this.v(i) - vect.v(i)
  ));

  // Addition operator
  // Returns newly allocated vector
  public operator this + (vect:Vector) = new Vector(new Rail[Long](
    VECTOR_SIZE, (i:Long) => this.v(i) + vect.v(i)
  ));

  //Element-wise copy all elements of @vec to this
  // Used instead of assignment operator
  public def copy(vec: Vector) {
    for (var p: Long = 0; p < VECTOR_SIZE; p++)
    {
      this.v(p) = vec(p);
    }
  }

  // Compares each element of @vec1 and @vec2
  // Stores highest element in this vector
  // Returns pair containing index of highest value and highest value respectively
  public def maxVector(vec1: Vector, vec2: Vector){
    var highest: Long = this.v(0);
    var index: Long = 0;
    for (var p: Long = 0; p < VECTOR_SIZE; p++)
    {
      if (vec1(p) < vec2(p))
        this.v(p) = vec2(p);
      else
        this.v(p) = vec1(p);

      if (this.v(p) > highest)
      {
        highest = this.v(p);
        index = p;
      }

    }

    return new Pair(index, highest);
  }

  // Subtract each element of vector with literal
  public def subtractLit(i: Long){
    for (var j: Long = 0n; j < VECTOR_SIZE; j++){
      this.v(j) = this.v(j) - i;
    }
  }

  // Add literal to each element
  public def addLit(i: Long){
    for (var j: Long = 0n; j < VECTOR_SIZE; j++){
      this.v(j) = this.v(j) + i;
    }
  }

  // Shift elements in vector by @param i steps to the right
  public def rightShift(i: Long) {
    for (var j: Long = 0n; j < i; j++)
    {
      rightShiftOne();
    }
  }

  //Shift elements in vector to the right by one
  def rightShiftOne(){
    for (var p: Long = 0; p < VECTOR_SIZE-2; p++)
    {
      v(p) = v(p+1);
    }
    v(VECTOR_SIZE-1) = 0;
  }

  // Shift elements in vector by @param i steps to the left
  public def leftShift(i: Long) {
    for (var j: Long = 0n; j < i; j++)
    {
      leftShiftOne();
    }

  }

  //Shift elements in vector to the left by one
  def leftShiftOne(){
    for (var p: Long = VECTOR_SIZE-1; p > 0; p--)
    {
      v(p) = v(p-1);
    }
    v(0) = 0;
  }
}


/**
 * Smith-Waterman algorithm.
 */
public class SmithWatermanParVect {

  val VECTOR_SIZE = 8;
  var n:Long; // Length of a
  var m:Long; // Length of b
  var a:String;
  var b:String;
  var u:Long; // Gap extension penalty
  var v:Long; // Gap opening penalty
  var alphabet:String; // Amino acids
  var H:Array_2[Cell]{self!=null};
  var HH:Rail[Vector];
  var EE:Rail[Vector];
  var S:Array_2[Int]{self!=null};
  var maxH:Cell;

  // A cell is an element of a matrix.
  static struct Cell(score:Long, x:Long, y:Long) {}



  // Initialize cells of the H matrix.
  def initH() {
    val vectorCount = n / VECTOR_SIZE + 1;
    EE = new Rail[Vector](vectorCount+1);
    HH = new Rail[Vector](vectorCount+1);
    H = new Array_2[Cell](n+1, m+1);

  }

  // Fill in each cell of H using Rognes' Smith-Waterman algorithm.
  def fillH() {
    //Initialize temporary varibles for algorithm
    var currentChar:Char;                //One database sequence symbol (c)
    val vectorCount = (n / VECTOR_SIZE);  //Vectors (H, X, E, F, T1, T2) and Array (HH, EE)
    var tempH: Vector = new Vector(0);
    var tempE: Vector = new Vector(0);
    var temp1: Vector = new Vector(0);
    var temp2: Vector = new Vector(0);
    var F: Vector = new Vector(0);
    var tempF: Vector = new Vector(0);
    var X: Vector = new Vector(0);
    var zeroVect: Vector = new Vector(0);
    var score: Vector = new Vector(0);

    //Default vectors containing gap penalty variables
    var gapOpen: Vector = new Vector(v);  //Gap open penalty vector (Vector.subtractLit(v) is used instead)
    var gapExt: Vector = new Vector(u);   //Gap extension penalty (Vector.subtractLit(u) is used instead)

    for (i in 0..vectorCount){
      HH(i) = new Vector(0); //Initialize HH-array of H-values from previous column
      EE(i) = new Vector(0); // Initialize EE-array of E-values from previous column
    }

    for (var j: Long = 0; j < m; j++) //For each symbol in the DB sequence
    {
      X.copy(zeroVect); //Initialize X-vector for 1. round
      F.copy(zeroVect); //Initialize F-vector for 1. round

      currentChar = b(j as Int); //Get one database symbol

      for (var i: Long = 0; i < vectorCount; i++) // For each vector of 8 matrix cells along query sequence
      {
        tempH.copy(HH(i)); // Load previous H-vector from HH-array
        tempE.copy(EE(i)); // Load previous E-vector from EE-array

        temp1(0) = tempH(7);// Save previous H[7] for use below
        tempH.leftShift(1); // Shift H-vector and OR with H[7] from previous round
        tempH(0) = X(0);    // Shift H-vector and OR with H[7] from previous round
        X(0) = temp1(0);    //save old H[7] in X for next round

        // Add score profile vector to H (and subtract base)
        for(var vInd: Long = 0; vInd < VECTOR_SIZE; vInd++)
        {
          tempH(vInd) = (tempH(vInd) +
          S(alphabet.indexOf(a.charAt((i*VECTOR_SIZE+vInd) as Int)),
            alphabet.indexOf(currentChar)
          ));
        }

        tempH.maxVector(tempH, tempE); //Check if score with database gap is better

        tempF.copy(F);      // Calculate initial F-vector by shifting H and previous F
        F.copy(tempH);      // Calculate initial F-vector by shifting H and previous F
        F.leftShift(1);     // Calculate initial F-vector by shifting H and previous F
        F(0) = tempF(7);    // Calculate initial F-vector by shifting H and previous F
        F.subtractLit(u+v); //Subtract single gap penalty

        //Check if vertical gaps are possible
        var fAboveZero: Boolean = false;
        for(var vInd: Long = 0; vInd < VECTOR_SIZE; vInd++){
          if (F(vInd) > 0)
          {
            fAboveZero = true;
          }
        }


        if (fAboveZero) // Compute correct F-vector if necessary
        {
          temp2.copy(F); //T2 is initial F-vector

          while (fAboveZero) //Repeat while any element of T2 is nonzero
          {
            fAboveZero = false;

            temp2.leftShift(1);    // Shift and subtract gap extension penalty
            temp2.subtractLit(u);  // Shift and subtract gap extension penalty
            F.maxVector(F, temp2); // Update F if new score is higher

            for(var vInd: Long = 0; vInd < VECTOR_SIZE; vInd++){
              if (temp2(vInd) > 0)
              {
                fAboveZero = true; // If any value is still above 0, set to true
              }
            }
          }

          tempH.maxVector(tempH, F); //Update H if vertical gap is better
          F.addLit(v);               // Update F for use in next round
          F.maxVector(F, tempH);     // Update F for use in next round
        }
        else{
          F.copy(tempH); // Update F for use in next round
        }

        tempH.maxVector(tempH, zeroVect); //Maybe this statement is actually not supposed to be here :O


        HH(i).copy(tempH);                       // Store H-vector in HH-array
        EE(i).maxVector(tempH - gapOpen, tempE); // Store E-vector in EE-array
        EE(i).subtractLit(u);                    // Store E-vector in EE-array

        // Save elements in H-matrix
        for (vInd in 0..7) {
          val row: Long = i*VECTOR_SIZE+vInd;
          H(row+1, j+1) = new Cell(tempH(vInd), j+1, row+1);
        }

        //Find highest value found so far
        val max = score.maxVector(tempH, score);

        // Update maxH cell if score of this cell exceeds that of maxH.
        if (max.second > maxH.score){
          maxH = new Cell(max.second, i*VECTOR_SIZE+max.first+1, j+1);
        }
      }
    }

    //Print highest found values
    Console.OUT.println("----3----");
    for (vInd in 0..7) {
      Console.OUT.printf("(%d)", score(vInd));
    }
    Console.OUT.println("\n----3----");
  }

  //Public constructor, zero initialize variables
  public def this() {
    S = new Array_2[Int](0, 0);
    H = new Array_2[Cell](0, 0);
    maxH = Cell(0, 0, 0);
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
  // Also pads first sequnce until it is divisible by vector size
  def parseSWithPadding(fr:FileReader) {
    Console.OUT.println("LENGTH OF ALPHABET = " + alphabet.length());

    S = new Array_2[Int](alphabet.length() + 1, alphabet.length() + 1);
    val ip = new IntParser(fr);
    for (i in 0..(alphabet.length()-1)) {
      for (j in 0..(alphabet.length()-1)) {
        S(i, j) = ip.next();
      }
    }
    var padding:Long = (a.length() % VECTOR_SIZE); //Pad data if a-sequence
                                                   // Does not scale with vector size
    if (padding != 0)
    {
      alphabet = alphabet + "0";

      for (i in 0..padding){
        a = a + '0';
      }

      // Add large negative value to padding character in substitution matrix
      // To avoid padding characters affecting alignment result
      for (i in 0..(alphabet.length() -1))
      {
        S(alphabet.length()-1, i) = -2000n;
        S(i, alphabet.length()-1) = -2000n;
      }

      n = a.length();
    }


    Console.OUT.println("Length of n: " + n);
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

  //Finds highest value of three integers and returns index of highest
  def maxThree(a: Long, b: Long, c: Long) : Long
  {
    val m1: Long = maxTwo(a, b);
    val m2: Long = maxTwo(m1, c);

    if (m2 == a)
      return 0;
    else if (m2 == b)
      return 1;
    else
      return 2;

  }

  // Return highest long of two longs, i and j.
  def maxTwo(i:Long, j:Long) {
    if (i > j) {
      return i;
    } else {
      return j;
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

    val c1 = H(i-1, j-1);
    val c2 = H(i-1, j);
    val c3 = H(i, j-1);

    // Get coordinates of largest scoirng neighbor in new temp variables.
    var k: Long = i-1;
    var l: Long = j-1;
    val maxInd = maxThree(c1.score, c2.score, c3.score);
    if (maxInd == 1)
    {
      l = j;
    }
    else if (maxInd == 2)
    {
      k = i;
    }

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

  // Reverse a string builder.
  // Defining this as a standalone function avoids affecting performance time,
  // Since we can call this function outside of timed code.
  def reverse(sb:StringBuilder):StringBuilder {
    val sbNew:StringBuilder = new StringBuilder();

    // First convert string builder into a rail.
    // Then swap the elements of the rail.
    // Define a temporary buffer to perform the swap.
    val r:Rail[Char] = sb.result().chars();
    val n:Long = sb.length()-1;
    var tmp:Char;
    for (val i in (0..(n/2))) {
      tmp = r(i);
      r(i) = r(n-i);
      r(n-i) = tmp;
    }
    for (val c:Char in r) {
      sbNew.add(c);
    }
    return sbNew;
  }

  // Program entry point with @param args like in the usage prompt.
  public static def main(args:Rail[String]):void {

    // Usage prompt.
    if (args.size != 5) {
      Console.OUT.println("Usage: SmithWatermanTrans
        fileSeqA
        fileSeqB
        fileSubst
        openPenalty
        extendPenalty");
      return;
    }

    // Create new SW object.
    val sw = new SmithWatermanParVect();

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
    sw.parseSWithPadding(frS);
    sw.printS();

    // Time filling H and backtracking across H.
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

    // Print alignment results and time values.
    Console.OUT.printf("%s\n%s\n", sw.reverse(pair.first),
      sw.reverse(pair.second));
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
