# Smith-Waterman Algorithm
In the x10 language, for the National University of Singapore's CS3211 module.

## Branch Notes
This is the master branch.

## Usage
- To avoid badly formatted stdout, redirect to file.
  ```
  x10c++ SmithWaterman.x10 \
  && ./a.out sequences/p53_human.fasta,v sequences/p53_mouse.fasta,v \
  matrices/BLOSUM62,v 5 1 > out \
  && cat out
  ```
  This command short circuits if there are any compilation errors.
- From [Wikipedia's gap example].
  ```
  x10c++ SmithWaterman.x10 \
  && ./a.out sequences/wiki_gap1.fasta,v sequences/wiki_gap2.fasta,v \
  matrices/EDNAFULL,v 5 1 > out \
  && cat out
  ```
- From [Wikipedia's subst example].
  ```
  x10c++ SmithWaterman.x10 \
  && ./a.out sequences/wiki_subst1.fasta,v sequences/wiki_subst2.fasta,v \
  matrices/wiki_DNA3,v 0 2 > out \
  && cat out
  ```

## Roadmap
1. Implement sequential algorithm.

[Wikipedia's gap example]: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Gap_penalty_example
