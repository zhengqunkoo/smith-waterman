# Smith-Waterman Algorithm
In the x10 language, for the National University of Singapore's CS3211 module.

## Usage
To avoid badly formatted stdout, redirect to file.
```
x10c++ SmithWaterman.x10 && ./a.out matrices/BLOSUM62,v > out && cat out
```
This command short circuits if there are any compilation errors.  Show these
compilation errors by running `cat out`.

## Roadmap
1. Implement sequential algorithm.
