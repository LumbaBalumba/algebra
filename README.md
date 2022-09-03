# algebra
A library that includes complex numbers, vectors and matrices
### *complex*:
 This class contains four constructors (default constructor, copy constructor, double-complex cast and a basic constructor for initialising both real and imagined components), basic arithmetic operations with complex numbers, real-checker function, module and argument functions, and i/o.
 
### *vec*:
 A complex vector class. Contains basic arithmetic operations, real-checker function, fill function to set all vector components to a given complex value, length function (it only works for real vectors), simple constructors (copy and size constructor).
 
### *matrix*:
  This class contains four constructors (default, rows and cols, copy and scalar matrix constructor), rows and cols function, basic arithmetic operations (matrix multiplication function works as O(n^3), later it will be implemented as Strassen method), rows and cols swapping, adding and multiplying functions, matrix trace, transposed matrix function, upper-triangle and lower-triangle forms, conjugate matrix, matrix determinant, angebraic complement to a matrix element, adjusted matrix and inverted matrix (adjusted and inverted matrices are working as O(n^5), Jordan-Gauss method will be implemented later and will be O(n^3)).
