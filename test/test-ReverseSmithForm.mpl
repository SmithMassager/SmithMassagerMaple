read "SmithMassager/reverseSmithForm.mpl";
macro(Equal=LinearAlgebra:-Equal);
macro(Mod=LinearAlgebra:-Modular:-Mod);
macro(Multiply=LinearAlgebra:-Modular:-Multiply);
macro(DiagonalMatrix=LinearAlgebra:-DiagonalMatrix);
kernelopts(assertlevel=1);

#------------------------ Begin ReverseSmithForm --------------------------------
checkReverseSmithForm := proc()
  local A, S, U, V, SNew, UNew, VNew, s;
  s := 15;
  A := Mod(s, Matrix([[1,2,3],[4,5,6],[7,8,9]]), integer);
  S := Mod(s, Vector([1,3,0]), integer);
  U := Mod(s, Matrix([[1,0,0],[4,14,0],[1,13,1]]), integer);
  V := Mod(s, Matrix([[1,13,1],[0,1,13],[0,0,1]]), integer);

  SNew, UNew, VNew := reverseSmithForm(s, S, U, V, 3, 3);
  ASSERT(Equal(SNew, Vector([0,3,1])));
  ASSERT(Equal(UNew, Matrix([[1,13,1],[4,14,0],[1,0,0]])));
  ASSERT(Equal(VNew, Matrix([[1,13,1],[13,1,0],[1,0,0]])));
  ASSERT(Equal(Multiply(s, Multiply(s, UNew, A), VNew), DiagonalMatrix(SNew)));
end proc;
#------------------------ End ReverseSmithForm ----------------------------------

#------------------------ Begin ReconstructM --------------------------------
checkReconstructM := proc()
  local A, S, U, V, s, M;
  s := 15;
  A := Mod(s, Matrix([[1,2,3],[4,5,6],[7,8,9]]), integer);
  S := Mod(s, Vector([1,3,0]), integer);
  U := Mod(s, Matrix([[1,0,0],[4,14,0],[1,13,1]]), integer);
  V := Mod(s, Matrix([[1,13,1],[0,1,13],[0,0,1]]), integer);

  M := reconstructM(s, S, A, V, 3, 3);
  ASSERT(Equal(Multiply(s, A, V), modp(M . DiagonalMatrix(S), s)));
end proc;
#------------------------ End ReconstructM ----------------------------------


checkReverseSmithForm();
checkReconstructM();
