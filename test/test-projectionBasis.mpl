read "SmithMassager/projectionBasis.mpl";
macro(RandomMatrix=LinearAlgebra:-RandomMatrix);
kernelopts(assertlevel=1);


#------------------------ Begin ProjectionBasis ---------------------------------
checkProjectionBasis := proc()
  local B, s, J, P, S, U, M;
  B := Matrix([[1, 8, 0, 0],[2, 8, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
  s := 8;
  J := Matrix([[0, 0], [3, 4], [0, 0], [0, 0]]);
  P := modp(s * MatrixInverse(B) . J, s);

  S, U, M := computeProjBasis(P, 2, 2, 8); 

  ASSERT(Equal(S, Vector([1,8])));
  ASSERT(gcd(Determinant(U), s) = 1);
end proc;
#------------------------ End ProjectionBasis ---------------------------------

checkProjectionBasis();
