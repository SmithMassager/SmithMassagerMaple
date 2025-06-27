macro(HermiteForm=LinearAlgebra:-HermiteForm);
macro(MatrixInverse=LinearAlgebra:-MatrixInverse);
macro(RandomMatrix=LinearAlgebra:-RandomMatrix);
with(LinearAlgebra);
with(ArrayTools);

macro(MX=1000);
kernelopts(assertlevel=1);
unfoundAns := Array([]);
unfoundCount := 0;
tot := 0;

# randomUnimodularMatrix(n)
#
# Input:
#   n, positive integer.
#
# Output:
#   n x n unimodular matrix.
# 
randomUnimodularMatrix := proc(n)
  local A;
  while (true) do
    A := RandomMatrix(n, n);
    if (Determinant(A) <> 0) then break; end if;
  end do;
  return A . MatrixInverse(HermiteForm(A));
end proc;

# checkUnitUpperTriangular(T)
# Input:
#   T, n x n matrix.
#
# Output: True iff T is unit upper triangular.
#
checkUnitUpperTriangular := proc(T, n)
  local i, j;
  for i from 1 to n do
    for j from 1 to i do
      if (j = i and T[i,j] <> 1) then return false;
      elif (j <> i and T[i,j] <> 0) then return false end if;
    end do;
  end do;
  return true;
end proc;

# checkSmithMassager(A, n, U, M, T, S)
#
# Input:
#   A, n x n integer matrix nonsingular.
#   U, M, T, S, n x n integer matrix.
#
# Output: True iff (U, M, T, S) is a (0, n)-index massager for A.
#
checkValidMassager := proc(A, n, U, M, T, S)
  local B, botRight, AMSInv, Ans;
  ASSERT(checkUnitUpperTriangular(T, n));
  B := Identity(2, 2 * n, integer);
  B[1..n, 1..n] := A;
  B[n+1..2*n, 1..n] := U;
  AMSInv := A . M . MatrixInverse(S);
  B[1..n, n+1..2*n] := AMSInv;
  botRight := (T + U . M) . MatrixInverse(S);
  B[n+1..2*n, n+1..2*n] := botRight;

  ASSERT(type(B, 'Matrix(integer)'));
  ASSERT(abs(Determinant(B)) = 1);
  Ans := SmithForm(A);
  ASSERT(Equal(S, Ans));
  # Check AM == 0 cmod S and -UM = T cmod S. Since (T + UM) S^-1 is integer
  # and T is upper triangular.
  ASSERT(Equal(cmod(A . M, S), Matrix(n, n)));
  ASSERT(Equal(cmod(-U . M, S), cmod(T, S)));
end proc;

recordMatrices := proc(A, n)
  global unfoundCount, unfoundAns;
  Append(unfoundAns, A, inplace = true);
  ++unfoundCount;
end proc;

randomNonsingular := proc(n)
  local roll, A;
  roll := rand(1..MX);
  A := Matrix(n, n, (i,j) -> if i = j then roll() else 0 end if);

  A := randomUnimodularMatrix(n) . A . randomUnimodularMatrix(n);

  return A;
end proc;
