macro(MatrixInverse=LinearAlgebra:-MatrixInverse);
macro(ModInverse=LinearAlgebra:-Modular:-Inverse);
macro(Identity=LinearAlgebra:-Modular:-Identity);
macro(Multiply=LinearAlgebra:-Multiply);
macro(IntegerLinearSolve=LinearAlgebra:-Modular:-IntegerLinearSolve);

# IntCertificate(s, A, B, n, m)
# 
# Input:
#   A, n x n matrix over Z, nonsingular.
#   s, positive integer.
#   B, n x m matrix over Z(s).
#   n, m, positive integer.
#
#
# Output:
#   If sA^-1B is integral then return modp(sA^-1B, s) otherwise return false.
#
naiveIntCertificate := proc(s, A, B, n, m):
  local res;
  if (not type(A, 'Matrix'(integer))) then
    return false;
  end if;

  res := LinearAlgebra:-LinearSolve(A, s*B);
  if (type(res, 'Matrix'(integer))) then
    return modp(res, s);
  else
    return false;
  end if;
end proc;

# UniCert(A, n)
#
# Input:
#   A, n x n integer matrix in Z, nonsingular.
#
# Output:
#   True if |det(A)| = 1 else false
#
naiveUniCert := proc(A, n)
  if (type(naiveIntCertificate(1, A, Identity(2, n, integer), n, n), Matrix)) then
    return true;
  else 
    return false;
  end if;
end proc; 
