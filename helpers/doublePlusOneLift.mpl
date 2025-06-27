macro(ModInverse=LinearAlgebra:-Modular:-Inverse);
macro(Multiply=LinearAlgebra:-Multiply);
macro(ModMultiply=LinearAlgebra:-Modular:-Multiply);
macro(MatrixPower=LinearAlgebra:-MatrixPower);
macro(Identity=LinearAlgebra:-Modular:-Identity);
macro(Dimension=LinearAlgebra:-Dimension);
macro(AddMultiple=LinearAlgebra:-Modular:-AddMultiple);
macro(ZeroMatrix=LinearAlgebra:-ZeroMatrix);
macro(Equal=LinearAlgebra:-Equal);
macro(Mod=LinearAlgebra:-Modular:-Mod);
macro(ModDeterminant=LinearAlgebra:-Modular:-Determinant);

# doublePlusOne(X, k)
# 
# Input: X, k positive.
#
# Output: X^(2^(k+1) - 1).
#
doublePlusOne := proc(X, k):
  return X^(2^(k+1) - 1);
end proc;

# DoublePlusOneLift(A, X, n, k)
#
# Input:
#   A, n x n matrix in Z, nonsingular.
#   X, gcd(x, det(A)) = 1.
#   n, positive integer.
#   k, positive integer.
#
# Pre:
#   X ≥ max(10000, 3.61n^2||A||), Where ||A|| denotes the maximum absolute
#   entry.
#
# Output:
#   A_0, R_0, ..., R_(k), M_0, ..., M_(k-1):
#     A_0 = Rem(A^-1, X)
#     X_i = X^{2^(i+1) - 1}
#     B_(i+1) = B_(i) (I + R_i X_i) + M_i X_i^2  = A^-1 mod X_(i+1)
#
# NOTE: We append an extra term R_(k) which is used in our IntCertificate().
#
DoublePlusOneLift := proc(A, X, n, k)
  local A0, R, M, i, RSquared;
  ASSERT(X >= max(10000, 3.61 * n*n * max(A)));
  A0 := mods(ModInverse(X, A), X);
  R := Array(0..k);
  M := Array(0..k-1);
  R[0] := (Identity(2, n, integer) - Multiply(A, A0)) / X;

  for i from 0 to (k-1) do
    RSquared := MatrixPower(R[i], 2);
    M[i] := mods(ModMultiply(X, A0, RSquared), X);
    R[i+1] := (RSquared - Multiply(A, M[i])) / X;
  end do;

  return A0, R, M;
end proc;

# ComputeDoublePlusOneLift(n, k, X, A0, R, M, F)
#
# Input:
#   n, positive integer.
#   k, positive integer.
#   X, positive integer with, X ≥ max(10000, 3.61n^2||A||).
#      Where ||A|| denotes the maximum absolute entry.
#   A0, R, M, Computed from DoublePlusOneLift(A, X, n, k).
#   F, (n x m) matrix over Z.
#
# Output:
#   Rem(A^-1F, X_k), where X_i = X^{2^(i+1) - 1}
#
ComputeDoublePlusOneLift := proc(n, k, X, A0, R, M, F := Identity(2, n, integer))
  local p, i, Xk, Xi, E, m, Fc, _;
  _, m := Dimension(F);
  E := Matrix(n, m, 0);
  Xk := doublePlusOne(X, k);
  Xi := Xk;
  Fc := Matrix(F);

  for i from (k-1) to 0 do
    Xi := iquo(Xi, X);
    AddMultiple(Xk, Xi, E, ModMultiply(Xk, M[i], Fc), E);
    Xi := sqrt(Xi);
    AddMultiple(Xk, Xi, Fc, ModMultiply(Xk, R[i], Fc), Fc);
  end do;

  AddMultiple(Xk, 1, E, ModMultiply(Xk, A0, Fc), E);
  return mods(E, Xk);
end proc;

# SpecialSolve(A, B, d, n, m, p)
#
# Input:
#   A, n x n matrix over Z, nonsingular.
#   B, n x m matrix over Z for some m.
#   d, positive integer.
#   n, positive integer.
#   m, positive integer.
#   p, positive prime relatively prime to det(A).
#
# Output:
#   Rem(A^-1 B, 2^d)
#
SpecialSolve := proc(A, B, d, n, m, p)
  local mx, X, k, A0, R, M;
  mx := max(abs(A));
  X := p^(ceil(log(max(10000, 3.61 * n^2 * mx, p))));
  k := ceil(log2(d/log(X, p) + 1)) - 1;
  if (k < 1) then k := 1 end if;
  A0, R, M := DoublePlusOneLift(A, X, n, k);

  return mods(ComputeDoublePlusOneLift(n, k, X, A0, R, M, B), 2^d);
end proc;

# nearestPowerTwo(x)
#
# Input: x, positive integer
# 
# Output: y, where y is the smallest integer s.t. 2^y >= x.
#
nearestPowerTwo := proc(x)
  return 2^ceil(log2(x));
end proc;

# IntCertificate(s, A, B, n, m)
# 
# Input:
#   A, n x n matrix over Z, nonsingular.
#   s, positive integer.
#   B, n x m matrix over Z(s).
#   n, m, positive integer.
#
# Pre:
#   gcd(2, det(A)) = 1
#
# Output:
#   If sA^-1B is integral then return modp(sA^-1B, s) otherwise return false.
#
IntCertificate := proc(s, A, B, n, m):
  local bMx, aMx, X, h, A0, R, M, sRB, RBar, l, F;
  bMx := max(abs(B));
  aMx := max(abs(A));
  X := 2^(ceil(log2(max(10000, 3.61 * n^2 * aMx))));
  h := max(ceil(log2(log(2*s*n^(n/2)*aMx^(n-1)*bMx, X)/2+1)), 1);
  # First compute A^-1 = D + A^-1 Rbar 2^h for 2^(2^(h+1) -1) > 2*sn^n/2 ||A||^(n-1)||B||
  A0, R, M := DoublePlusOneLift(A, X, n, h);
  #RBar := (R[h] * X) + Multiply(A, M[h-1]);
  RBar := R[h-1] . R[h-1];
  sRB := s * Multiply(RBar, B);
  # Compute some l s.t. 2^l > 2n ||A||(0.6sn||B||)
  l := ceil(log2(1.2 * s * n^2 * aMx * bMx + 1));
  F := SpecialSolve(A, sRB, l, n, m);
  print(F);
  print(LinearAlgebra:-MatrixInverse(A).sRB);
  print(0.6 * s * n * bMx);
  # NOTE: F = sA^-1 Rbar B iff ||F|| < 0.6sn||B||
  if (max(abs(F)) < 0.6 * s * n * bMx) then
    modp(iquo(doublePlusOne(X, h), X) * F, s);
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
UniCert := proc(A, n)
  local mx, e, X, expr, k, A0, R, M, zero, det;

  if (ModDeterminant(2, A) = 0) then return false end if;

  zero := Matrix(n, n);
  mx := max(abs(A));
  e := ceil(log2(3.61 * n^2 * mx));
  X := max(2^e, 10000);
  expr := n^((n-1)/2) * mx^(n-1) / (n^2 * mx);

  if (X > convert(expr, float)) then
    k := 1;
  else
    k := ceil(log2(log(expr, X) + 2) - 1);
  end if;

  A0, R, M := DoublePlusOneLift(A, X, n, k);

  for i from k-1 to 0 do
    if (Equal(zero, R[i])) then return true end if;
  end do;
  return false;
end proc; 
