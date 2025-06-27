macro(ScalarMultiply=LinearAlgebra:-ScalarMultiply);
macro(ModMultiply=LinearAlgebra:-Modular:-Multiply);
macro(Concatenate=ArrayTools:-Concatenate);
macro(ModAddMultiple=LinearAlgebra:-Modular:-AddMultiple);

# ReplaceDiagonal(A, n, x, y)
#
# Input:
#   A: n x n matrix
#   x, y: a number
#
# Post:
#   Every diagonal entries with value x replaced by y.
#
ReplaceDiagonal := proc(A, n, x, y)
  local i;
  for i from 1 to n do
    if (A[i,i] = x) then A[i,i] := y end if;
  end do;
  return A;
end proc;

# DiagInv(A, n)
#
# Input:
#   A, n x 1 vector reprenseting n x n diagonal matrix.
# 
# Output: A^-1
#
DiagInv := proc(A, n)
  return Vector(n, (i) -> 1/A[i]);
end proc;

# ModDiagInv(p, A, n)
#
# Input:
#   A, n x 1 vector reprsenting n x n diagonal matrix.
#   p, prime.
# 
# Output: A^-1
#
ModDiagInv := proc(p, A, n)
  return Vector(n, (i) -> modp(1/A[i], p));
end proc;


# ConcatenateDiag(A, B, n, m, p, q)
#
# Input:
#   A: n x n matrix.
#   B: p x p matrix.
#
# pre: p = q, n = m.
#
# Output: (n+p) x (n+p) matrix:
#   [A ]
#   [ B]
#
ConcatenateDiag := proc(A, B, n, m, p, q)
  local r, l;

  ASSERT(n = m);
  ASSERT(p = q);

  if (n = 0) then
    return B;
  elif (p = 0) then
    return A;
  end if;

  l := Concatenate(1, A, Matrix(p, n));
  r := Concatenate(1, Matrix(n, p), B);
  return Concatenate(2, l, r);
end proc;

# ConcatenateDiagVector(A, B, n, p)
#
# Input:
#   A: n x 1 matrix.
#   B: p x 1 matrix.
#
# Output: (n+p) x (n+p) matrix:
#   [A ]
#   [ B]
#
ConcatenateDiagVector := proc(A, B, n, p)
  local res;

  if (n = 0) then
    return B;
  elif (p = 0) then
    return A;
  end if;

  res := Vector(n+p);
  res[1..n] := A;
  res[n+1..-1] := B;
  return res;
end proc;


# DiagMult(A, B, m, n)
#
# Input:
#   A, m x n matrix.
#   B, n x 1 vector representing n x n diagonal matrix.
#   n, positive integer.
#
# Output: E
#   E, E = A * B.
DiagMult := proc(A, B, m, n):
  local E, i;
  E := Matrix(m, n);

  for i from 1 to n do
    E[1..-1,i] := ScalarMultiply(A[1..-1, i], B[i]);
  end do;

  return E;
end proc;

# ModDiagMult(p, A, B, n)
#
# Input:
#   A, m x n matrix.
#   B, n x 1 diagonal matrix.
#   n, positive integer.
#   p, positive integer.
#
# Output: E
#   E, E = A * B mod s.
ModDiagMult := proc(p, A, B, m, n):
  return modp(DiagMult(modp(A, p), modp(B, p), m, n), p);
end proc;

# DiagDet(A, n)
#
# Input:
#   A, n x 1 diagonal matrix.
# 
# Output: Det(A)
#
DiagDet := proc(A, n)
  local i, ret;
  ret := 1;
  for i from 1 to n do
    ret .= A[i];
  end do;

  return ret;
end proc;
