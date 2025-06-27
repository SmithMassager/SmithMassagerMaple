read("helpers/diagHelpers.mpl");
macro(ScalarMultiply=LinearAlgebra:-ScalarMultiply);
macro(ModAddMultiple=LinearAlgebra:-Modular:-AddMultiple);
macro(ModMultiply=LinearAlgebra:-Modular:-Multiply);

# smallestPower(thres, p)
#
# Input:
#   p, thres, positive integer.
#
# Output: X
#   Where X is the smallest power of p > thres.
#
smallestPower := proc(thres, p)
  local X, i;

  X := p;
  i := 1;
  while(X <= thres) do
    X *= p;
    ++i
  end do;

  return X, i;
end proc;

# modScaleVector(w, u, X, si, Y, k)
#
# Input:
#   w, n vector.
#   u, scalar.
#   si, k, X, Y, positive integer.
#
# Pre:
#   All the undefined variables are referring to the variable in parent
#   function: CompactScaledMatVecProd.
#   X^k > s
#   X > 2h log det S
#   Y > X^2, relatively prime to s.
#
# Output: v
#   w' := (w1, w2, ..., wn)^T where wi := (wi0, wi1, .., wik-1)^T 
#   with wi := wi0 + wi1 X + ... + wik-1 X^k-1 mod X^k.
#   u' := u * (X^0, X^1, ..., X^k-1)^T mod s
#   v := 1/si (w' u' ) mod Y
#
modScaleVector := proc(w, u, X, si, Y, k)
  local n, res, a, q, r, tmpQ, tmpR, mxIdx, j;

  n := Dimensions(w);
  res := Vector(n);
  a := u;
  q := Vector(n, (i) -> w[i]);
  r := Vector(n);

  for j to k do
    # Compute a linear lift of w wrt to X.
    for i to n do
      tmpQ := iquo(q[i], X, 'tmpR');
      q[i] := tmpQ;
      r[i] := tmpR;
    end do;

    ModAddMultiple(Y, a, res, r, res);
    a := modp(a * X, si);
  end do;

  ModMultiply(Y, modp(1/si, Y), res, res);
  return res;
end proc;

# CompactScaledMatVecProd(n, M, SEntry, SNum, u, h, p)
#
# Input:
#   SEntry, SNum, compact representation of a diagonal ascending n x n matrix, S.
#                 let s denote the largest invariant factor.
#   M, n x SNum integer matrix in compact form wrt cmod S.
#   u, Snum vector in compact form wrt rmod S.
#   h, positive divisor of s s.t. (s/h)^-1 M sS^-1 u is integer.
#   p, odd prime with gcd(p, s) = 1 and p in O(loglog det(S)).
#
#
# Ouptut: v
#   v, n vector in mod h s.t.
#   s/h v = M sS^-1 u mod s
#
CompactScaledMatVecProd := proc(n, M, SEntry, SNum, u, h, p)
  local ret, s, detS, X, Y, i, k, _, logDetS;
  s := SEntry[-1];
  res := Vector(n);

  detS := 1;
  for i to SNum do
    detS *= SEntry[i];
  end do;

  logDetS := convert(log2(detS), 'float');
  X, _ := smallestPower(2 * h * logDetS, 2);
  Y, _ := smallestPower(X * X, p);

  for i to SNum do
    _, k := smallestPower(SEntry[i], X);
    tmp := modScaleVector(M[1..-1, i], u[i], X, SEntry[i], Y, k);
    ModAddMultiple(Y, 1, res, tmp, res);
  end do;

  ModMultiply(Y, modp(h, Y), res, res);
  return res;
end proc;
