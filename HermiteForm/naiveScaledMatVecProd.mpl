read("helpers/diagHelpers.mpl");
macro(ScalarMultiply=LinearAlgebra:-ScalarMultiply);
macro(ModAddMultiple=LinearAlgebra:-Modular:-AddMultiple);
macro(ModMultiply=LinearAlgebra:-Modular:-Multiply);

# ScaledMatVecProd(n, M, S, u, h, p)
#
# Input:
#   S, n x n smith form of a nonsingular integer matrix, let s := S[n,n].
#   M, n x n integer matrix s.t. M = cmod(M, S)
#   u, n vector s.t. u = rmod(u, S)
#   h, positive divisor of s s.t. (s/h)^-1 M sS^-1 u is integer.
#   p, odd prime with gcd(p, s) = 1 and p in O(loglog det(S)).
#
# NOTE: p is not used in this naive approach.
#
# Ouptut: v
#   v, n vector in mod h s.t.
#   s/h v = M sS^-1 u mod s
#
ScaledMatVecProd := proc(n, M, S, u, h, p := 0)
  local ret, s;
  s := S[n,n];

  return (modp(M . (s * diagInv(S, n)) . u, s)) * h/s;
end proc;

# CompactScaledMatVecProd(n, M, SEntry, SNum, u, h, p)
#
# Input:
#   SEntry, SNum, compact representation of a diagonal ascending n x n matrix, S.
#   M, n x SNum integer matrix in compact form wrt cmod S.
#   u, Snum vector in compact form wrt rmod S.
#   h, positive divisor of s s.t. (s/h)^-1 M sS^-1 u is integer.
#   p, odd prime with gcd(p, s) = 1 and p in O(loglog det(S)).
#
# NOTE: p is not used in this naive approach.
#
# Ouptut: v
#   v, n vector in mod h s.t.
#   s/h v = M sS^-1 u mod s
#
CompactScaledMatVecProd := proc(n, M, SEntry, SNum, u, h, p := 0)
  local s, res;
  s := SEntry[-1];
  res := Vector(n);

  for i to SNum do
    ModAddMultiple(s, modp(iquo(s * u[i], SEntry[i]), s), res, M, 1..-1, i, res);
  end do;
  return h * res / s;
end proc;
