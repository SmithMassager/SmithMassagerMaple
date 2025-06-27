read "SmithMassager/projectionBasis.mpl";
macro(DiagonalMatrix=LinearAlgebra:-DiagonalMatrix);
kernelopts(assertlevel=1);

# trivialMassager(n, m, r)
# Input: 
#   n, m, r positive integer.
#
# Output: 0_(r x n), 0_(r x n), I_r, I_r
#
trivialMassager := proc(n, m, r):
  return Matrix(r, n), Matrix(n, r), Identity(2, r, integer), Vector(r, (i) -> 1);
end proc;

# IndexMassager(B, n, m, r, s, eps)
#
# Input: 
#   B, (2n x 2n) integer matrix with shape:
#         [A              *]
#         [    I_{n-m}     ]
#         [*              *]
#     Last n rows and columns of B^{-1} are integral.
#   n, m, r, positive integer with n+m <= r
#   s, is positive integer multiple of the largest invariant factor of B.
#   eps, in (0,1).
#   Q, an optinal projection matrix.
#
# Output: U, M, T, S
#   U, r x n integer Mod s matrix.
#   M, n x r integer Mod s matrix.
#   T, r x r integer matrix.
#   S, r x r integer matrix, nonsingular and in smith form, with S_rr a divisor
#   of s.
#
#   Index-(m,r) Smith massager (U, M, T, S) for B with T == I_r.
#   With at least 1 - eps probability a maximal index-(m,r) smith massager for
#   B is returned. i.e. S is compromised of the r largest invariant factors of
#   B.
#
IndexMassager := proc(B, n, m, r, s, eps, Q := false, kk := 0):
  local F, U, M, k, S, J1, J2, J, P;
  global timesIMLSolve;

  ASSERT(m + r <= n);

  if (s = 1) then
    return trivialMassager(n, m, r);
  end if;

  #k := ceil(log2(1/eps)) + 1;
  k := 5;
  if (Q = false) then
    J1 := RandomMatrix(n, r+k, generator = 0 .. (s-1));
    #J1 := RandomMatrix(n, r+k, generator = 0 ..1);
    J2 := Matrix(n,r+k);
    J  := Concatenate(1, J1, J2);
    start := time[real]();
    P := SpecialIntCertificate(s, B, J, n, r+k, m);
    endTime := time[real]();
    timesIMLSolve := [op(timesIMLSolve), endTime-start];

    if (P = false) then
      WARNING(sprintf("sB^{-1}J is not integral"));
      return false;
      #return trivialMassager(n, m, r);
    end if;

    # P2 should be all zero matrix if not return the trivial index massager.
    P2 := SubMatrix(P, n+1..2*n, 1..r+k);
    P1 := SubMatrix(P, 1..n, 1..r+k);
    if (not Equal(P2, Matrix(n, r+k))) then
      WARNING(sprintf("P2 not zero"));
      return false;
    end if;
  else
    k := kk;
    P1 := Q;
  end if;

  #(S, U, M) := computeProjBasis(P, n, r+k, s);
  ret := computeProjBasis(P1, n, r, s, k);
  if (ret = false) then 
    WARNING(sprintf("computeProjBasis returned false"));
    return false;
  end if;
  S, U, M, T := ret;

  # TODO: Change this to a vector.
  #S := DiagonalMatrix(S);

  # Extract the r largest invariant factors and the corresponding multiplier.
  #S := SubMatrix(S, -r..-1, -r..-1);
  #U := SubMatrix(-U, -r..-1, 1..-1);
  #M := SubMatrix(-M, 1..-1, -r..-1);
  #U := -U;
  M := cmod(M, S, true);
  #ReplaceDiagonal(T, r, 0, 1);

  return U, M, T, S 
  #return U, M, Identity(2, r, integer), S 
end proc;
