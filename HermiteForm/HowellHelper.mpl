read("helpers/basic.mpl");
read("SmithMassager/SmithMassager.mpl");
read("HermiteForm/hcol.mpl");
read("helpers/mod.mpl");
read("HermiteForm/naiveScaledMatVecProd.mpl");
#read("HermiteForm/scaledMatVecProd.mpl");
macro(Identity=LinearAlgebra:-Modular:-Identity);
macro(Multiply=LinearAlgebra:-Multiply);
macro(ModMultiply=LinearAlgebra:-Modular:-Multiply);
macro(ModAddMultiple=LinearAlgebra:-Modular:-AddMultiple);
kernelopts(assertlevel=1);

# HermiteDiagonals(n, M, S, SNum)
#
# Input:
#   S, SNum, compact representation of some n x n Smith form of A.
#   M, compact cmod form of n x n reduced smith massager of A.
#
# Output: h
#   h, n vector, where h_1, ..., h_n corresponds to the Hermite diagonal of A.
#
HermiteDiagonals := proc(n, M, S, SNum)
  local k, Hi, MM;
  MM := Copy(M);
  diagonals := Vector(n, 1);

  l := 1;
  for i from SNum to 1 by -1 do
    Hi, idices, numIdx := compactHcol(MM[1..n, i], S[i]);
    for j from 1 to numIdx do
      diagonals[idices[j]] *= Hi[idices[j], j];
      if (i <> 1 and Hi[idices[j], j] >= S[i-1]) then
        Hi[1..idices[j], j] := modp(Hi[1..idices[j], j], S[i-1]);
      end if;
    end do;
    MM := cmodCompactHermiteMult(Hi, idices, numIdx, MM, n, S, SNum - l);
    ++l;
  end do;
  return diagonals;
end proc;

# SpecialHowellTransform(n, M, SEntry, Snum, h, p)
#
# Input:
#   SEntry, SNum, compact representation of some n x n Smith form of A.
#   M, compact form of a reduced smith massager of A, i.e. AM = 0 cmod S.
#   h, n vector containing the Hermite diagonal entries of A.
#   p, a prime with gcd(p, s) = 1 and p in O(loglog(s))
#   n, a positive integer.
#
# Output: U in compact form wrt to rmod S.
#   U, n x n s.t. M sS^-1 U is a Howell form of M sS^-1 in mod s.
#
SpecialHowellTransform := proc(n, M, SEntry, SNum, h, p)
  local a, i, s, j, tmp;

  if (SNum = 0) then
    return Matrix(0, n);
  end if;

  s := SEntry[-1];
  k := n - SNum;
  cols := n + SNum;
  U := Matrix(SNum, cols, (i, j) -> if (i = j - SNum - k) then 1 else 0 end if);

  for i from 0 to n-1 do
    ePivot := cols-i;
    sPivot := ePivot - SNum + 1;
    if (h[n-i] <> 1) then  
      # Calculate a^T := M[n-i] S* U[1..-1, sPivot..ePivot].
      a := CompactScaledMatVecProd(SNum, ModTranspose(s, U[1..-1, sPivot..ePivot]), SEntry, SNum, ModTranspose(s, M[n-i]), h[n-i], p);
      c := ExtendedStab(SNum, a, h[n-i], SNum);
      # Calculate U := UC_i as defined in the paper.
      U[1..-1, ePivot] := rmodCompactMult(U[1..-1, sPivot..ePivot], n, convert(c, 'Matrix'), 1, SEntry, SNum);
      d, u, g := reScale(Multiply(Transpose(a), c) * iquo(s, h[n-i]), s);
      e := modp(u + d * s / g, s);
      U[1..-1, ePivot] := compactRmodScalarMultiply(U[1..-1, ePivot], e, SEntry, SNum);

      # Calculate U := UW_i as defined in the paper.
      for j from 1 to SNum - 1 do
        tmp := compactRmodScalarMultiply(U[1..-1, ePivot], modp(-a[-j - 1], s), SEntry, SNum);
        U[1..-1, ePivot - j] := compactRmodAdd(U[1..-1, ePivot - j], tmp, SEntry, SNum);
      end do;
    end if;

    tmp := compactRmodScalarMultiply(U[1..-1, ePivot], modp(h[n - i], s), SEntry, SNum);
    U[1..-1, sPivot - 1] := compactRmodAdd(U[1..-1, sPivot - 1], tmp, SEntry, SNum);
  end do;

  return U[1..-1, -n...-1];
end proc;

# HermiteViaHowell(n, M, SEntry, SNum, U, h, p)
#
# Input:
#   SEntry, SNum, compact representation of some n x n Smith form of A.
#   M, n x SNum, compact form of a reduced smith massager of A, i.e. AM = 0 cmod S.
#   U, SNum x n in compact form.
#      With T := MS*U, where T is the Howell form of MS* in Z(s). 
#      Note: T is also Howell form of sA^-1, since vA^-1 is integral 
#            iff vMS^-1 is integral.
#   h, n vector containing the Hermite diagonal entries of A.
#   p, a prime with gcd(p, s) = 1 and p in O(loglog(s))
#   n, a positive integer.
#
# Output: H
#   H, n x n Hermite form of A.
#
HermiteViaHowell := proc(n, M, SEntry, SNum, U, h, p)
  local j;
  s := SEntry[-1];
  H := Matrix(n, n, (i,j) -> if (i = j) then h[i] else 0 end if);
  if (s = 1) then return H end if;

  for j from 1 to n do:
    if (h[j] = 1) then next; end if;
    v := CompactScaledMatVecProd(n, M, SEntry, SNum, U[1..-1, j], h[j], p);

    v := modp(-v, h[j]);

    H[1..j-1, j] := v[1..j-1];
    M[1..-1,1..-1] := cmodCompactHermiteMult1(H[1..-1, j..j], Vector([j]), 1, M, n, SEntry, SNum);
  end do;

  return H;
end proc;

# HermiteMassager := proc(n, A)
#
# Input:
#   A, n x n nonsingular integer matrix.
#
# Output: H
#   H, n x n upper triangular Hermite form of A, i.e. H = XA for some
#      unimodular X.
#
HermiteMassager := proc(n, A)
  local ret, M, S, _, s, h;

  ret := SmithMassager(A, n);

  if (ret = false) then return false; end if;

  _, M, _, S := ret;
  SS := DiagonalMatrix(S);

  return HermiteMassagerHelper(n, M, SS);
end proc;

# HermiteMassagerHelper := proc(n, A, M, S)
#
# Input:
#   M, S, the smith massager and smith form of some non-singular n x n matrix A.
#
# Output: H
#   H, n x n upper triangular Hermite form of A, i.e. H = XA for some
#      unimodular X.
#
HermiteMassagerHelper := proc(n, M, S)
  local ret, s, h;

  SnumNonTrivial := 0;
  SnonTrivialEntry := [];
  for i from 1 to n do
    if (S[i,i] <> 1)  then
      ++SnumNonTrivial;
      SnonTrivialEntry := [op(SnonTrivialEntry), S[i,i]];
    end if;
  end do;
  SnonTrivialEntry := convert(SnonTrivialEntry, 'Vector');

  s := S[n,n];
  p := 2;
  while (irem(s, p) = 0) do
    # TODO: Change this to nextprime(p) and test it.
    p := nextprime(s);
  end do;

  k := n - SnumNonTrivial;
  MM := Copy(M[1..-1, k+1..n]);

  h := HermiteDiagonals(n, MM, SnonTrivialEntry, SnumNonTrivial);
  U := SpecialHowellTransform(n, MM, SnonTrivialEntry, SnumNonTrivial, h, p);
  return HermiteViaHowell(n, MM, SnonTrivialEntry, SnumNonTrivial, U, h, p);
end proc;

A := Matrix([[1,2],[5,4]]);
H := HermiteMassager(2, A);
ASSERT(Equal(HermiteForm(A), H));

A := Matrix([[1,2,3],[1,2,10],[5,4,8]]);
H := HermiteMassager(3, A);
ASSERT(Equal(HermiteForm(A), H));
HermiteForm(A);

A := Matrix([[1,9,3],[4,3,10],[5,4,8]]);
H := HermiteMassager(3, A);
ASSERT(Equal(HermiteForm(A), H));

A := Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 4]]);
H := HermiteMassager(3, A);
ASSERT(Equal(HermiteForm(A), H));
