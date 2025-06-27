read("helpers/basic.mpl");
read("config.mpl");

macro(Map=Threads:-Map);

#
# MatrixRadix(A, n, m, X)
#
# Input:
#   A, n x m integer matrix.
#   X, positive modulus
#
# Output: Abar, q 
#   Abar := [A_0 | A_1 | A_2 | ... | A_q] s.t. A = A_0 + A_1X + ... + A_qX^q
#
MatrixRadix := proc(A, n, m, X)
  local Apos, Aneg, q, ret;
  q := MinLargestPower(max(abs(A)), X);
  res := Matrix(n, m * q);

  MatrixRadixHelper := proc(i)
    for j from 1 to m do
      XExpansion := convert(A[i,j], base, X);
      l := nops(XExpansion);
      for k from 1 to l do
	res[i, j + (k-1)*m] := XExpansion[k];
      end do;
    end do;
  end proc;

  Map[tasksize=TASKSIZE](MatrixRadixHelper, [seq(1..n)]);

  return res, q;
end proc;

#
# InverseRadix(Abar, n, m, q, X)
#
# Input:
#   Abar, n x q*m integer matrix.
#   X, positive modulus
#
# Pre:
#   Abar := [A_0 | A_1 | ... | A_q]
#
# Output: A 
#   A = A_0 + A_1X + ... + A_qX^q
#
InverseRadix := proc(Abar, n, m, X, q)
  return InverseRadixHelper(Abar, n, m, X, 1, q);
end proc;

InverseRadixHelper := proc(Abar, n, m, X, l, r)
  local mid;

  if (l = r) then return Abar[1..n, (l-1)*m+1..l*m] end if;

  mid := iquo(l + r, 2);
  return InverseRadixHelper(Abar, n, m, X, l, mid) + X^(mid+1-l) * InverseRadixHelper(Abar, n, m, X, mid+1, r);
end proc;

#
# RadixMult(A, B, l, n, m, a, p)
#
# Input:
#   A, l x n integer matrix.
#   B, n x m integer matrix.
#
# Output: mods(A . a * B, p)
#
RadixMult := proc(A, B, l, n, m, a := 1, p := 0)
  local X, BRadix, q, aABRadix, aAB;

  X := max(2^MinLargestPower(max(abs(A)), 2), 2^64);
  BRadix, q := MatrixRadix(B, n, m, X);

  if (p <> 0) then
    #aABRadix := imlMult(A, BRadix, a, p);
    aABRadix := fmpzMult(A, BRadix, a, p);
  else
    #aABRadix := imlMult(A, BRadix, a);
    aABRadix := fmpzMult(A, BRadix, a);
  end if;

  aAB := InverseRadix(aABRadix, l, m, X, q);

  return aAB;
end proc;
