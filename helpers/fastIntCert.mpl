read "helpers/partialLinearization.mpl";
read "helpers/radix.mpl";
read "helpers/time.mpl";
read "helpers/LinearSolver.mpl";
read "helpers/basic.mpl";
read "macros.mpl";


macro(Multiply=LinearAlgebra:-Multiply);
macro(ModMultiply=LinearAlgebra:-Modular:-Multiply);
macro(Inverse=LinearAlgebra:-Modular:-Inverse);
macro(Mod=LinearAlgebra:-Modular:-Mod);

# IntCertificate(s, A, B, n, m)
# 
# Input:
#   A, n x n matrix over Z, nonsingular.
#   s, positive integer.
#   B, n x m matrix over Z(s).
#   n, m, positive integer.
#
# Output:
#   If sA^-1B is integral then return modp(sA^-1B, s) otherwise return false.
#
fastIntCertificate := proc(s, A, B, n, m)
  global timesIMLSolve, IMLDim, timesHighOrder, timesRadix, SDim, ResDim, SRBDim;

  local res, d, Abar, nbar, _, Bbar, ret, l;
  if (not type(A, 'Matrix'(integer))) then
    return false;
  end if;
  
  #l := 2;
  thres := 1.2 * s * n^2 * max(abs(A)) * max(abs(B));
  #while (l < thres) do
  #  l := l * l;
  #end do;

  #Abar, nbar, _ := ColPL(A, n, n);
  Abar := A;
  nbar := n;
  start := time[real]();
  R := highOrderResidue(Abar);
  endTime := time[real]();
  timesHighOrder := [op(timesHighOrder), endTime-start];

  Bbar := Matrix(nbar, m);
  Bbar[1..n, 1..m] := B;
  start := time[real]();
  #sRB := RadixMult(R, Bbar, nbar, nbar, m, s, l);
  sRB := imlMult(R, Bbar, s);
  endTime := time[real]();
  timesRadix := [op(timesRadix), endTime-start];
  IMLDim := [op(IMLDim), m];
  start := time[real]();
  if (getDT(s) <> float[8]) then
    res, d := imlSolveHelper(Abar, sRB);
    print(res, d);
    if (d <> 1) then return false; end if;
  else
    res := ChineseLinearSolve(Abar, sRB, thres, s);
    if (res = false) then return false end if;
  end if;
  endTime := time[real]();
  timesIMLSolve := [op(timesIMLSolve), endTime-start];

  SDim := [op(SDim), convert(log2(s), float)];
  #ResDim := [op(ResDim), convert(log2(max(R)), float)];
  #SRBDim := [op(SRBDim), convert(log2(max(sRB)), float)];

  ret := Matrix(n, m);
  ret[1..n, 1..m] := res[1..n, 1..m];

  return modp(ret, s);
end proc;

# SpecialIntCertificate(s, A, B, n, m, r)
# 
# Input:
#   A, 2n x 2n matrix over Z, nonsingular.
#   s, positive integer.
#   B, 2n x m matrix over Z(s).
#   n, m, r, positive integer.
#
# Pre:
#   NOTE: Empty space means zero entries.
#
#   A = [ nxn    nxr ]
#       [            ]
#       [            ]
#       [ rxn    rxr ]
#
#   B = [ n x m ]
#       [       ]
#
# Output:
#   If sA^-1B is integral then return modp(sA^-1B, s) otherwise return false.
#
SpecialIntCertificate := proc(s, A, B, n, m, r)
  local res, d, Abar, nbar, _, Bbar, ret, l;
  if (not type(A, 'Matrix'(integer))) then
    return false;
  end if;

  AA := Matrix(n+r, n+r);
  AA[1..n, 1..n] := A[1..n, 1..n];
  if (r <> 0) then
    AA[-r..-1, 1..n] := A[-r..-1, 1..n];
    AA[1..n, -r..-1] := A[1..n, -r..-1];
    AA[-r..-1, -r..-1] := A[-r..-1, -r..-1];
  end if;

  BB := Matrix(n + r, m);
  BB[1..n, 1..m] := B[1..n, 1..m];
  res := IntCertificate(AA, BB, s);

  if (res = false) then return false end if;

  ret := Matrix(2*n, m);
  ret[1..n, 1..m] := res[1..n, 1..m];
  if (r <> 0) then
    ret[-r..-1, 1..m] := res[-r..-1, 1..m];
  end if;

  return ret;
end proc;
