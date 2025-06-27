read "SmithMassager/indexMassager.mpl";
read "helpers/diagHelpers.mpl";
read "helpers/largestInvariantFactor.mpl";
read "helpers/mod.mpl";
read "helpers/radix.mpl";
read "helpers/naiveCert.mpl";
read "helpers/fastUniCert.mpl";
read "macros.mpl";

macro(Multiply=LinearAlgebra:-Multiply);
macro(ModMult=LinearAlgebra:-Modular:-Multiply);
macro(Add=LinearAlgebra:-Add);
macro(AddMultiple=LinearAlgebra:-Modular:-AddMultiple);
macro(Determinant=LinearAlgebra:-Determinant);
macro(SmithForm=LinearAlgebra:-SmithForm);
macro(DeleteColumn=LinearAlgebra:-DeleteColumn);
macro(DeleteRow=LinearAlgebra:-DeleteRow);
macro(Inverse=LinearAlgebra:-Modular:-Inverse);
macro(AreCoprime=NumberTheory:-AreCoprime);

# SmithMassager(A, n)
#
# Input:
#  A, n x n integer matrix, nonsingular.
#  n, positive integer.
#
# Output: U, M, T, S
#   Tries to find a (0, n)-index Smith Massager for A or false.
#   U, M, T, S, n x n integer matrix.
#   S, in smith form of A.
#   M, with property AM = 0 cmod S
#   T, unit upper triangular with T + UM = 0 cmod S.
#
SmithMassager := proc(A, n)
  local s, B, U, M, T, S, sumR, i, m, r, Up, Mp, Tp, Sp, UpMp, SpInv, sLargest, k,
  MpSpInv, NegUpM, a, b, l;
  startDim := 5;
  k := 5;
  sLargest, Q := LargestInvariantFactor(A, n, startDim+k);
  s := sLargest;
  B := Identity(2, 2 * n, integer);
  B[1..n, 1..n] := A;
  U := Matrix(0, n);
  M := Matrix(n, 0);
  T := Matrix(0, 0);
  S := Vector(0);
  l := 2*n;
  sumR := 0;

  # Pre: 
  # At the beginning of the loop: (U, M, T, S) is a (0, m) index
  # massager for diag(A, I_n).
  # At the end of the loop (except the last iteartion): (U, M, T, S) 
  # becomes a (0, m+r) index massager with high probability for diag(A, I_n).
  # At the last iteration: we will have a (0, n) index massager for diag(A, I_n).
  for i while sumR < n do
  #i := 0;
  #while (sumR <> n) do
  #  ++i;

    m := sumR;
    #r := min(2^(i-1), n - m, 256);
    #r := min(2^(i-1), n - m, 64);
    #r := max(8^(i-1), startDim);
    r := min(max(ceil((log(s) - 25) * 1.5), startDim), 500);
    # 500: 453, 250: 320, 125: 274, 50: 242 
    #r := max(500, startDim);
    r := max(2^(i-1), r);
    if (i = 1) then r := startDim end if;
    #r := max(r, log2(s) - 25);
    r := min(r, n - m);
    #R := highOrderResidue(B);
    #if (s <> 1) then 
    #  #r := min(ceil(n*log2(max(abs(R), 2))/log2(s)), n - m);
    #  r := min(ceil(n*(1 - log2(s)/log2(sLargest))) + 1, n - m);
    #else
    #  r := n-m;
    #end if;

    dt := getDT(s);
    if (i <> 1 and dt = float[8]) then
      r := n - m;
    end if;
    if (i = 1) then Q := Q[1..-1, 1..r+k] end if;
    print(r);

    # Compute an (m,r) index massager for B.
    # (rxn),(nxr),(rxr),(rxr) matrices.
    print("Calling Index massager", i);
    ret := IndexMassager(B, n, m, r, s, 1/(8*(r)), Q, k);
    Q := false;
    k := 0;
    if (ret = false) then return false end if;
    Up, Mp, Tp, Sp := ret;

    # In case dt was float 8 we need to recalculate T.
    if (dt = float[8]) then
      # Negative is needed to cancel out stuff later on.
      Tp := Matrix(r, r);
      Tp := -LeftSparseMult(Up, Mp, r, n, r);
      # Clear the lower triangular.
      for j to r do for k from j to r do Tp[k, j] := 0 od od;
      # Replace 0 by 1.
      for k to r do Tp[k,k] := 1 od;
    end if;
    

    print("Finish Index massager", i);
    #Up := rmod(Up, Sp);
    #Mp := cmod(Mp, Sp);
    #UpMp := Multiply(Up, Mp);
    #UpMp := imlMult(Up, Mp);
    #Tp := lowerCmod(-UpMp, Sp);
    #ReplaceDiagonal(Tp, r, 0, 1);
    #if (not Equal(Tp, ReplaceDiagonal(lowerCmod(-UpMp, Sp), r, 0, 1))) then ASSERT(false); end if;

    # Apply the (m,r)-index massager on B and update (U, M, T, S) to become 
    # (0, r+m)-index massager for diag(A, I_n).
    # Pre:
    # B looks like:
    # [A       A M S^-1   ]
    # [   I               ]
    # [U    (T+ U M) S^-1 ]
    #
    # i.e. B has a (0, sumR)-index smith massager (U, M, T, S).
    #
    # Post:
    # B looks like:
    # [A         A Mp Sp^-1               A M S^-1        ]
    # [   I                                               ]
    # [Up    (Up Mp + Tp) Sp^-1                           ]
    # [U          U Mp Sp^-1          (T + U M)S^-1       ]

    P := max(abs(A), abs(U), abs(Up)) * n * 2 + 1;
    p := 1709;
    while (not AreCoprime(p, Sp[-1])) do p := nextprime(p); end do;
    while (p < P) do p *= p end do;
    B[l-sumR-r+1..l-sumR, 1..n] := Up;
    SpInv := ModDiagInv(p, Sp, r);
    MpSpInv := ModDiagMult(p, Mp, SpInv, n, r);
    Update := imlMult(Concatenate(1, A, Concatenate(1, Up, U)), MpSpInv);
    #Update := fmpzMult(Concatenate(1, A, Concatenate(1, Up, U)), MpSpInv);
    #Update := RadixMult(Concatenate(1, A, Concatenate(1, Up, U)), MpSpInv, n+sumR+r, n, r);
    TpSpInv := ModDiagMult(p, Tp, SpInv, r, r);
    B[1..n, l-sumR-r+1..l-sumR] := mods(Update[1..n, 1..r], p);
    B[l-sumR+1..l, l-sumR-r+1..l-sumR] := mods(Update[n+r+1..-1, 1..r], p);
    B[l-sumR-r+1..l-sumR, l-sumR-r+1..l-sumR] := mods(Update[n+1..n+r, 1..r] + TpSpInv, p);

    # Update U, M, T, S to be (0, sumR+r)-index massager for A.
    T := ConcatenateDiag(Tp, T, r, r, sumR, sumR);
    if (sumR <> 0) then
      #NegUpM := RadixMult(-Up, M, r, n, sumR);
      #NegUpM := imlMult(-Up, M);
      NegUpM := LeftSparseMult(-Up, M, r, n, sumR);
      #NegUpM := -Up . M;
      T[1..r, -sumR..-1] := NegUpM;
    end if;
    U := Concatenate(1, Up, U);
    M := Concatenate(2, Mp, M);
    S := ConcatenateDiagVector(Sp, S, r, sumR);
    sumR := sumR + r;
    s := S[1];
  end do;

  # Binary search for the index 1
  l := 0;
  r := n;
  while l < r do
    m := ceil((l+r)/2);
    if S[m] <> 1 then
      r := m-1;
    else
      l := m;
    end if;
  end do;
  if (l <> 0) then
   E := DeleteColumn(DeleteRow(B, n+1..n+l), n+1..n+l);
   correct := UniCert(E, 2*n-l);
  else
   correct := UniCert(B, 2*n);
  end if;
  #start := time[real]();
  #correct := UniCert(B, 2*n);
  #endTime := time[real]();
  #print("UniCert:", endTime - start);
  if (correct) then
    userinfo(1, SmithMassager, B);
    return U, M, T, S
  else
    userinfo(1, SmithMassager, B);
    userinfo(1, SmithMassager, S);
    userinfo(1, SmithMassager, Determinant(B));
    return false;
  end if;
end proc;
