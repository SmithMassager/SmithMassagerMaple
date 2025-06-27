read "SmithMassager/SNFMod.mpl";
read "SmithMassager/reverseSmithForm.mpl";
read "helpers/basic.mpl";
read "helpers/diagHelpers.mpl";
read "helpers/fastIntCert.mpl";
read "helpers/fastUniCert.mpl";
read "helpers/multiplyHelper.mpl";
read "helpers/naiveCert.mpl";
read "helpers/radix.mpl";
read "helpers/time.mpl";
read "macros.mpl";

macro(AddMultiple=LinearAlgebra:-Modular:-AddMultiple);
macro(Concatenate=ArrayTools:-Concatenate);
macro(DiagonalMatrix=LinearAlgebra:-DiagonalMatrix);
macro(Equal=LinearAlgebra:-Equal);
macro(Identity=LinearAlgebra:-Modular:-Identity);
macro(MatrixInverse=LinearAlgebra:-MatrixInverse);
macro(Mod=LinearAlgebra:-Modular:-Mod);
macro(Multiply=LinearAlgebra:-Modular:-Multiply);
macro(RandomMatrix=LinearAlgebra:-RandomMatrix);
macro(Random=LinearAlgebra:-Modular:-Random);
macro(SubMatrix=LinearAlgebra:-SubMatrix);

timeBase := 0;
timeExtract := 0;
timeVector := 0;
F := Matrix(0, 0);

# SNF(A, n, r, k, s)
# 
# Input:
#  A, n x (r + k) Z(s) matrix.
#  k, can be viewed as extra columns.
# 
# Output: S, U, M, T
#  S, largest r invariant factors of A
#  UAV = S', AV = MS', where S is reverse form of S' i.e. S := sS'^-1.
#  T := UM lowerCmod S, with diagonal entries of 0's replaced by 1.
#       That is the strictly upper part of T corresponds to U . M.
#       The diagonal entry are all 1.
#       The stirctly lower triangular part are U . M cmod S and they are all 0's.
#
SNF := proc(A, n, r, k, s)
  global F;
  global timesBase, timeBase, timesExtractMatrixGCDHelper, timesVectorGCD, timeExtract, timeVector;
  while(true) do
    F := Matrix(n, k+1);
    dt := getDT(s);
    AA := ifelse(dt <> float[8], copy(A[1..-1, 1..r]), Mod(s, A[1..-1, 1..r], dt));
    E := ifelse(dt <> float[8], copy(A[1..-1, r+1..-1]), Mod(s, A[1..-1, r+1..-1], dt));
    #E := copy(A[1..-1, r+1..-1]);
    U := Matrix(min(r, n), n);
    U := ifelse(dt <> float[8], U, Mod(s, U, dt));
    M := Matrix(n, min(r, n));
    M := ifelse(dt <> float[8], M, Mod(s, M, dt));
    T := Matrix(min(r, n), min(r, n), (i, j) -> if (i = j) then 1 else 0 end if);
    T := ifelse(dt <> float[8], T, Mod(s, T, dt));
    S := Vector(min(r, n), 1);
    timeBase := 0;
    timeExtract := 0;
    timeVector := 0;
    if (dt = float[8]) then
      res := SNFHelperFloat(AA, E, U, M, S, T, n, r, k, s, 1, r, s);
      AA := map(round, AA);
      U := map(round, U);
      M := map(round, M);
      T := map(round, T);
      T := ReplaceDiagonal(T, r, 0, 1);
    else
      res := SNFHelper(AA, E, U, M, S, T, n, r, k, s, 1, r, s);
    end if;
    timesBase := [op(timesBase), timeBase];
    timesExtractMatrixGCDHelper := [op(timesExtractMatrixGCDHelper), timeExtract];
    timesVectorGCD := [op(timesVectorGCD), timeVector];

    #if (res <> false) then
    #  T[1, 1..-1] := LinearAlgebra:-Multiply(U[1, 1..-1], M[1..-1, 1..-1]);
    #  T[1,1] := 1;
    #end if;
    # TODO: Delete this after.
    #T := U . M;
    #for i to min(r, n) do for j from 1 to i do T[i,j] := 0 end do end do;
    #T := ReplaceDiagonal(T, min(r, n), 0, 1);

    # TODO: Handle false case.
    if (res = false) then return false end if;
    return S, U, M, T;
    #if (res) then
    #  return S, U, M;
    #end if;
  end do;
end proc;

# Return S, U, M s.t. UAV = S', AV = MS', where S is reverse form of S'. So S := sS'^-1.
SNFHelperFloat := proc(A, E, U, M, S, T, n, r, k, s, lower, upper)
  global F;
  global timeBase, timesExtractMatrixGCDHelper, timesVectorGCD, timeExtract, timeVector;
  if (s = 1) then return true end if;
  if (lower = upper) then
    start := time[real]();
    F[1..-1, 1] := map(round, A[1..-1, lower]);
    F[1..-1, 2..-1] := map(round, E);
    _, p := extractMatrixGCDHelper(s, F, n, k+1);
    u := vectorGCD(s, p, n);
    startExtract := time[real]();
    c, b, g := reScale(u . p, s);
    e := modp(b + c * s / g, s);
    p := Mod(s, p, float[8]);
    u := Mod(s, u, float[8]);
    Multiply(s, e, p, p);
    up := convert(Multiply(s, u, p), 'integer');
    if (up = 0) then return true end if;
    m := Mod(s/up, p/up, float[8]);
    j := -lower;
    u := Mod(s/up, u, float[8]);
    U[j, 1..-1] := u;
    M[1..-1, j] := m;
    S[j] := s/up;
    muE := Multiply(s, M[1..-1, j], Multiply(s, U[j,1..-1], E));
    #muE := MultiplyHelper(s, m, LeftSparseMult(u, E, 1, n, k, s));
    AddMultiple(s, s - 1, E, muE, E);
    E[1..-1, 1..-1] := E / up;
    if (j <> -1) then
      startVectorGCD := time[real]();
      T[j, j+1..-1] := Multiply(S[-1], U, j, 1..-1, M, 1..-1, j+1..-1);
      endVectorGCD := time[real]();
      timeVector := timeVector + endVectorGCD - startVectorGCD;
    end if;

    endTime := time[real]();
    timeBase := timeBase + endTime - start;
    #if (not type(E, 'Matrix(integer)')) then print("something is wrong with up"); return false end if;
    return true;
  end if;

  mid := floor((lower + upper)/2);
  res := SNFHelperFloat(A, E, U, M, S, T, n, r, k, s, lower, mid);

  if (not res) then return false end if;
  Q := LinearAlgebra:-Modular:-Inverse(s, T[-mid..-lower, -mid..-lower]);
  MQUA := Multiply(s, M, 1..-1, -mid..-lower, Multiply(s, Q, Multiply(s, U, -mid..-lower, 1..-1, A, 1..-1, mid+1..upper)));
  AddMultiple(s, s-1, A, 1..-1, mid+1..upper, MQUA, A, 1..-1, mid+1..upper);
  A[1..-1, mid+1..upper] := A[1..-1, mid+1..upper] * (S[-mid] / s);

  return SNFHelperFloat(A, E, U, M, S, T, n, r, k, S[-mid], mid+1, upper);
end proc;



# Return S, U, M s.t. UAV = S', AV = MS', where S is reverse form of S'. So S := sS'^-1.
SNFHelper := proc(A, E, U, M, S, T, n, r, k, s, lower, upper)
  global F;
  global timeBase, timesExtractMatrixGCDHelper, timesVectorGCD, timeExtract, timeVector;
  if (s = 1) then return true end if;
  if (lower = upper) then
    start := time[real]();
    F[1..-1, 1] := A[1..-1, lower];
    F[1..-1, 2..-1] := E;
    startExtract := time[real]();
    _, p := extractMatrixGCDHelper(s, F, n, k+1);
    endExtract := time[real]();
    startVectorGCD := time[real]();
    u := vectorGCD(s, p, n);
    endVectorGCD := time[real]();
    timeExtract := timeExtract + endExtract - startExtract;
    timeVector := timeVector + endVectorGCD - startVectorGCD;
    #u := vectorGCD(s, -p, n);
    #u := Multiply(s, s-1, u);
    c, b, g := reScale(u . p, s);
    e := modp(b + c * s / g, s);
    Multiply(s, e, p, p);
    up := Multiply(s, u, p);
    if (up = 0) then return true end if;
    m := modp(p/up, s/up);
    j := -lower;
    u := modp(u, s/up);
    U[j, 1..-1] := u;
    M[1..-1, j] := m;
    S[j] := s/up;
    #muE := MultiplyHelper(s, m, MultiplyHelper(s, u, E)); 
    muE := MultiplyHelper(s, m, LeftSparseMult(u, E, 1, n, k, s));
    AddMultiple(s, s - 1, E, muE, E);
    E[1..-1, 1..-1] := E / up;
    if (j <> -1) then
      #T[j, j+1..-1] := U[j, 1..-1] . M[1..-1, j+1..-1];
      T[j, j+1..-1] := LeftSparseMult(U[j, 1..-1], M[1..-1, j+1..-1], 1, n, lower - 1);
    end if;

    endTime := time[real]();
    timeBase := timeBase + endTime - start;
    if (not type(E, 'Matrix(integer)')) then print("something is wrong with up"); return false end if;
    return true;
  end if;

  mid := floor((lower + upper)/2);
  res := SNFHelper(A, E, U, M, S, T, n, r, k, s, lower, mid);

  # Need to claculate product of (I - miui)
  #for i from lower to mid do
  #  muA := MultiplyHelper(s, M[1..-1, -i], MultiplyHelper(s, U[-i, 1..-1], A[1..-1, mid+1..upper])); 
  #  AddMultiple(s, s - 1, A, 1..-1, mid+1..upper, muA, A, 1..-1, mid+1..upper);
  #end do;

  #T[-mid..-lower, -mid..-1] := imlMult(U[-mid..-lower, 1..-1], M[1..-1, -mid..-1]);
  #T[-mid..-lower, -mid..-lower] := ReplaceDiagonal(lowerCmod(T[-mid..-lower, -mid..-lower], S[-mid..-lower]), mid-lower+1, 0, 1);

  if (not res) then return false end if;
  Q := LinearAlgebra:-Modular:-Inverse(s, T[-mid..-lower, -mid..-lower]);
  #Q := U[-mid..-lower, 1..-1] . M[1..-1, -mid..-lower];
  #for i to mid-lower+1 do for j from 1 to i do Q[i,j] := 0 end do end do;
  #ReplaceDiagonal(Q, mid-lower+1, 0, 1);
  #MQUA := MultiplyHelper(s, M[1..-1, -mid..-lower], MultiplyHelper(s, LinearAlgebra:-Modular:-Inverse(s, Q), MultiplyHelper(s, U[-mid..-lower, 1..-1], A[1..-1, mid+1..upper])));
  MQUA := MultiplyHelper(s, M[1..-1, -mid..-lower], MultiplyHelper(s, Q, LeftSparseMult(U[-mid..-lower, 1..-1], A[1..-1, mid+1..upper], mid-lower+1, n, upper-mid)));
  AddMultiple(s, s-1, A, 1..-1, mid+1..upper, MQUA, A, 1..-1, mid+1..upper);

  A[1..-1, mid+1..upper] := A[1..-1, mid+1..upper] * (S[-mid] / s);
  if (not type(A, 'Matrix(integer)')) then print("**something is wrong"); return false end if;

  return SNFHelper(A, E, U, M, S, T, n, r, k, S[-mid], mid+1, upper);
end proc;

# computeProjBasis(P, n, r, s)
#
# Input:
#   P, 2n x r matrix.
#   n, r, s, positive integer
#
# Pre: The last n rows of P are all 0.
#
# Output: S, U, M
#   S, r x r integer matrix.
#   U, r x n Mod s matrix.
#   M, n x r Mod s matrix.
#   Let P := [P1], Where P1 and P2 are n x r matrices.
#            [P2]
#   -UP1V = F mod s and P1V = MF mod s, where F is in reverse smith form of P1
#   over Mod s and S := sF^-1.
#
#   Calculates the projection basis for Proj(sI, P).
#   Proj(sI, P) := {v length 2n vector | vP/s is integral vector}. And the
#   basis is:
#   [I           ]   [I -M ]    [I    ]
#   [ I          ] * [  I  ]  * [  I  ]
#   [  sF^{-1}   ]   [   I ]    [-U I ]
#   [         I_m]   [    I]    [    I]
#
computeProjBasis := proc(P, n, r, s, k);
  global timesSNFModS;
  local P1, U, S, V, U1, M, F, UNew, VNew, v, p, u, m;

  snew := s;
  U := Matrix(min(r, n), n);
  M := Matrix(n, min(r, n));
  S := Vector(min(r, n), 1);
  P1 := SubMatrix(P, 1..n, 1..r+k);
  m := Vector(n);
  colNum := min(5, r);
  start := time[real]();
  ret := SNF(P1, n, r, k, s);
  if (ret = false) then return false end if;
  S, U, M, T := ret;
  endtime := time[real]();
  timesSNFModS := [op(timesSNFModS), endtime - start];
  return S, U, -M, T;
  
  # We should get back a S (n x r) in normal smith form, U (r x n), M (n x r).
  # s.t. UP1V = S mod s and P1V = MS mod s.
  for j from min(r, n) to 1 by -1 do
    v := Random(2, r, colNum, integer);
    randP := MultiplyHelper(snew, P1, v);
    #randP := P1[1..-1, j..min(j+colNum-1, r)];
    #_, p := extractMatrixGCDHelper(snew, randP, n, min(j+colNum-1, r) - j + 1);
    _, p := extractMatrixGCDHelper(snew, randP, n, colNum);
    u := vectorGCD(snew, -p, n);
    u := Multiply(snew, snew-1, u);
    up := u . p;
    c, b, g := reScale(u . p, snew);
    e := modp(b + c * snew / g, snew);
    Multiply(snew, e, p, p);
    up := modp(e * up, snew);
    if (up = 0) then break end if;
    if (not type(P1/up, 'Matrix(integer)')) then ++j; next; end if;
    #if (not type(P1/up, 'Matrix(integer)')) then break; end if;
    S[j] := up;
    m := modp(p/up, snew/up);
    U[j, 1..-1] := u;
    M[1..-1, j] := m;
    # Only update the necessary parts of P1.
    muP1 := MultiplyHelper(snew, m, MultiplyHelper(snew, u, P1)); 
    AddMultiple(snew, snew - 1, P1, muP1, P1);
    P1 := P1/up;
    snew := snew/up;
    S[j] := snew;
  end do;

  return S, U, -M;

  ##S, U, V := SmithFormModS(P1, n, r, s);
  #M := reconstructM(s, S, P1, V, n, r);
  ## Convert M into reverse smith form.
  #reverseCol(s, M, min(n, r));
  #F, UNew, VNew := reverseSmithForm(s, S, U, V, n, r);

  ## If the entry for smith form vanish then replace with `s`.
  #return map(a -> if (a <> 0) then s/a else 1 end if, F), UNew, -M;
end proc;

#computeProjBasis := proc(P, n, r, s, k);
#  global timesSNFModS;
#  local P1, U, S, V, U1, M, F, UNew, VNew, v, p, u, m;
#
#  snew := s;
#  U := Matrix(min(r, n), n);
#  M := Matrix(n, min(r, n));
#  V := Matrix(r+k, r);
#  S := Vector(min(r, n), 1);
#  #P1 := SubMatrix(P, 1..n, 1..r+k);
#  #P1 := P;
#  P1 := Mod(s, P, integer);
#  m := Vector(n);
#  #start := time[real]();
#  #ret := SNF(P1, n, r, k, s);
#  #if (ret = false) then return false end if;
#  #S, U, M, T := ret;
#  #endtime := time[real]();
#  #timesSNFModS := [op(timesSNFModS), endtime - start];
#  #return S, U, -M, T;
#  
#  # We should get back a S (n x r) in normal smith form, U (r x n), M (n x r).
#  # s.t. UP1V = S mod s and P1V = MS mod s.
#  colNum := min(5, r);
#  colNum := r+5;
#  for j from min(r, n) to 1 by -1 do
#    #v := Random(2, r, colNum, integer);
#    #randP := MultiplyHelper(snew, P1, v);
#    #randP := P1[1..-1, j..min(j+colNum-1, r)];
#    #_, p := extractMatrixGCDHelper(snew, randP, n, min(j+colNum-1, r) - j + 1);
#    v, p := extractMatrixGCDHelper(snew, P1[1..n, -j-5..-j], n, 5);
#    u := vectorGCD(snew, p, n);
#    #u := Multiply(snew, snew-1, u);
#    up := u . p;
#    c, b, g := reScale(up, snew);
#    e := modp(b + c * snew / g, snew);
#    Multiply(snew, e, p, p);
#    up := modp(u . p, snew);
#    if (up = 0) then break; end if;
#    if (not type(P1/up, 'Matrix(integer)')) then return false end if;
#    S[j] := up;
#    m := modp(p/up, snew/up);
#    u := modp(u, snew/up);
#    print("u . m: ", modp(u . m, s));
#    U[j, 1..-1] := u;
#    M[1..-1, j] := m;
#    V[1..-1, j] := modp(e*v, snew/up);
#    # Only update the necessary parts of P1.
#    print(modp(LinearAlgebra:-Determinant((LinearAlgebra:-IdentityMatrix(n) - m . u)), s/up));
#    muP1 := MultiplyHelper(snew, m, MultiplyHelper(snew, u, P1)); 
#    AddMultiple(snew, snew - 1, P1, muP1, P1);
#    P1 := P1/up;
#    S[j] := snew/up;
#    snew := snew/up;
#  end do;
#
#  T := U . M;
#  for i from 1 to r do for j from 1 to i do T[i, j] := 0 end do end do;
#  for i from 1 to r do T[i,i] := 1 end do;
#
#  Q := LinearAlgebra:-Modular:-Inverse(s, T);
#  interface(rtablesize=100);
#  print(S);
#  print(P);
#  print(U[-1]);
#  print(U . P);
#  print(modp(U . P, s));
#  print(P . V);
#  print(P . V[1..-1, r]);
#  print(modp(P . V, s));
#  print(modp(U . P . V, s));
#  print(modp(Q . U . P, s));
#  print(modp(Q . U . P . V, s));
#  print(modp(P - M . Q . U . P, s));
#  print(V);
#  #print(LinearAlgebra:-Determinant(V));
#
#
#  return S, U, -M, T;
#
#  ##S, U, V := SmithFormModS(P1, n, r, s);
#  #M := reconstructM(s, S, P1, V, n, r);
#  ## Convert M into reverse smith form.
#  #reverseCol(s, M, min(n, r));
#  #F, UNew, VNew := reverseSmithForm(s, S, U, V, n, r);
#
#  ## If the entry for smith form vanish then replace with `s`.
#  #return map(a -> if (a <> 0) then s/a else 1 end if, F), UNew, -M;
#end proc;
