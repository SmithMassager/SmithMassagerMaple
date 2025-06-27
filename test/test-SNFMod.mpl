read "SmithMassager/SNFMod.mpl";
macro(Mod=LinearAlgebra:-Modular:-Mod);
macro(Equal=LinearAlgebra:-Equal);
macro(SmithForm=LinearAlgebra:-SmithForm);
macro(MatrixInverse=LinearAlgebra:-MatrixInverse);
macro(RandomMatrix=LinearAlgebra:-RandomMatrix);
macro(Mod=LinearAlgebra:-Modular:-Mod);
macro(Determinant=LinearAlgebra:-Modular:-Determinant);
macro(DiagonalMatrix=LinearAlgebra:-DiagonalMatrix);
macro(HermiteForm=LinearAlgebra:-HermiteForm);
kernelopts(assertlevel=1);
randomize(3);

smallestPivot := proc(A, i)
  description "Find the smallest entry in column i of A that is nonzero";
  local n, m, j, res;
  n, m := Dimension(A);
  res := 1;
  for j from 2 to n do
    res := ifelse(A[res,i] <> 0 or (A[j,i] < A[res,i] and A[j, i] <> 0), j, res);
  end do;
  return res;
end proc;

# Returns gcd of v
vecGcd := proc(v)
  description "Given a vector v return gcd(v).";
  local rowV, colV, g, i, vNew;
  
  rowV := Dimension(v);
  vNew := v;

  g := vNew[1];
  for i from 2 to rowV do
    g := gcd(g, vNew[i]);
  end do;
  return g;
end proc;

# Returns gcd of A
matrixGcd := proc(A)
  description "Given a vector A return gcd(A).";
  local rowA, colA, g, i, rowG;
  rowA, colA := Dimension(A);
  rowG := 1..rowA;
  g := vecGcd(A[rowG, 1]);
  for i from 2 to colA do
    g := gcd(g, vecGcd(A[rowG, i]));
  end do;
  return g;
end proc;

stabVec := proc(v, s)
  description "Given a column vector v, return a row vector c s.t. c v = \\
  gcd(v, s), with c_1 = 1";
  local g, n, vNew, c, i;
  vNew := convert(v, Vector);
  n := Dimension(vNew);
  g := vNew[1];
  c := Vector[row](n);
  c[1] := 1;
  for i from 2 to n do
    c[i] := Stab(vNew[i], g, s);
  end do;
  return c;
end proc;

#------------------------ Begin Eliminate Col ---------------------------------
checkEliminateCol := proc()
  local A, U, s, Acopy, Vcopy;
  s := 6;
  A := Matrix([
    [2,2, 2],
    [8,4, 4],
    [2,0, 0]
  ]);
  U := Matrix([
   [1,1,0],
   [0,1,0],
   [0,0,1]
  ]);
  Acopy := Matrix(A);
  Ucopy := Matrix(U);
  eliminateCol(A, U, 3, 3, s, 1);
  ASSERT(Equal(A, Matrix([[2,2,2],[2,2,2],[5,4,4]])));
  ASSERT(Equal(U, Ucopy));
end proc;
#------------------------ End Eliminate Col -----------------------------------


#------------------------ Begin Extract GCD -----------------------------------
checkextractGCD := proc()
  local gen, i, s, v1, v2, c;
  gen := rand(2..100);
  for i from 0 to 1000 do
    s := gen();
    v1 := RandomMatrix(2,1, generator = 0 ..(s-1));
    v2 := RandomMatrix(2,1, generator = 0 ..(s-1));
    c := extractGCD(v1[1..2,1], v2[1..2,1], s);
    if (gcd(vecGcd(v1 + c * v2),s) <> gcd(gcd(vecGcd(v1), vecGcd(v2)),s)) then
      print("FOUND AN COUNTER EXAMPLE");
      print(vecGcd(v1));
      print(vecGcd(v2));
      print(gcd(vecGcd(v1), vecGcd(v2)));
      print(vecGcd(v1 + c * v2));
      print(s);
      print(c);
      print(v1);
      print(v2);
      ASSERT(true);
    fi;
  end do;
end proc;

checkExtractMatrixGCD := proc()
  local A, U, s, Acopy, Vcopy;
  s := 8;
  A := Matrix([
    [4, 2, 3],
    [8, 4, 2],
    [4, 0, 0]
  ]);
  U := Matrix([
   [1,0,0],
   [0,1,0],
   [0,0,1]
  ]);
  Acopy := Matrix(A);
  Ucopy := Matrix(U);
  extractMatrixGCD(s, A, U, 3, 3, 1);
  ASSERT(vecGcd(A[1..-1, 1]) = 1);

  A := Matrix([
    [4, 2, 3],
    [8, 4, 2],
    [4, 0, 0]
  ]);
  U := Matrix([
   [1,0,0],
   [0,1,0],
   [0,0,1]
  ]);
  extractMatrixGCD(s, A, U, 3, 3, 2);
  print(A);
  ASSERT(igcd(vecGcd(A[2..-1, 2]), s) = 2);

  A := Matrix([
    [4, 2, 3],
    [8, 3, 2],
    [4, 0, 0]
  ]);
  U := Matrix([
   [1,0,0],
   [0,1,0],
   [0,0,1]
  ]);
  extractMatrixGCD(s, A, U, 2, 3, 3, row);
  ASSERT(igcd(vecGcd(A[1, 1..-1]), s) = 1);
end proc;
#------------------------ End Extract GCD -------------------------------------

#------------------------ Begin Check Rescale ---------------------------------
checkRescale := proc()
  local a, s, c;
  a := 48;
  s := 60;
  c, u, g := reScale(a, s);
  ASSERT(igcd((u + c*s/g), s) = 1);
  ASSERT(igcd(g, s) = modp(a * (u + c*s/g), s));
end proc;
#-------------------------- End Check Rescale ---------------------------------


#------------------------ Begin Check Stab ------------------------------------
checkStab := proc()
  local gen, i, a, b, s, c, c1;
  gen := rand(2..100);
  for i from 0 to 1000 do
    s := gen();
    a := gen() mod s;
    b := gen() mod s;
    c := Stab(a, b, s);
    if (igcd(modp(a + c * b, s), s) <> igcd(a, igcd(b, s))) then
      printf("Found a counter example a, b, s, c: %d %d %d %d \n", a, b, s, c);
    fi;
  end do;
end proc;
#-------------------------- End Check Stab ------------------------------------

#------------------------ Begin Check SmithFormModS ---------------------------
# Check that UAV = S is a valid smith normal form decomposition under Z(s).
verifySmithFormModS := proc(A, U, V, S, s)
  local Sz, Uz, Vz, i, n, m, res, j;
  res := true;
  n, m := Dimension(A);
  r := min(n, m);
  Sz, Uz, Vz := SmithForm(A, output=[':-S', ':-U', ':-V']);
  ASSERT(Equal(HermiteForm(<MTM:-transpose(U), s * Identity(2,r,integer)>),
              <Identity(2,r,integer), Matrix(n,r)>));
  ASSERT(gcd(Determinant(s, V), s) = 1, "V not unimodular.");
  ASSERT(Equal(Mod(s, U . A . V, integer), DiagonalMatrix(S, min(n,m), m)), "UAV != S.");
  for i from 1 to min(n, m) do
    ASSERT(gcd(Sz[i,i], s) mod s = S[i], "Wrong invariant factor.");
    ASSERT(gcd(S[i], s) mod s = S[i], "Unnormalized entries.");
  end do;
end proc;

# Generate a random unimodular (n x n) matrix under Z(s).
randUniMatrix := proc(s, n)
  local m;
  while(true) do
    m := RandomMatrix(n, n, generator = 0 ..(s-1));
    if (gcd(Determinant(s, m), s) = 1) then
      return m;
    end if;
  end do;
end proc;

checkSmithNormalFormModSRec := proc()
  local s, S, U, V, A, i;
  s := 6;
  P1 := Matrix([[0], [1]]);
  S, U, V := SmithFormModS(P1, 2, 1, s);
  verifySmithFormModS(P1, U, V, S, s);
  s := 2^5 * 3^4 * 5^6 * 11^4;
  A := Matrix([
    [2, 0, 0,  0],
    [0, 6, 0,  0],
    [0, 0, 6*13, 0],
    [0, 0, 0,  0],
    [0, 0, 0,  0],
    [0, 0, 0,  0]
  ]);
  for i from 1 to 50 do
    A := Multiply(s, A, randUniMatrix(s, 4));
    A := Multiply(s, randUniMatrix(s, 6), A);
  end do;
  S, U, V := SmithFormModS(A, 6, 4, s);
  verifySmithFormModS(A, U, V, S, s);

  s := 3^4 * 5^6 * 11^4;
  A := Matrix([
    [2, 0, 0,  0],
    [0, 6, 0,  0]
  ]);
  for i from 1 to 30 do
    A := Multiply(s, A, randUniMatrix(s, 4));
    A := Multiply(s, randUniMatrix(s, 2), A);
  end do;
  S, U, V := SmithFormModS(A, 2, 4, s);
  verifySmithFormModS(A, U, V, S, s);
end proc;

checkSmithNormalFormModSSqr := proc()
  local s, S, U, V, A, i;
  s := 3^4 * 5^6 * 11^4;
  A := Matrix([
    [2, 0, 0,  0],
    [0, 6, 0,  0],
    [0, 0, 6 * 11,  0],
    [0, 0, 0, 6 * 11 * 33]
  ]);
  for i from 1 to 30 do
    A := Multiply(s, A, randUniMatrix(s, 4));
    A := Multiply(s, randUniMatrix(s, 4), A);
  end do;
  S, U, V := SmithFormModS(A, 4, 4, s);
  verifySmithFormModS(A, U, V, S, s);

  s := 3^4 * 5^6 * 11^4;
  A := Matrix([
    [2, 0, 0,  0],
    [0, 6, 0,  0],
    [0, 0, 6 * 11,  0],
    [0, 0, 0, 6 * s]
  ]);
  for i from 1 to 30 do
    A := Multiply(s, A, randUniMatrix(s, 4));
    A := Multiply(s, randUniMatrix(s, 4), A);
  end do;
  S, U, V := SmithFormModS(A, 4, 4, s);
  verifySmithFormModS(A, U, V, S, s);

end proc;

checkSmithNormalFormModSRnd := proc()
  local s, S, U, V, A, gen, n, m, dimGen, i;
  gen := rand(2..10000);
  dimGen := rand(2..100);
  s := gen();
  n := dimGen();
  m := dimGen();
  A := RandomMatrix(n, m, generator = 0 ..(s-1));
  for i from 1 to 30 do
    A := Multiply(s, A, randUniMatrix(s, m));
    A := Multiply(s, randUniMatrix(s, n), A);
  end do;

  S, U, V := SmithFormModS(A, n, m, s);
  verifySmithFormModS(A, U, V, S, s);
end proc;

checkSmithNormalFormModLarge := proc()

  A := Matrix([[0,    3992435006,             0,    0,             0,    3992435006,    0],
               [0,    3992435006,             0,    0,             0,    3992435006,    0],
               [0,             0,    5988652509,    0,    5988652509,    3992435006,    0],
               [0,    3992435006,             0,    0,             0,    3992435006,    0]]);
  s := 7984870012;
  S, U, V := SmithFormModS(A, 4, 7, s);
  verifySmithFormModS(A, U, V, S, s);

end proc;

checkSmithNormalFormModS := proc()
  checkSmithNormalFormModSSqr();
  checkSmithNormalFormModSRec();
  checkSmithNormalFormModLarge();
  checkSmithNormalFormModSRnd();
end proc;
#-------------------------- End Check SmithFormModS ---------------------------

checkExtractGCD();
checkRescale();
checkStab();
checkSmithNormalFormModS();
checkEliminateCol();
checkExtractMatrixGCD();
