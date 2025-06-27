read "helpers/doublePlusOneLift.mpl";
macro(Equal=LinearAlgebra:-Equal);
macro(Mod=LinearAlgebra:-Modular:-Mod);
macro(MatrixInverse=LinearAlgebra:-MatrixInverse);
kernelopts(assertlevel=1);

#------------------------ Begin Check --------------------------------
checkSpecialSolve := proc()
  local A, B, d, AInverseB;
  A := Matrix([[1,2], [3,5]]);
  B := Matrix([[1,2,3], [4,5,6]]);
  d := 6;
  AInverseB := SpecialSolve(A, B, d, 2, 3);
  ASSERT(Equal(AInverseB, MatrixInverse(A) . B));
end proc;

checkDoublePlusOneLift := proc()
  local A0, R, M, A;
  A := Matrix([777]);
  A0, R, M := DoublePlusOneLift(A, 1000, 1, 1);
  ASSERT(Equal(A0, Matrix([-287])));
  ASSERT(Equal(R[0], Matrix([223])));
  ASSERT(Equal(M[0], Matrix([-223])));
end proc;

checkComputeDoublePlusOneLift := proc()
  local A0, R, M, A, B, ans;
  A := Matrix([777]);
  A0, R, M := DoublePlusOneLift(A, 10000, 1, 1);
  B := Matrix([1]);
  ans := ComputeDoublePlusOneLift(1, 1, 10000, A0, R, M, B);
  ASSERT(ans[1,1] = mods(777^(-1), doublePlusOne(10000, 1)));
end proc;

checkIntCertificate := proc()
  local A, B, s;
  A := Matrix([[3, 0], [0, 21]]);
  B := Matrix([[1,2,3], [4,5,6]]);
  s := 21;
  ASSERT(Equal(IntCertificate(s, A, B, 2, 3), modp(s * A^(-1) . B, s)));

  s := 1;
  ASSERT(not IntCertificate(s, A, B, 2, 3));

  s := 105;
  A := Matrix([
         [-13, 10, -20, 27, 0, 0, 0, 0 ],
         [ 27, 30, 15,  30, 0, 0, 0, 0 ],
         [ 0,  15, 15,  6,  0, 0, 0, 0 ],
         [ -21, 0, -15, 9,  0, 0, 0, 0 ],
         [ 0,   0, 0,   0,  1, 0, 0, 0 ],
         [ 0,   0, 0,   0,  0, 1, 0, 0 ],
         [ 0,   0, 0,   0,  0, 0, 1, 0 ],
         [ 0,   0, 0,   0,  0, 0, 0, 1 ]
       ]);

  B := Matrix([[95], [41], [73], [37], [0], [0], [0], [0]]);
  ASSERT(Equal(IntCertificate(s, A, B, 8, 8), modp(s * A^(-1) . B, s)));
end proc;

checkUniCert := proc()
  local A;
  A := Matrix([[1,2], [2, 3]]);
  ASSERT(UniCert(A, 2));
  A := Matrix([[2,2],[2,2]]);
  ASSERT(not UniCert(A, 2));
  A := Matrix([[1,3],[1,5]]);
  ASSERT(not UniCert(A, 2));
  A := Matrix([
              [ 1,     0,     0,    0,     0,     0],
              [ 0,     9,     0,    0,     1,     8],
              [ 0,     0,    27,    0,    24,     4],
              [ 0,     0,     0,    1,     0,     0],
              [ 8,     8,     0,    0,     1,    -7],
              [26,    26,    26,    0,    26,    27]
             ]);
  ASSERT(UniCert(A, 6));
end proc;
#------------------------ End Check ----------------------------------

checkSpecialSolve();
checkDoublePlusOneLift();
checkComputeDoublePlusOneLift();
checkIntCertificate();
checkUniCert();
