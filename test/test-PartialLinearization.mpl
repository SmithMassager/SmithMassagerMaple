read "test/test-Helpers.mpl";
read "helpers/partialLinearization.mpl";
with(LinearAlgebra);

checkBasic := proc()
  local A, n, m, Abar, P, Q, d, e, ebar;
  A := Matrix([[100]]);
  Abar, n, m := ColPLHelper(A, 1, 1, Vector([1]), 3);
  ASSERT(Equal(Abar, Matrix([[4, 12], [-8, 1]])));

  A := Matrix([[2, 4, 44199, 3061969404], [4, 8, 19644, 765492351], [7, 8, 44199, 5358446457], [7, 5, 9822, 765492351]]);
  P, Q := getPerm(A, 4);
  d, e, ebar := getPermPLParam(A, 4);
  ASSERT(Equal(P, Vector([3, 1, 2, 4])));
  ASSERT(Equal(Q, Vector([4, 3, 2, 1])));
  ASSERT(d = 14);
  ASSERT(Equal(e, Vector([0, 0, 1, 2])));
  ASSERT(Equal(ebar, Vector([1, 0, 2, 0, 0, 0, 0])));
end proc;

checkPartial := proc(A, Abar)
  local n, m, nbar, mbar;
  n, m := Dimension(A);
  nbar, mbar := Dimension(Abar);
  ASSERT(Determinant(A) = Determinant(Abar));
  ASSERT(Equal(MatrixInverse(A), MatrixInverse(Abar)[1..n, 1..m]));
end proc;

check := proc()
  local A, n, m, _, Abar, d, e;
  A := Matrix([[2, 4, 44199, 3061969404], [4, 8, 19644, 765492351], [7, 8, 44199, 5358446457], [7, 5, 9822, 765492351]]);
  n := 4; m := 4;
  d, e := getColPLParam(A, n, m);
  Abar, _, _ := ColPLHelper(A, n, m, e, d);
  print(Abar);
  checkPartial(A, Abar);
  Abar, _, _ := ColPL(-A, n, m);
  checkPartial(-A, Abar);
  Abar, _, _ := RowPL(A, n, m);
  checkPartial(A, Abar);
  Abar, _, _ := RowPL(-A, n, m);
  checkPartial(-A, Abar);
  Abar, _, _ := PermPL(A, n, m);
  checkPartial(A, Abar);
  Abar, _, _ := PermPL(-A, n, m);
  checkPartial(-A, Abar);
end proc;

checkAll := proc()
  check();
  checkBasic();
end proc;

checkAll();
