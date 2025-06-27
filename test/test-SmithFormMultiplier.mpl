read "SmithMultiplier/SmithFormMultiplier.mpl";
read "test/test-Helpers.mpl";

macro(SmithForm=LinearAlgebra:-SmithForm);

kernelopts(assertlevel=1);

unfoundAns := Array([]);
unfoundCount := 0;
tot := 0;

checkValidMultiplier := proc(A, n, S, U, V)
  ASSERT(abs(Determinant(U)) = 1);
  ASSERT(abs(Determinant(V)) = 1);
  SCorrect := SmithForm(A);
  ASSERT(Equal(SCorrect, S));
  ASSER(Equal(A . V, U . S));
end proc;

checkValid := proc(A, n, S := false, U := false, V := false)
  global tot;
  ++tot;
  if (S <> false) then
    SS := DiagonalMatrix(S);
    checkValidMultiplier(A, n, SS, U, V);
  else
    recordMatrices(A, n);
  end if;
end proc;

recordMatrices := proc(A, n)
  global unfoundCount, unfoundAns;
  Append(unfoundAns, A, inplace = true);
  ++unfoundCount;
end proc;

checkRandom := proc()
  local A, n, roll, ret;
  roll := rand(2..30);
  for i from 1 to 100 do
    n := roll();
    A := randomNonsingular(n);
    ret := SmithFormMultipliers(A, n);
    checkValid(A, n, ret);
  end do;
end proc;

check := proc()
  checkRandom();
  print("unfoundAns:", unfoundAns);
  print("total test cases:", tot);
  print("total unfound:", unfoundCount);
  print("Success %:", 1 - unfoundCount/tot);
end proc;

check();
