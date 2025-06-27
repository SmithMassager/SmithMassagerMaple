read "config.mpl";

macro(Dimension=LinearAlgebra:-Dimensions);

#### These are naive implementation of the c external library calls ####
# Solve for Ax = b
naiveIMLSolveHelper := proc(A, b)
  local res, d, n, m, i, j;

  res := LinearAlgebra:-LinearSolve(A, b);
  n, m := Dimensions(res);
  d := 1;
  for i from 1 to n do
    for j from 1 to m do
      d := lcm(d, denom(res[i,j]));
    end do;
  end do;

  return d * res, d;
end proc;


# Calculate A(yb) mods m
naiveIMLMult := proc(A, b, y := -1, m := -1)
  local res;
  res := A . b;
  if (y <> -1) then res := y * res end if;
  if (m <> -1) then res := mods(res, m) end if;
  return res;
end proc;

# Macros for CLIB and naive implementation in maple.
if (CLIB) then
  read "lib/extern.mpl";
  macro(MultiplyHelper=fastMultiplyHelper);
  macro(IntCertificate=fastIntCertC);
  macro(UniCert=fastUniCert);
else
  macro(MultiplyHelper=naiveMultiplyHelper);
  macro(IntCertificate=naiveIntCertificate);
  macro(UniCert=naiveUniCert);
  imlSolveHelper := naiveIMLSolveHelper;
  imlMult := naiveIMLMult;
  fmpzMult := naiveIMLMult;
end if;
