read "macros.mpl";

macro(RandomVector=LinearAlgebra:-RandomVector);
macro(RandomMatrix=LinearAlgebra:-RandomMatrix);
macro(IntegerLinearSolve=LinearAlgebra:-Modular:-IntegerLinearSolve);

# ExtractDen(V, n)
#
# Input:
#   V, a vector of size n
#
# Output: Res
#   Res, a vector containing denominator of V.
#
ExtractDen := proc(V, n)
  local res, i;
  res := Vector(n);
  for i from 1 to n do
    res[i] := denom(V[i]);
  end do;
  return res;
end proc;

# VecLCM(V)
#
# Input:
#   V, a vector of size n
#
# Output, LCM(V)
#
VecLCM := proc(V, n)
  local res, i;
  res := 1;
  for i from 1 to n do
    res := lcm(res, V[i]);
  end do;
  return res;
end proc;

# LargestInvariantFactor(A, n)
#
# Input:
#   A, nxn matrix in Z.
#   n, positive integer.
#
# Output: s
#   s, With high probability s is the largest invariant factor of A.
#
LargestInvariantFactor := proc(A, n, startDim := 5)
  local mx, M, s, b, x, xDen, i;
  mx := max(abs(A));
  # Change this to c style.
  M := 10 * ceil(6 + 2 * n * (log2(n) + log2(mx)));
  s := 1;

  B := RandomMatrix(n, startDim, generator = 0..M-1);
  #X, d := dixonSolveHelper(A, B);
  X, d := imlSolveHelper(A, B);
  Y := X/d;
  for i from 1 to 5 do
    xDen := ExtractDen(Y[1..-1, i], n);
    s := lcm(s, VecLCM(xDen, n));
  end do;

  Y := s * Y;
  if (not type(Y, 'Matrix'(integer))) then
    print("Failed to find the largest invariant factor, retrying...");
    return LargestInvariantFactor(A, n, startDim);
  end if;

  return s, modp(Y, s);
end proc;
