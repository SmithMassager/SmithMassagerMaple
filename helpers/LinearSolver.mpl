macro(Create=LinearAlgebra:-Modular:-Create);
macro(Dimensions=LinearAlgebra:-Dimensions);
macro(Equal=LinearAlgebra:-Equal);
macro(ModLinearSolve=LinearAlgebra:-Modular:-LinearSolve);
macro(Mod=LinearAlgebra:-Modular:-Mod);
macro(AddMultiple=LinearAlgebra:-Modular:-AddMultiple);

read "helpers/basic.mpl";

# LinearSolve(A, B)
#
# Pre: Assumes AX = B have integer solution X.
#
# Input: 
#   A, B matrix of suitable size
#   thres, a threshold on the answer A^(-1) B
#
# Output: X, where AX = B, or false if no solution is found.
#
ChineseLinearSolve := proc(A, B, thres, s)
  local primes, p, X, ABp, n1, m1, n2, m2, l, prod, Ans, oldAns;
  n1, m1 := Dimensions(A);
  n2, m2 := Dimensions(B);

  Assert(n1=n2);
  Assert(n1=m1);

  l := 0;
  p := 2292392;
  prod := 1;
  primes := [];
  X := [];

  print("thres", thres);

  while(prod < thres) do
    p := prevprime(p);
    if (igcd(s, p) > 1) then next end if;
    ++l;
    prod *= p;
    print("prod", prod);
    primes := [op(primes), p];
    start := time[real]();
    ABp := Create(p, n1, m1+m2, float[8]);
    ABp[1..n1, 1..m1] := Mod(p, A, float[8]);
    ABp[1..n1, -m2..-1] := Mod(p, B, float[8]);
    endTime := time[real]();
    start := time[real]();
    ModLinearSolve(p, ABp, m2);
    endTime := time[real]();
    print("time for mod solve", endTime - start);
    X := [op(X), map(round, ABp[1..m1, -m2..-1])];
  end do;

  start := time[real]();
  if (l = 1) then
    return mods(X[1], primes[1]);
  end if;

  Ans := ChineseRemainder(primes, X, l, m1, m2, prod);
  endTime := time[real]();
  print("time for chinese remainder", endTime - start);

  return Ans;
end proc;

# ChineseRemainder(p, X, k, n, m, prod)
#
# Input: 
#   X, list of k matrices of size n x m, where X_i is mod p_i.
#   p, list of k coprime numbers.
#   prod, product of p.
#
# Output: Y where Y <-> (X_1 mod p_1, ..., X_n mod p_n)
#
ChineseRemainder := proc(p, X, k, n, m, prod)
  local Y;

  Y := Matrix(n, m);
  print(n, m);
  for i from 1 to k do
    q := iquo(prod, p[i]);
    L := modp(modp(1/q, p[i]) * q, prod);
    Xi := X[i];
    Y := add_maple(Y, Xi, L, prod);
  end do;

  return Y;
end proc;
