macro(Dimension=LinearAlgebra:-Dimension);
macro(Reverse=ArrayTools:-Reverse);
macro(Swap=LinearAlgebra:-Modular:-Swap);
macro(Copy=LinearAlgebra:-Modular:-Copy);
macro(Create=LinearAlgebra:-Modular:-Create);
macro(Multiply=LinearAlgebra:-Modular:-Multiply);


# reverseCol(s, A, n)
#
# Input:
#   s,n positive integer.
#   A, n column Mod s matrix.
#
# Post: A modified with column reversed.
#
# Output: None.
#
reverseCol := proc(s, A, n):
  for i from 1 to floor(n/2) do 
    Swap(s, A, 1..-1, i, A, 1..-1, n-i+1);
  end do;
end proc;

# reverseRow(s, A, n)
#
# Input:
#   s,n positive integer.
#   A, n row Mod s matrix.
#
# Post: A modified with row reversed.
#
# Output: None.
#
reverseRow := proc(s, A, n):
  for i from 1 to floor(n/2) do 
    Swap(s, A, i, A, n-i+1);
  end do;
end proc;


# reverseSmithForm(s, S, U, V, n, m)
#
# Input: 
#   n, m, s, Positive integer.
#   S, min(n, m) vector Smith normal form with entry in Z(s).
#   U, the top r x n rows of a unimodular Mod s matrix, say W.
#      Where r := min(n,m)
#   V, an m x m Mod s unimodular matrix
#   Note: W A V = Diag(S[1],S[2],...,S[r]) is in Smith form mod s 
#
# Output:
#   S, U, V but in reverse smith form.
#
reverseSmithForm := proc(s, S, U, V, n, m)
  local SNew, VNew, UNew, i, r;
  r := min(n,m);
  SNew := Vector(S);
  for i from 1 to floor(min(n,m)/2) do
    Swap(s, SNew, i, SNew, min(n,m)-i+1);
  end do;

  UNew := Create(s, r, n, integer);
  Copy(s, U, UNew);
  VNew := Create(s, m, m, integer);
  Copy(s, V, VNew);

  reverseCol(s, VNew, m);
  reverseRow(s, UNew, r);

  return SNew, UNew, VNew;
end proc;

# reconstructM(s, S, P, V, n, r)
#
# Input:
#   s, positive integer.
#   S, min(n, r), congruent to the reverse smith form of P over Z(s).
#   P, n x r Mod s matrix.
#   V, unimodular Mod s matrix, with dimension (r x r).
#   There exists some unimodular matrix U with size (n x n) s.t.
#   UPV = S mod s
#
# Output:
#   M, s.t. PV = M Diag(S) mod s
#
reconstructM := proc(s, S, P, V, n, r)
  local Q, i, j, M;
  Q := Create(s, n, r, integer);
  M := Create(s, n, min(r,n), integer);
  Multiply(s, P, V, Q);
  for i from 1 to min(n, r) do
    if S[i] = 0 then next end if;

    for j from 1 to n do
      M[j, i] := modp(iquo(Q[j,i], S[i]), s);
    end do;
  end do;
  
  return M;
end proc;
