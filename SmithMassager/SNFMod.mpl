read "helpers/basic.mpl";
read "helpers/time.mpl";
kernelopts(assertlevel=1);
macro(Dimension=LinearAlgebra:-Dimension);
macro(Multiply=LinearAlgebra:-Modular:-Multiply);
macro(AddMultiple=LinearAlgebra:-Modular:-AddMultiple);
macro(Identity=LinearAlgebra:-Modular:-Identity);

#mymodp := proc(a, b)
#  print(b);
#  print(a);
#  return modp(convert(a, integer), b);
#end proc;
#
#myquo := proc(a, b)
#  return iquo(convert(a, integer), convert(b, integer));
#end proc;
#
#
#
#macro(modp=mymodp);
#macro(igcd=gcd);
#macro(iquo=myquo);

# MatrixGCD(s, A, n, m)
#
# Input:
#  A, n x m mod s matrix.
#
# Output: GCD(A)
#
#
MatrixGCD := proc(s, A, n, m)
  local g, i, j;
  g := s;

  if (m = false) then
    for i to n do
      g := igcd(A[i], g);
      if (g = 1) then return g end if;
    end do;
  else
    for i to n do
      for j to m do
	g := igcd(A[i,j], g);
	if (g = 1) then return g end if;
      end do;
    end do;
  end if;

  return g;
end proc;

# createDiag(n, m, a)
#
# Input:
#   n, m positive integer
#   a, integer
#
# Output:
#   Output a matrix of n x m with shape:
#   [a      ]
#   [ a     ]
#   [  a    ]
#   [   a   ]
#
createDiag := proc(n, m, a := 1)
  return Matrix(n, m, (i,j) -> if i = j then a else 0 end if);
end proc;

# reScale(a, s)
#
# Input:
#   a, mod s integer.
#   s, positive integer.
#
# Output: c, u, g
# Such that (a * e) = (a, s). Where e = u + c * s/g and e is coprime to s.
#
# Post: None
#
# Proof: Note if gcd(a, s) = 1 then we are done. So suppose gcd(a, s) != 1.
# By EEA: ua + vs = g. Consider (u, s/g, s) = 1, by definition of g = gcd(s,a).
# For c to exists we need to prove gcd(u, s/g) = 1.
# u(a/g) + v(s/g) = 1 <-> (a/g) u + v (s/g) = 1. So by EEA u and s/g are
# coprime.
# a * e = a * (u + c s/g) = au + cs (a/g) = g mod s.
#
reScale := proc(a, s)
  local g, u, v;
  g := igcdex(a, s, 'u', 'v');
  return Stab(u, s/g, s), u, g;
end proc;

# reScaleEntry(A, n, m, V, s, i)
#
# Input:
#   A, an n x m Mod s matrix
#   V, an m x m Mod s matrix, unimodular.
#   s, positive integer.
#   i, 1 <= i <= min(m, n) integer.
#
#   Further assume that V records all the unimodular column operations applied
#   to A.
#   A has shape:
#   [s_1  *  *     *  *]           [s_1  *  *     *  *]
#   [    ... *     *  *]           [    ... *     *  *]
#   [        a_i   *  *]           [        s_i   *  *]
#   [         *    *  *]   ->      [         *    *  *]
#   [         ...  *  *]           [         ...  *  *]
#   [         *    *  *]           [         *    *  *]
#   Where all s_j/gcd(s_j, s) = 1 for j < i.
#
# Output: None
#
# Post:
#   Modify A and V accordingly s.t. A_ii/gcd(A_ii, s) = 1 with V still
#   unimodular in Z(s).
#
reScaleEntry := proc(A, n, m, V, s, i)
  local c, u, g, e;
  c, u, g := reScale(A[i,i], s);
  e := modp(u + c * s / g, s);
  Multiply(s, e, A, 1..-1, i, A, 1..-1, i);
  Multiply(s, e, V, 1..-1, i, V, 1..-1, i);
end proc;


# extractGCD(v, w, s)
#
# Input:
#   v, w, vector Mod s of same size.
#   s, positive integer.
#
# POST: None.
#
# Output: c
#   s.t. gcd(v + cw, s) = gcd(v, w, s).
#
# Follow this description from section 5: T. Mulders and A. Storjohann.
# The modulo N extended gcd problem for polynomials.
# Idea: Compute a diagonal matrix T = V[v, w]. Where V is unimodular. And
# T looks like:
# [a g]
# [* 0]
# [* 0]
# [...]
# [* 0]
#
# With g = gcd(w). (This done via stab.)
# Compute c s.t. (a + cg, s) = (a, g, s).
# Then gcd(T, s) = gcd(T_1 + cT_2, s). Since unimodular matrices preserve gcd
# We have gcd([v, w], s) = gcd(v + cw, s).
#
extractGCD := proc(v, w, s)
  local dV, dW, m, c, i, g;
  dV := Dimension(v);
  dW := Dimension(w);

  if (dV <> dW) then
    ERROR(sprintf("dimension of given v and w are not the same size")):
  end if;

  g := MatrixGCD(s, w, dW, false);

  # Find a and g given above.
  m := Vector([igcd(v[1], s), igcd(w[1],s)]);
  for i from 2 to dV do
    if (m[2] = g) then break end if;
    c := Stab(m[2], w[i], s);
    AddMultiple(s, c, m, Vector([v[i], w[i]]), m);
  end do;

  return Stab(m[1], m[2], s);
end proc;

# Return a matrix of size half of A s.t. gcd(A, s) = gcd(A', s)
extractHalfMatrixGCD := proc(s, A, n, m, l, r)
  local B, c;
  mid := ceil((l+r)/2);
  colNum := mid-l+1;
  B := Matrix(n, colNum);

  for i from 0 to mid-l do
    if (mid+i+1 > r) then break end if;
    c := extractGCD(A[1..-1, l+i], A[1..-1, mid+i+1], s);
    AddMultiple(s, c, A, 1..-1, l+i, A, 1..-1, mid+i+1, A, 1..-1, l+i);
  end do;

  return colNum;
end proc;

# extractMatrixGCDHelper(s, A, n, m, idx)
#
# Input:
#   A, n x m mod s matrix.
#   n, m, s, positive integer.
#   idx, 1 <= idx <= min(m,n) an integer.
#
# Output: v, g s.t. A'v = gcd(A'), g = A'v.
#   A' = A[idx..-1, idx..-1]
#
extractMatrixGCDHelper := proc(s, A, n, m, idx := 1)
  local v, w, rG, ci, i, g, b;
  rG := idx..-1;
  g := A[rG, idx];
  b := MatrixGCD(s, A[idx..-1, idx..-1], n-idx+1, m-idx+1);
  v := Vector(m - idx + 1);
  v[1] := 1;

  for i from idx+1 to m do
    if (MatrixGCD(s, g, n-idx+1, false) = b) then break end if;
    w := A[rG, i];
    ci := extractGCD(g, w, s);
    AddMultiple(s, ci, g, w, g);
    v[i-idx+1] := ci;
  end do;

  return v, g;
end proc;

# extractMatrixGcd(s, S, V, n, m, idx)
#
# Input:
#   S, n x m mod s matrix.
#   V, m x m mod s matrix, unimodular.
#   n, m, s, positive integer.
#   idx, 1 <= idx <= min(m,n) an integer.
#
# Assume V records the column operations applied to S.
#
# Post:
# Consider the rectangular submatrix of S formed from (idx, idx) to (n, m),
# call it S'.
# Modify S s.t. the idx-th column vector v, gcd(v) = gcd(S').
# Record the column operation in V.
# If o is row instead then the idx-th row vectors v, gcd(v) = gcd(S').
#
# Output: None.
#
extractMatrixGCD := proc(s, S, V, n, m, idx)
  local v, g;
  v, g := extractMatrixGCDHelper(s, S, n, m, idx);
  V[idx..-1, idx] := v;
  S[1..-1, idx] := S . V[1..-1, idx];
end proc;

# vectorGCD(s, w, n)
#
# Input:
#   w, nx1 mod s vector
#   s, positive integer
#
# Output: v
#   <v, w> = gcd(w)
#
vectorGCD := proc(s, w, n)
  local v, rG, ci, i, g;
  g := modp(w[1], s);
  v := Vector[row](n);
  v[1] := 1;

  for i from 2 to n do
    if (w[i] = 0) then next end if;
    if (g = 1) then break end if;
    ci := Stab(g, w[i], s);
    g := modp(g + ci * w[i], s);
    v[i] := ci;
  end do;

  return v;
end proc;

# extractRowGcd(s, S, V, n, m, idx)
#
# Input:
#   S, n x m mod s matrix.
#   U, m x n mod s matrix.
#   n, m, s, positive integer.
#   idx: 1 <= idx <= min(m,n) an integer.
#
# Pre:
#   Assume U records the row operations applied to S.
#   Let the submatrix from (idx, idx) to (n, m) be S'.
#   The first column of S' records the gcd(S', s).
#
# Post:
#   Modify S s.t. gcd(S_(idx,idx)) = gcd(S').
#   Record the row operation in U.
#
# Output: None.
#
extractRowGcd := proc(s, S, U, n, m, idx)
  local v, w, rG, ci, i, g;
  v := vectorGCD(s, S[idx..-1, idx], n-idx+1);
  U[idx, idx..-1] := v;
  S[idx, 1..-1] := Multiply(s, U, idx, 1..-1, S);
end proc;

# eliminateCol(A, U, n, m, s, idx)
#
# Input:
#   A, n x m Mod s matrix
#   U, m x n Mod s matrix, unimodular.
#   s, positive integer.
#   idx, 1 <= idx <= n, integer.
#
# Post:
#   Modify: A.
#   Assuming U records the row operations corresponding to the upper triangular
#   part applied on A.
#   A also records itself for lower triangular row operations. (See
#   SmithFormModS for more details).
#   Assume that entry A_{idx, idx} is the ith smith invariant factor under
#   Z(s).
#   Apply row operations to A s.t. A have the following shape (ignoring the
#   recorded lower triangular row operations, i.e. imagine the 'x' are 0's):
#
#   [s_1  *  *     *  *]         [s_1  *  *     *  *]
#   [ x  ... *     *  *]         [ x  ... *     *  *]
#   [ x   x  s_i   *  *]   ->    [ x   x  s_i   *  *]
#   [ x   x   *    *  *]         [ x   x   x    *  *]
#   [ x   x   ...  *  *]         [ x   x   ...  *  *]
#   [ x   x   *    *  *]         [ x   x   x    *  *]
#
#   The recorded lower triangular row operations are represented by 'x'.
#
# Output: None.
#
eliminateCol := proc(A, U, n, m, s, idx)
  local i, x;

  if (A[idx, idx] <> 0) then
    for i from idx+1 to n do
      x := modp(-1 * A[i, idx]/A[idx, idx], s);
      AddMultiple(s, x, A, i, 1..-1, A, idx, 1..-1, A, i, 1..-1);
      A[i,idx] := x;
    end do;
  end if;
end proc;

# SmithFormModS(A, n, m, s)
#
# Input:
#
#   A, an n x m Mod s matrix
#
# Output: S, U, V
#
#   S, a Vector of length r = min(n,m)
#   U, the top r x n rows of a unimodular Mod s matrix, say W
#   V, an m x m Mod s unimodular matrix
#
#   Note: W A V = Diag(S[1],S[2],...,S[r]) is in Smith form mod s 
#
# Idea:
#   Start with U, V = I(n x n), I(m x m).
#   Apply column operations s.t. the first column of A contains the gcd of the
#   whole matrix A w.r.t. s.
#   Similarly apply row operations so that entry (1,1) now contains the gcd of
#   the whole matrix of A w.r.t. s.
#   Now the (1,1) entry may not be normalized. So we rescale it and record it
#   in V, let this rescaling factor be r.
#   Thus we have:
#
#   [ 1 *  *  *]    [ r      ]       [ s_1 *  *]
#   [   1      ]    [ * 1    ]       [  *  *  *]
#   [     1    ] A  [ *  ... ]   =   [  *  *  *]
#   [      ... ]    [ *     1]       [  *  *  *]
#   [         1]                     [  *  *  *]
#
#   Note by definition of S_1 the entries in column 1 and row 1 are all
#   divisible by s_1.
#   We apply the row operations to zero out the first column.
#     [ 1        ]  [ s_1 *  *]     [ s_1 *  *]
#     [ x 1      ]  [  *  *  *]     [     *  *]
#     [ x   1    ]  [  *  *  *]  =  [     *  *]
#     [ x    ... ]  [  *  *  *]     [     *  *]
#     [ x       1]  [  *  *  *]     [     *  *]
#   NOTE: That this operation is lower triangular. Thus we can record these
#   entries in A it self. Thus our updated matrices are:
#        [ s_1 *  *]      [ 1 *  *  *]      [ r      ]
#        [ x   *  *]      [   1      ]      [ * 1    ]
#   A =  [ x   *  *], U = [     1    ], V = [ *  ... ]
#        [ x   *  *]      [      ... ]      [ *     1]
#        [ x   *  *]      [         1]
#
#   Now we repeat this process iteratively for the submatrix starting at (2,2)
#   to (n, m). Continue this process we will get A into the following form:
#         [ s_1 *   *]
#         [ x  s_2  *]
#    A =  [ x   x s_3]
#         [ x   x   x]
#         [ x   x   x]
#   We can then extract the lower triangular out to get the actual matrix
#   U required.
#   To get the actual S and V just apply column operations to zero out all the
#   off diagonals.
#
SmithFormModS := proc(A, n, m, s)
  local S, U, V, c, w, i, Up, J, Vp, SVec, x;

  S := Matrix(n, m, A);
  V := Identity(s, m, integer);
  U := createDiag(m, n);
  r := min(n,m);

  for i from 1 to r do
    # Apply row and column operations s.t. GCD(S_ii) = GCD(S').
    # Where S' is the submatrix of S from (i, i) to (n, m).
    extractMatrixGCD(s, S, V, n, m, i);
    extractRowGcd(s, S, U, n, m, i);
    reScaleEntry(S, n, m, V, s, i);
    eliminateCol(S, U, n, m, s, i);
  end do;


  # Reconstruct the matrix U: Recall that S embeds the lower triangular part of
  # U. So extract that part out and multiply with U to get the left unimodular
  # matrix.
  Up := createDiag(r, m);

  for i from 1 to r do
    Up[i, i] := 1;
    for j from i+1 to r do
      Up[j, i] := S[j, i];
      S[j, i] := 0;
    end do;
  end do;

  U := Multiply(s, Up, U);
  #for i from r+1 to n do
  #  U[i,i] := U[i,i] + 1;
  #end do;

  # Eliminate all the off diagonals of the upper triangular matrix S. And
  # update the right unimodular matrix V.
  Vp := Identity(s, m, integer);
  SVec := Vector(r);
  for i from 1 to r do
    SVec[i] := S[i,i];
    for j from i+1 to m do
      if (S[i, i] <> 0) then
        x := modp(-1 * S[i, j]/S[i, i], s);
        AddMultiple(s, x, Vp, 1..-1, j, Vp, 1..-1, i, Vp, 1..-1, j);
        S[i,j] := 0;
      end if;
    end do;
  end do;
  V[1..-1, 1..-1] := Multiply(s, V, Vp);

  return SVec, U, V;
end proc;
