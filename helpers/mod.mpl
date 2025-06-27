macro(Copy=LinearAlgebra[Copy]):
macro(Dimensions=LinearAlgebra[Dimensions]):
macro(ModTranspose=LinearAlgebra:-Modular:-Transpose);
macro(Transpose=LinearAlgebra:-Transpose);

# cmod(A, S)
#
#  Input: A - an n x m integer matrix
#         S - an m x m nonsingular diagonal matrix
#
# Output: A cmod S
#             
cmod := proc(A, S, inplace := false)
    local B,i,j,n,m;

    if (type(A, Vector[column])) then
      n:= Dimensions(A);
      m := 1;
    elif (type(A, Vector[row])) then
      m := Dimensions(A);
      n := 1;
    else
      n,m := Dimensions(A);
    end if;

    if (inplace) then
      B := A;
    else
      B := convert(Copy(A), 'Matrix');
    end if;

    if (type(S, Vector[column])) then
      for j to m do for i to n do B[i,j] := modp(B[i,j],S[j]) od od;
    else
      for j to m do for i to n do B[i,j] := modp(B[i,j],S[j,j]) od od;
    end if;

    return B;
end:

# lowerCmod(A, S)
#
#  Input: A - an n x m integer matrix
#         S - an m x m nonsingular diagonal matrix
#
# Output: A'
#  A' = perform A cmod S on only the lower triangular part of A (including diagnols)
#             
lowerCmod := proc(A, S, inplace := false)
    local B,i,j,n,m;

    if (type(A, Vector[column])) then
      n:= Dimensions(A);
      m := 1;
    elif (type(A, Vector[row])) then
      m := Dimensions(A);
      n := 1;
    else
      n, m := Dimensions(A);
    end if;

    if (inplace) then
      B := A;
    else
      B := convert(Copy(A), 'Matrix');
    end if;

    if (type(S, Vector)) then SS := S;
    else SS := Vector(m, (i) -> S[i, i]) end if;

    for j to m do for i to n do if (i >= j) then B[i,j] := modp(B[i,j],SS[j]) end if od od;

    return B;
end:

#  Input: A - an n x m integer matrix
#         S - an n x n nonsingular diagonal matrix
#
# Output: A rmod S
#
# TODO: Buggy for vector inplace updates.
rmod := proc(A, S, inplace := false)
    local B,i,j,n,m;

    if (type(A, Vector[column])) then
      n:= Dimensions(A);
      m := 1;
    elif (type(A, Vector[row])) then
      m := Dimensions(A);
      n := 1;
    else
      n,m := Dimensions(A);
    end if;

    if (inplace) then
      B := A;
    else
      B := convert(Copy(A), 'Matrix');
    end if;

    for j to m do for i to n do B[i,j] := modp(B[i,j],S[i,i]) od od;

    return B;

end:

# Define: The compact representation of A cmod S to be A', S', sNum.
#         Where S is a n x n diagonal matrix nonsingular with ascending order. A is n x n matrix.
#         sNum is the number of nontrivial entries of S (i.e. S[i,i] != 1).
#         S' represents the last sNum nontrivial entry.
#         A' represents the last sNum nontrivial column.

# Similarly for A rmod S.

# Define: The compact representation of H to be H', V, hNum.
#         Where H is a n x n matrix in Hermite form.
#         hNum represent the number of nontrivial columns (i.e. H[i,i] != 1).
#         H' is a n x hNum matrix containing the nontrivial columns of H.
#         V is the indices for those hNum nontrivial columns of H, i.e. col(H, V[i]) = col(H', i).

# compactHermiteModMult(H, idices, numIdx, v, n, s)
#
# Input:
#   H, idices, numIdx in compact representation of a n x n Hermite form matrix.
#   v, n column vector.
#   s, positive integer.
#
# Output: w
#   w = Hv mod s
#
compactHermiteModMult := proc(H, idices, numIdx, v, n, s)
  local res, j, i;
  res := Vector(n);
  j := 1;

  for i to n do
    #if (j <= numIdx and i <> 1 and H[idices[j], j] >= s) then
    #  H[1..idices[j], j] := modp(H[1..idices[j], j], s);
    #end if;

    if (j <= numIdx and i = idices[j]) then
      ModAddMultiple(s, v[i], res, H[1..-1, j], res);
      ++j;
    else
      res[i] := modp(res[i] + v[i], s);
    end if;
  end do;
  return res;
end proc;

# cmodCompactHermiteMult(H, idices, numIdx, v, n, S, SNum)
#
# Input:
#   H, idices, numIdx in compact representation of a n x n Hermite form matrix.
#   A, n x SNum Matrix i.e. cmod S matrix.
#   S, SNum, compact representation of a diagonal ascending n x n matrix.
#
# Output: W in compact representation wrt cmod S.
#   W = HA cmod S
#
cmodCompactHermiteMult := proc(H, idices, numIdx, A, n, S, SNum)
  local res, i;
  res := Matrix(n, SNum);
  for k from SNum to 1 by - 1 do
    res[1..-1, k] := compactHermiteModMult(H, idices, numIdx, A[1..-1, k], n, S[k]);
    #j := 1;
    #sk := S[k];
    #for i to n do
    #  if (j <= numIdx and i = idices[j]) then
    #    ModAddMultiple(sk, A[i,k], res, 1..-1, k, H, 1..-1, j, res, 1..-1, k);
    #    ++j;
    #  else
    #    res[i, k] := modp(res[i, k] + A[i, k], sk);
    #  end if;
    #end do;
  end do;

  return res;
end proc;

cmodCompactHermiteMult1 := proc(H, idices, numIdx, A, n, S, SNum)
  local res, i;
  res := Matrix(n, SNum);
  for k from SNum to 1 by - 1 do
    res[1..-1, k] := compactHermiteModMult(H, idices, numIdx, A[1..-1, k], n, S[k]);
    #j := 1;
    #sk := S[k];
    #for i to n do
    #  if (j <= numIdx and i = idices[j]) then
    #    ModAddMultiple(sk, A[i,k], res, 1..-1, k, H, 1..-1, j, res, 1..-1, k);
    #    ++j;
    #  else
    #    res[i, k] := modp(res[i, k] + A[i, k], sk);
    #  end if;
    #end do;
  end do;

  return res;
end proc;


# compactRmodScalarMultiply(v, e, SEntry, SNum)
#
# Input:
#   v, Snum vector that is a compact representation of a n vector.
#   e, scalar.
#   SEntry, SNum, compact representation of a diagonal ascending n x n matrix.
#
# Output: w in compact representation.
#   w = ev rmod S
#
compactRmodScalarMultiply := proc(v, e, SEntry, SNum)
  local res;
  res := Vector(SNum);

  for i to SNum do
    res[i] := modp(e * v[i], SEntry[i]);
  end do;
  return res;
end proc;

# compactRmodAdd(v, w, SEntry, SNum)
#
# Input:
#   v, Snum vector that is a compact representation of a n vector.
#   w, Snum vector that is a compact representation of a n vector.
#   SEntry, SNum, compact representation of a diagonal ascending n x n matrix.
#
# Output: u in compact representation.
#   u = v + w rmod S
#
compactRmodAdd := proc(v, w, SEntry, SNum)
  res := Vector(SNum);
  for i to SNum do
    res[i] := modp(v[i] + w[i], SEntry[i]);
  end do;
  return res;
end proc;

# cmodCompactMult(A, m, B, n, SEntry, SNum)
#
# Input:
#   A, m x n matrix.
#   B, n x SNum matrix that is a compact representation of B' cmod S. 
#   SEntry, SNum, compact representation of a diagonal ascending n x n matrix, S.
#
# Output: C m x SNum in compact representation.
#   C = AB cmod S.
#
cmodCompactMult := proc(A, m, B, n, SEntry, SNum)
  local res, i;
  res := Matrix(m, SNum);
  for i to SNum do
    ModMultiply(SEntry[i], A, B, 1..-1, i, res, 1..-1, i);
  end do;

  return res;
end proc;


# rmodCompactMult(A, m, B, n, SEntry, SNum)
#
# Input:
#   A, SNum x m matrix that is a compact representation of a n vector.
#   B, m x n matrix that is a compact representation of B' cmod S. 
#   SEntry, SNum, compact representation of a diagonal ascending n x n matrix, S.
#
# Output: C SNum x n in compact representation.
#   C = AB rmod S.
#
rmodCompactMult := proc(A, n, B, m, SEntry, SNum)
  return Transpose(cmodCompactMult(Transpose(B), m, Transpose(A), n, SEntry, SNum));
end proc;

