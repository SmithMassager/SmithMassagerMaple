read "macros.mpl";

macro(Multiply=LinearAlgebra:-Modular:-Multiply);

#
# MultiplyHelper(m, A, B, dt)
#
# Input:
#   m, positive integer
#   A, B mod m matrix
#   dt, datatype
#
# Output: modp(A . B, m)
#
fastMultiplyHelper := proc(m, A, B)
  local AA, BB, dt;

  dt := getDT(m);

  if (dt <> float[8] and not (type(A, 'Vector') or type(B, 'Vector'))) then
    #return modp(imlMult(A, B), m);
    return modp(fmpzMult(A, B), m);
  end if;

  AA := Matrix(A, datatype=dt);
  BB := Matrix(B, datatype=dt);

  return ifelse(dt = float[8], map(round, Multiply(m, AA, BB)), Multiply(m, AA, BB));
end proc;

naiveMultiplyHelper := proc(m, A, B)
  AA := Matrix(A, datatype=integer);
  BB := Matrix(B, datatype=integer);
  return Multiply(m, AA, BB);
end proc;

#
# LeftSparseMult(A, B, n, m, l)
#
# Input:
#   n, m, l positve integers
#   A, n x m integer matrix. A should be sparse.
#   B, m x l integer matrix
#
# Output: A . B
#
LeftSparseMult := proc(A, B, n, m, l, s := 0)
  local res;
  res := Matrix(n, l);
  AA := A;
  BB := B;
  if (n = 1) then AA := Matrix(A); end if;
  if (l = 1) then BB := Matrix(B); end if;

  for i from 1 to n do
    for j from 1 to m do
      if (AA[i,j] <> 0) then
	res[i, 1..-1] += AA[i,j] * BB[j, 1..-1];
	if (s <> 0) then
	  res := mods(res, s);
	end if;
      end if;
    end do;
  end do;

  return res;
end proc;
