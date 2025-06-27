# See the paper: A fast algorithm for computing the Smith normal form with 
#                multipliers for a nonsingular integer matrix.

read("config.mpl");

macro(Transpose=MTM:-transpose);
macro(Dimension=LinearAlgebra:-Dimension);
macro(Swap=LinearAlgebra:-Modular:-Swap);
macro(Map=Threads:-Map);
macro(Seq=Threads:-Seq);
macro(Create=Threads:-Create);
macro(Wait=Threads:-Wait);

# Swap rows and columns inplace.
SwapRow := proc(A, r1, r2)
  Swap(2, A, r1, A, r2);
end proc;

SwapCol := proc(A, c1, c2)
  Swap(2, A, 1..-1, c1, A, 1..-1, c2);
end proc;

# Length(x)
#
# Input:
#   x, some integer.
#
# Output: l
#   l := length(x)
#
Length := proc(x)
  local bits, a, absx;

  if (x = 0) then return 1; end if;

  bits := 0;
  a := 1;
  absx := abs(x);
  while a <= absx do
    ++bits;
    a *= 2;
  end do;

  return bits;
end proc;

# VecLength(v)
#
# Input:
#   v, vector of integers.
#
# Output: l
#   l := length(||v||)
#
VecLength := proc(v)
  return Length(max(abs(v)));
end proc;

# AvgColLength(A, n, m)
#
# Input:
#   A, n x m integer matrix
#
# Output: l, vLen
#   l := (length(A_1) + ... + lengt(A_m))/m, where A_i is ith col of A.
#   vLen := [length(A_1), ..., length(A_m)]
#
AvgColLength := proc(A, n, m)
  local res, vLen, i;
  vLen := Vector(m);
  res := 0;
  for i from 1 to m do
    vLen[i] := VecLength(A[1..-1, i]);
    res += vLen[i];
  end do;
  return ceil(res/m), vLen;
end proc;

# AvgRowLength(A, n, m)
#
# Input:
#   A, n x m integer matrix
#
# Output: l
#   Where l := sum_{i=1}^{m} length(A_i)/m, where A_i is ith row of A.
#
AvgRowLength := proc(A, n, m)
  return AvgColLength(Transpose(A), m, n);
end proc;


# getColPLParam(A, n, m)
#
# Input:
#   A, n x m integer matrix
#
# Output: d, e
#   d, avg col length of A.
#   e, s.t. the col partial LP of A has entry <= 2^d and dimension grows at most by m.
#      See Theorem 27.
#
getColPLParam := proc(A, n, m)
  local e, d, i, vLen;
  e := Vector(m);
  d, vLen := AvgColLength(A, n, m);
  d := max(d, 64);

  for i from 1 to m do
    e[i] := ceil(vLen[i]/d) - 1;
  end do;

  return d, e;
end proc;

# Similar to getColPLParam but wrt to rows now.
getRowPLParam := proc(A, n, m)
  return getColPLParam(Transpose(A), m, n);
end proc;

# Ed(e, d)
#
# Input:
#   e, d, positive integer.
#
# Output: e x 1 vector,
#   [ -2^d ]
#   [      ]
#   [  ..  ]
#   [      ]
Ed := proc(e, d)
  local ret;
  ret := Vector(e);
  if (e <> 0) then ret[1] := -2^d end if;
  return ret;
end proc;

# Fd(e, d)
#
# Input:
#   e, d, positive integer.
#
# Output: e x e vector, let x := -2^d ,
#   [ 1               ]
#   [ x  1            ]
#   [    x   1        ]
#   [     .    .      ]
#   [       .    .    ]
#   [             1   ]
#   [             x 1 ]
#
Fd := proc(e, d)
  local ret, x;
  x := -1 * (2^d);
  ret := Matrix(e, e, (i, j) -> if (i = j) then 1 elif (i - j = 1) then x else 0 end if);
  return ret;
end proc;

# getbAdic(a, e, b)
#
# Input:
#   a, integer.
#   e, non-negative integer.
#   d, positive integer.
#
# Output: v of length e+1 s.t.
#   each entry in v is mod 2^d except the last entry.
#   a == v[0] + v[1]2^d + ... + v[e-1]2^((e-1)d) + v[e]2^(ed)
#
getbAdic := proc(a, e, b)
  local v, r, bAdic, n;
  v := Vector(e+1);
  r := 0;

  v[e+1] := iquo(a, b^e, 'r');
  bAdic := convert(r, 'base', b);
  n := nops(bAdic);
  v[1..n] := convert(bAdic, 'Vector')[1..n];

  return v;
end proc;

# VecPL(v, n, e, d)
#
# Input:
#   v, integer vector, length n.
#   e, non-negative.
#   d, positive.
#
# Pre: e must be large enough to satisfy the output conditon.
#
# Output: c, C
#   v := c + col(C, 1)*2^d + col(C, 2)*2^(2d) + ... + col(C, e)*2^(ed)
#   c, n x 1 vector mod 2^d.
#   CC, n x e Matrix mod 2^d.
#
VecPL := proc(v, n, e, d)
  local c, CC, b, i;
  c := Vector(n);
  CC := Matrix(n, e);

  b := 2^d;

  VecPLHelper := proc(i)
    local bAdic, tmp, m;
    if (v[i] = 0) then return end if;
    if (e = 0) then 
      c[i] := v[i];
      return;
    end if;
    #bAdic := getbAdic(v[i], e, b);
    bAdic := Vector(e+1);
    r := 0;

    bAdic[e+1] := iquo(v[i], b^e, 'r');
    tmp := convert(r, 'base', b);
    m := nops(tmp);
    bAdic[1..m] := convert(tmp, 'Vector')[1..m];

    c[i] := bAdic[1];
    CC[i, 1..e] := bAdic[2..-1];
  end proc;

  map(VecPLHelper, [seq(1..n)]);

  return c, CC;
end proc;

# getPosPart(v, n)
#
# Input:
#   v, n x 1 vector
#
# Output: v'
#   v', is v but only with the negative entries set to 0.
#
getPosPart := proc(v, n)
  local res;
  res := Vector(n, i -> if(v[i] > 0) then v[i] else 0 end if);
  return res;
end proc;

# getNegPart(v, n)
#
# Input:
#   v, n x 1 vector
#
# Output: v'
#   v', is v but only with the positive entries set to 0, and negative entries *-1.
#
getNegPart := proc(v, n)
  return getPosPart(-v, n)
end proc;

# ColPLHelper(A, n, m, e, d)
#
# Input:
#   A, n x m integer matrix.
#   e, vector of integers.
#   d, positive integer.
#
# Output: D,
#   D, is the col partial linearization of A wrt to parameter e and d.
#      See corollary 25.
#
ColPLHelper := proc(A, n, m, e, d)
  local s, partialSum, DD, updateColPLHelper;
  s := add(e);
  DD := Matrix(n + s, m + s);
  # Calculate: partialSum[i] = e[1] + .. + e[i-1]
  partialSum := Vector(m+1);
  for i from 2 to m+1 do
    partialSum[i] := partialSum[i-1] + e[i-1];
  end do;

  updateColPLHelper := proc(i)
    local c, CC, F, E;
    c, CC := VecPL(A[1..-1, i], n, e[i], d);
    E := Ed(e[i], d);
    F := Fd(e[i], d);

    DD[1..n, i] := c;
    DD[1..n, m+1+partialSum[i]..m+partialSum[i+1]] := CC;
    DD[n+1+partialSum[i]..n+partialSum[i+1], i] := E;
    DD[n+1+partialSum[i]..n+partialSum[i+1], m+1+partialSum[i]..m+partialSum[i+1]] := F;
  end proc;

  h := [seq(1..m)];
  #Map[tasksize=TASKSIZE](updateColPLHelper, h);
  map(updateColPLHelper, h);

  return DD, n+s, m+s;
end proc;

# Similar to ColPLHelper but wrt to rows.
#
RowPLHelper := proc(A, n, m, e, d)
  local ret, nbar, mbar, DD;
  DD, mbar, nbar := ColPLHelper(Transpose(A), m, n, e, d);
  return Transpose(DD), nbar, mbar;
end proc;

# ColPL(A, n, m)
#
# Input:
#   A, n x m integer matrix.
#
# Output: D,
#   D, is the col partial linearization of A wrt to parameter e and d.
#      See corollary 25 and Theorem 27.
#
ColPL := proc(A, n, m)
  local e, d;
  d, e := getColPLParam(A, n, m);
  return ColPLHelper(A, n, m, e, d);
end proc;

# Similar to ColPL but wrt to rows.
#
RowPL := proc(A, n, m)
  local e, d;
  d, e := getRowPLParam(A, n, m);
  return RowPLHelper(A, n, m, e, d);
end proc;


# greaterTuple(x,y)
#
# Input:
#   x, y, tuple of integers.
#
# Output: x > y elementwise compairson.
#
greaterTuple := proc(x, y)
  if (type(x, integer)) then 
    evalb(x > y);
  end if;

  if (x[1] > y[1]) then 
    return true;
  elif (x[1] = y[1]) then 
    return greaterTuple(x[2..-1], y[2..-1]);
  else 
    return false;
  end if
end proc;

# getPerm(A, n)
#
# Input:
#  A, n x n integer matrix.
#
# Output: A', P, Q
#   A' = PAQ, where P, Q are permutation matrices.
#   A' have the property that the length A[i,i] bounds the length of submatrix A[i..n, i..n].
#
getPerm := proc(A, n)
  local a, elimR, elimC, RMap, CMap, idx, P, Q, RInvMap, CInvMap, i, j, r, c;
  # RMap, CMap is mapping col/row of A to col/row of the modified A.
  # RInvMap, CInvMap does the opposite.
  RMap := Vector(n, i -> i);
  CMap := Vector(n, i -> i);
  RInvMap := Vector(n, i -> i);
  CInvMap := Vector(n, i -> i);

  # elimR, elimC records the rows and cols in A that are not used from previous iteartions.
  elimR := {};
  elimC := {};
  P := Vector(n, i -> i);
  Q := Vector(n, i -> i);
  a := Array(1..n*n);

  for i to n do
    for j to n do
      a[(i-1)*n + j] := [Length(A[i,j]), i, j];
    end do;
  end do;

  sort[inplace](a, greaterTuple);

  idx := 1;
  for i to n do
    while(a[idx][2] in elimR or a[idx][3] in elimC) do ++idx end do;
    r := a[idx][2];
    c := a[idx][3];
    elimR := elimR union {r};
    elimC := elimC union {c};
    SwapRow(P, i, RMap[r]);
    SwapRow(Q, i, CMap[c]);
    SwapRow(RMap, RInvMap[i], r);
    SwapRow(CMap, CInvMap[i], c);
    SwapRow(RInvMap, i, RMap[r]);
    SwapRow(CMap, i, CMap[c]);
  end do;
  
  return P, Q;
end proc;

# getPermPLParam(A, n)
#
# Input:
# 	A, n x n integer matrix
#
# Output: d, e, eBar
# 	d, e, eBar, specified in Corollary 30.
#
getPermPLParam := proc(A, n)
  local e, d, i, v, a, idx, P, Q, ebar;
  e := Vector(n);
  v := Vector(n);

  P, Q := getPerm(A, n);

  for i to n do
    v[i] := Length(A[P[i], Q[i]]);
  end do;

  d := ceil(add(v)/n);

  for i to n do
    e[Q[i]] := ceil(v[i]/d) - 1;
  end do;

  eBar := Vector(n + add(e));
  for i to n do
    eBar[P[i]] := e[Q[i]];
  end do;

  return d, e, eBar;
end proc;

# PermPL(A, n, m)
#
# Input:
#   A, n x m integer matrix.
#
# Output: D,
#   D, is the permutation partial linearization of A wrt to parameter d, e, eBar.
#      See section 4.2.
#
PermPL := proc(A, n);
  local d, e, ebar, nbar, mbar, colPL;
  d, e, ebar := getPermPLParam(A, n);
  colPL, nbar, mbar := ColPLHelper(A, n, n, e, d);
  return RowPLHelper(colPL, nbar, mbar, ebar, d);
end proc;
