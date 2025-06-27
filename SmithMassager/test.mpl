read "SmithMassager/SNFMod.mpl";
read "SmithMassager/reverseSmithForm.mpl";
read "helpers/basic.mpl";
read "helpers/naiveCert.mpl";
read "helpers/fastIntCert.mpl";
read "helpers/fastUniCert.mpl";
read "helpers/radix.mpl";
read "helpers/multiplyHelper.mpl";
read "helpers/time.mpl";
read "helpers/mod.mpl";
read "macros.mpl";
read "SmithMassager/projectionBasis.mpl";
read "SmithMassager/SmithMassager.mpl" ;
read "helpers/partial.mpl";

macro(AddMultiple=LinearAlgebra:-Modular:-AddMultiple);
macro(Concatenate=ArrayTools:-Concatenate);
macro(DiagonalMatrix=LinearAlgebra:-DiagonalMatrix);
macro(Equal=LinearAlgebra:-Equal);
macro(Identity=LinearAlgebra:-Modular:-Identity);
macro(MatrixInverse=LinearAlgebra:-MatrixInverse);
macro(Mod=LinearAlgebra:-Modular:-Mod);
macro(Multiply=LinearAlgebra:-Modular:-Multiply);
macro(RandomMatrix=LinearAlgebra:-RandomMatrix);
macro(Random=LinearAlgebra:-Modular:-Random);
macro(SubMatrix=LinearAlgebra:-SubMatrix);
macro(HermiteForm=LinearAlgebra:-HermiteForm);
macro(SmithForm=LinearAlgebra:-SmithForm);
macro(Random=LinearAlgebra:-Modular:-Random);

# computeProjection(s, T, S, Sstar, n, size)
#
# Input: T is cmod S
#        Sstar = s S^{-1}
#        s = S[-1]
#        S, Sstar are diagonal (n x n)
#        NOTE: All input matrices are (n x n)
# 
# Output: Compute modp(T Sstar P, s), where P is (n x size) matrix from Z(s) picken u@r.
#
computeProjection := proc(s, T, S, Sstar, n, size)
  ret := Create(s, n, size, integer);
  if (S[1, 1] <> 1) then
    Multiply(s, T, 1..-1, 1, Sstar[1,1]/Sstar[2,2] * Random(S[1,1], 1, size, integer), 1 , 1..-1, ret);
  end if;
  for i from 2 to n-1 do
    if (S[i, i] = 1) then next end if;
    AddMultiple(s, ret, Multiply(s, T, 1..-1, i, Random(S[i, i], 1, size, integer), 1, 1..-1), ret);
    Multiply(s, Sstar[i,i]/Sstar[i+1,i+1], ret, ret);
  end do;
  if (S[n, n] <> 1) then
    AddMultiple(s, ret, Multiply(s, T, 1..-1, n, Random(S[n, n], 1, size, integer), 1, 1..-1), ret);
    Multiply(s, Sstar[n,n], ret, ret);
  end if;
  return ret;
end proc;

naiveComputeProjection := proc(s, T, S, Sstar, n, size)
  P := RandomMatrix(n, size, generator=0..s-1);
  E := MultiplyHelper(s, MultiplyHelper(s, T, Sstar), P);
  return E;
end proc;

computeProjectionHelper := proc(s, T, S, Sstar, n, size)
  print(ceil(log2(s)));
  if (s < 2^500) then
    return naiveComputeProjection(s, T, S, Sstar, n, size);
  else
    return computeProjection(s, T, S, Sstar, n, size);
  end if;
end proc;

DiagInvHelper := proc(s, A, n)
  return Matrix(n, n, (i, j) -> if (i=j) then s/A[i,i] else 0 end if);
end proc;

helper := proc(TT, SS, n)
  global timesPostUpdate, timesComputeProjection, timesSNFModS;
  local P1, U, S, V, U1, M, F, UNew, VNew, v, p, u, m;

  S := Vector(n);
  s := SS[-1, -1];
  snew := s;
  U := Matrix(n, n);
  M := Matrix(n, n);
  #Sstar := s . MatrixInverse(SS);
  Sstar := DiagInvHelper(s, SS, n);
  T := Matrix(TT);
  m := Vector(n);
  #size := floor(n/3);
  size := 1;
  k := 10;
  p := 1;
  it := 1;
  while p < n do
    break if s = 1;
    dt := getDT(s);
    size := min(n - p + 1, 2^it);
    it := it + 1;
    if (dt = float[8]) then
      size := n-p+1;
    end if;
    #print(size);
    #print(s);
    #print(Sstar);
    #print(SmithForm(T[1..-1, 1..-p] . Sstar[1..-p, 1..-1]));
    #print(SmithForm(T[1..-1, 1..-1] . Sstar[1..-1, 1..-1]));
    #E := modp(T[1..-1, 1..-p] . Sstar[1..-p, 1..-p] . RandomMatrix(n-p+1, size+k, generator=0..s-1), s);
    #P := RandomMatrix(n-p+1, k, generator=0..s-1);
    #B := Matrix(n-p+1, size);
    #for i from 1 to size do B[-i, -i] := 1; end do;
    #for j from 1 to min(size+k,n-p+1) do for i from 1 to j do if (i = j) then P[-i,-j] := 1 else P[-i,-j] := 0 end if end do end do;
    #P := ArrayTools:-Concatenate(2, B, P);
    #print(P);
    #print(SmithForm(T[1..-1, 1..-1] . Sstar[1..-1, 1..-1] . P));
    #print(T);
    #print(Sstar);
    #P := RandomMatrix(n, size+k, generator=0..s-1);
    #E := MultiplyHelper(s, MultiplyHelper(s, T, Sstar), P);
    start := time[real]();
    E := computeProjectionHelper(s, T, SS, Sstar, n, size+k);
    endTime := time[real]();
    timesComputeProjection := [op(timesComputeProjection), endTime-start];

    #E := modp(T[1..-1, 1..-1] . Sstar[1..-1, 1..-1] . P, s);
    #E := modp(T[1..-1, 1..-p] . Sstar[1..-p, 1..-p] . P, s);
    print("calling SNF");
    start := time[real]();
    ret := SNF(E, n, size, k, s);
    endTime := time[real]();
    timesSNFModS := [op(timesSNFModS), endTime-start];
    print("finished");
    if (ret <> false) then
      Snew, Unew, Mnew, Qnew := ret;
    else 
      #print("iteration: ", it);
      return false;
    end if;
    start := time[real]();
    S[-p-size+1..-p] := Snew;
    #print(Snew);
    U[-p-size+1..-p, 1..-1] := Unew;
    M[1..-1, -p-size+1..-p] := Mnew;
    #print(p + size);
    #print(n);
    if (p + size > n) then break; end if;
    QInv := LinearAlgebra:-Modular:-Inverse(s, Qnew);
    #if (s < 2^500) then
      #MQUT := MultiplyHelper(s, Mnew, MultiplyHelper(s, QInv, LeftSparseMult(Unew, T, size, n, n)));
    #else
      #MQUT := cmodmul(Mnew, cmodmul(QInv, cmodmul(Unew, T, SS), SS), SS);
      start1 := time[real]();
      MQUT := cmodmulHelper(Mnew, cmodmulHelper(QInv, cmod(LeftSparseMult(Unew, T, size, n, n), SS), SS), SS);
      end1 := time[real]();
      print("MQUT time: ", end1 - start1);
    #end if;
    start1 := time[real]();
    AddMultiple(s, s-1, T, MQUT, T);
    #modp(add_maple(T, MQUT, s-1, s), s);
    end1 := time[real]();
    print("AddMultiple time: ", end1 - start1);

    #Sstar[1..-1, 1..-1] := Sstar * Snew[1]/s;
    s := SS[-p-size+1,-p-size+1];
    #Sstar[1..-1, 1..-1] := Sstar / min(Sstar[-p-size+1, -p-size+1], s/Snew[1]);
    Sstar[1..-1, 1..-1] := Sstar / Sstar[-p-size+1, -p-size+1];
    #Sstar[1..-1, 1..-1] := Sstar / (s/Snew[1]);
    #s := max(Snew[1], SS[-p-size+1,-p-size+1]);
    #s := Snew[1];
    #print(Sstar);
    #print(SS);
    #s := Snew[1];
 
    for i from 1 to p+size-1 do
      #lprint(T[1..-1, -i]*Sstar[-i, -i]);
      T[1..-1, -i] := modp(T[1..-1, -i]*Sstar[-i, -i], s);
      #T[1..-1, -i] := modp(T[1..-1, -i], s);
      Sstar[-i, -i] := 1;
    end do;
    for i from 1 to p+size-1 do
      #SS[-i, -i] := Snew[1];
      SS[-i, -i] := s;
    end do;
    T := cmod(T, SS);
    #print("SS", SS);
    #print("T", T);
    p := p + size;
    endTime := time[real]();
    timesPostUpdate := [op(timesPostUpdate), endTime-start];
  end do;
  return S, U, M;
end proc;

randomInput := proc(n)
  local SS, T, A;
  roll := rand(1..10);
  roll1 := rand(1..2);
  SS := Matrix(n, n);
  SS[1,1] := 2;
  for i from 2 to n do
    if (roll1() = 1) then
      SS[i,i] := SS[i-1,i-1] * 2;
    else 
      SS[i,i] := SS[i-1,i-1] * roll();
    end if;
  end do;
  A := RandomMatrix(n, n, density=0.2, generator=1..1000);
  T := HermiteForm(Concatenate(1, A, SS))[1..n, 1..n];
  return T, SS;
end proc;

Success := 0;
Total := 0;
check := proc(T, S, n)
  global Total, Success;
  #T, S := randomInput(n);
  s := S[-1, -1];
  Stmp := SmithForm(T . s . MatrixInverse(S));
  Sret, Uret, Mret := helper(T, S, n);
  SS := Vector(n);
  Sans := Vector(n);
  for i from 1 to n do
    Sans[i] := Stmp[i, i];
    if (Sret[i] = 0) then
      SS[-i] := s;
    else 
      SS[-i] := s/Sret[i];
    end if;
  end do;
  Total := Total + 1;
  if (Equal(SS, Sans)) then Success := Success + 1;
  else print("not equal") end if;
end proc;

#n := nextprime(1007);
#for i from 1 to 10 do
#  n := nextprime(100);
#  T, S := randomInput(n);
#  check(T, S, n);
#end do;
#print("Success prob", Success/Total);

#resetprofile();
#profile(helper, naiveComputeProjection, computeProjection);
#exprofile();
#n := nextprime(30);
#T, S := randomInput(n);
#helper(T, S, n);
S := ImportMatrix("S-1000");
T := ImportMatrix("B-1000");
helper(T, S, 1000);
#showprofile();
printGlobalTime();
