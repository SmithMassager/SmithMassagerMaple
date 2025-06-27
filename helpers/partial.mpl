macro(Dimensions=LinearAlgebra[Dimensions]):
macro(ColumnDimension=LinearAlgebra[ColumnDimension]):
macro(Mod=LinearAlgebra[Modular][Mod]):
macro(Multiply=LinearAlgebra[Multiply]):
read "helpers/mod.mpl":

cmodmulHelper := proc(A, B, F)
  #return cmod(plMultiply(A, B), F);
  return cmodMulColPL(A, B, F);
end proc;

#
#  Input:  F - an o x o diagonal matrices with positive entries
#          B - an m x o integer matrix with B = cmod(B,F)
#          A - an n x m integer matrix with nonnegative entries
#
# Output: cmod(A B,F)
#
cmodmul := proc(A,B,F)
   local n,m,X,e,f,P,i,j,k,AL,BE,BL,s,mx,AJ,ed,z,e_sum,BJ,fd,f_sum,CL,C,L,num,Y,tmp,o,CLf,dtype;
   
   n,m := Dimensions(A);
   o := ColumnDimension(B);

   # Note: Choice of modulo X and datatype are somewhat arbitrary
   #       For datatype=integer any X>=2 is OK.  
   #       For datatype in {integer[8],float[8]} need X small enough.
   #      

   # Choose modulus X
   # The choice 2^16 is motivated by the fact that
   #  
   #      2^16 is two bytes (4 X-adic cofficients for 64 bit integer)
   #      Thus maybe quick to get X-adic coefficients of GMP mpz_t.
   #      and to convert X-adic expansion back to mpz_t.
   #
   #      Largest m with m (2^16-1)^2 <= 2^53 -1 yields m < 2 million,
   #      so result of dot products should fit into the float mantissa
   #      for abitrarily large matrices
   #
   X := 2^39:
   dtype := float[8];
   #X := 2^100:
   #dtype := integer;

   # Get linearization paramaters e=(e1,..,em) 
   #print("Get linearization paramaters e=(e1,..,em)");
   e := Vector(m);
   AJ := Vector(m);
   s,ed := 1,0;
   for i to m do 
       mx := max(seq(A[z,i],z=1..n))+1;
       P := 1; 
       for e[i] from 0 while P < mx do P := P*X od; 
       ed := ed + e[i];
       AJ[i] := s..ed;
       s := s+e[i];
   od;
   e_sum := add(e[z],z=1..m);

   # Construct partial linearization of A
   #print("Construct partial linearization of A",n,e_sum);
   AL := Matrix(n,e_sum,datatype=dtype);
   for i to m do 
      if e[i]=0 then 
         next
      elif e[i]=1 then
         AL[1..-1,AJ[i]] := Matrix(n,e[i],A[1..-1,i],datatype=dtype); 
      else
         AL[1..-1,AJ[i]] := Matrix(n,e[i],map(convert,[seq(A[z,i],z=1..n)],base,X),datatype=dtype); 
      fi;
   od;
   print("AL", AL);

   # Construct expansion of B
   #print("Construct expansion of B",e_sum,o);
   BE := Matrix(e_sum,o);
   for i to m do
       if e[i]=0 then next fi;
       s := add(e[z],z=1..i-1);
       for j to o do BE[s+1,j] := modp(B[i,j],F[j,j]) od;
       for k from 2 to e[i] do
          for j to o do BE[s+k,j] := modp(X*BE[s+k-1,j],F[j,j]) od;
       od;
   od;
   print("Be", BE);
   
   # Get linearization paramaters f=(f1,...,fm)
   #print("Get linearization paramaters f=(f1,...,fm)");
   f := Vector(o);
   BJ := Vector(o);
   s,fd := 1,0;
   for i to o do 
       mx := max(seq(BE[z,i],z=1..e_sum))+1;
       P := 1; 
       for f[i] from 0 while P < mx do P := P*X od; 
       fd := fd + f[i];
       BJ[i] := s..fd;
       s := s+f[i];

      #P := 1; 
      #for f[i] from 0 while P<F[i,i] do P := P*X od; 
      #fd := fd + f[i];
      #BJ[i] := s..fd;
      #s := s + f[i];
   od;
   f_sum := add(f[z],z=1..o);

   # Construct partial linearization of BE
   #print("Construct partial linearization of BE",e_sum,f_sum);
   BL := Matrix(e_sum,f_sum,datatype=dtype);
   for i to o do 
      if f[i]=0 then 
         next 
      elif f[i]=1 then
         BL[1..-1,BJ[i]] := Matrix(e_sum,f[i],BE[1..-1,i],datatype=dtype);
      else
         BL[1..-1,BJ[i]] := Matrix(e_sum,f[i],map(convert,[seq(BE[z,i],z=1..e_sum)],base,X),datatype=dtype);
      fi;
   od;

   lprint("BL", BL);
   print("BL", BL);
   #print("AL:",Dimensions(AL),"BL:",Dimensions(BL));

   # Do product
   #print("Do product");
   CL := Multiply(AL,BL);
   print("CL", CL);
   
   # Compress back
   #print("Compress back");
   C := Matrix(n,o);
   for j to o do 
      if f[j] = 0 then next fi;
      for i to n do 
         L := map(ceil,[seq(CL[i,k],k=BJ[j])]);
         num,Y := nops(L),X;
         while num>1 do
            tmp := seq(L[2*z-1]+L[2*z]*Y,z=1..floor(num/2));
            if type(num,even) then L := [tmp] else L := [tmp,L[-1]] fi;
           num,Y:= ceil(num/2),Y^2;
         od;
         C[i,j] := modp(op(L),F[j,j]);
      od;
   od;
   
   #print("From",n,m,o,"To",Dimensions(AL),ColumnDimension(BL));
  
   return C;

end:

#with(LinearAlgebra):
#n := 1000:
#m := 1000:
#o := 1000:
#F := DiagonalMatrix([seq(1112221,j=1..o-1),1112221^(floor(o/50))]):
#f := rand(0..F[o,o]):
#A := map(abs,RandomMatrix(n,m)):
#for i to n do A[i,m] := f() od:
#B := RandomMatrix(m,o):
#for i to m do B[i,o] := f() od:
#B := cmod(B,F):
#print("Doing mul");
#C1 := cmod(A.B,F):
#print("Doing cmodmul");
#C2 := cmodmul(A,B,F):
#ArrayTools[IsEqual](C1,C2);
A := Matrix(4, 3, [174834750381,    354694120149,    1071471131108,
	    359397246838,    294769546183,      59265239510,
	    203378573253,    823099973867,     864314217739,
	    501330370363,    335388017177,     649171192121]);

B := Matrix(3,4, [155104849056,    187495169546,    825349823857,    1028426420281,
 312588803210,    383645827161,    729272400364,     185251129945,
    422706148896,    836789973312,    137615175382,     224031165462]);

#A := Matrix(4, 3, [
#186300155373410424057423,    809752245738409921813824,    471140202548962213029411,
#770308930516590027558466,     43702733057225910732464,    516589905810855916338210,
#206581292101041027114646,    771547136512582194805116,    436869002880112199654615,
#522493890169331510297197,    423185157151460671796098,    297075648780560562175426
#]);
#B := Matrix(3, 4, [
#465567712523,     183890914282,    553878770491,    206486682327,
#743388979582,    1077069463355,    137895274931,    450954183126,
#772524068084,     584127404289,    843708932012,    719563942545
#]);
mx := max(A . B):
F := Matrix(4, 4, (i, j) -> mx):

cmodmul(A, B, F);

#43296307643824213501477 1130305200449246346171553 703964104873426102649143 1079837531178409018889240 804458824879477327798075 802138307540900208217137 573673739543066951780359 120012274384783527725995 45022151820652663082170 484325930353226026092997 813666216737553974662182 1196795648157559437485620344226705267607452028 1 00000000000000436869002880112199654615 00000000000000522493890169331510297197 0 00000000000000423185157151460671796098 0 00000000000000297075648780560562175426
#[43296307643824213501477, 1130305200449246346171553, 703964104873426102649143, 1079837531178409018889240, 804458824879477327798075, 802138307540900208217137, 573673739543066951780359, 120012274384783527725995, 45022151820652663082170, 484325930353226026092997, 813666216737553974662182, 1196795648157559437485620344226705267607452028, 1, 00000000000000436869002880112199654615, 00000000000000522493890169331510297197, 0, 00000000000000423185157151460671796098, 0, 00000000000000297075648780560562175426]
