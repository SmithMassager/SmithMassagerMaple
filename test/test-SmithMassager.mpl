read "SmithMassager/SmithMassager.mpl";
read "test/test-Helpers.mpl";
read "HermiteForm/genherm.mpl";


checkValidMassagerWrapper := proc(A, n, U := false, M := false, T := false,
S := false)
  global tot;
  ++tot;
  if (U <> false) then
    SS := DiagonalMatrix(S);
    checkValidMassager(A, n, U, M, T, SS);
  else
    recordMatrices(A, n);
  end if;

end proc;

# Check the smithMassager function for odd determinant matrices.
checkOdd := proc()
  local A, ret;
  A := Matrix([[1,0], [0, 3]]);
  ret := SmithMassager(A, 2);
  checkValidMassagerWrapper(A, 2, ret);

  A := Matrix([[1,1], [4, 3]]);
  ret := SmithMassager(A, 2);
  checkValidMassagerWrapper(A, 2, ret);
  A := Matrix([[1,1], [4, 9]]);
  ret := SmithMassager(A, 2);
  checkValidMassagerWrapper(A, 2, ret);

  A := Matrix([[119, 0, 0], [0, 119 * 217, 0], [0, 0, 119 * 217 * 317]]);
  ret := SmithMassager(A, 3);
  checkValidMassagerWrapper(A, 3, ret);

  A := Matrix([[1,0,0], [0, 9, 0], [0, 0, 27]]);
  ret := SmithMassager(A, 3);
  checkValidMassagerWrapper(A, 3, ret);
  A := Matrix([[-13, 10, -20, 27], [27, 30, 15, 30], [0, 15, 15, 6], [-21, 0,
  -15, 9]]);
  ret := SmithMassager(A, 4);
  checkValidMassagerWrapper(A, 4, ret);
end proc;

checkEven := proc()
  local A;
  A := Matrix([[1, 0], [0, 2]]);
  ret := SmithMassager(A, 2);
  checkValidMassagerWrapper(A, 2, ret);

  A := Matrix([[1, 0], [0, 2]]);
  A := randomUnimodularMatrix(2) . A . randomUnimodularMatrix(2);
  ret := SmithMassager(A, 2);
  checkValidMassagerWrapper(A, 2, ret);


  A := Matrix([[2, 0, 0], [0, 2^3 * 3, 0], [0, 0, 2^3 * 3^2 * 11]]);
  A := randomUnimodularMatrix(3) . A . randomUnimodularMatrix(3);
  randomize(10);
  ret := SmithMassager(A, 3);
  checkValidMassagerWrapper(A, 3, ret);

  A := Matrix([[3, 0, 0], [0, 2^3 * 3 * 11, 0], [0, 0, 2^3 * 3^2 * 11 * 117]]);
  A := randomUnimodularMatrix(3) . A . randomUnimodularMatrix(3);
  ret := SmithMassager(A, 3);
  checkValidMassagerWrapper(A, 3, ret);
end proc;

checkRandom := proc()
  local A, n, roll, ret;
  #roll := rand(20..55);
  roll := rand(5..10);
  #for i from 1 to 10000 do
  for i from 1 to 100 do
    n := roll();
    n := nextprime(n);
    #A := randomNonsingular(n);
    A := interestHerm(n);
    ret := SmithMassager(A, n);
    checkValidMassagerWrapper(A, n, ret);
  end do;
end proc;


checkLarge := proc()
  A := Matrix([
    [-1017815489582682768096124315856800979599350 ,

    -83499703026404485333497439554236660593996 ,

    -216964060174020782059856472340773344821956 ,

    1023351860683022463443294819273572115738628],

    [-465389445052925418820542715467584949457670 ,

    -39026689603315303854057453818109423864948 ,

    -100830144515595522621139103810850518260396 ,

    468517777427407095493689510897691942452678],

    [590707025178377499071766107465989286440003 ,

    42776075453454746696806968315458427672144 ,

    115014733269835135211799270161464315255896 ,

    -589914446698171202372388304387598202782637],

    [-3159501274678740886973764884900136730563243 ,

    -256181127660637692225024725778183479461060 ,

    -667709235185318564979681260635896452374104 ,

    3174560178143433955199431730504531586375607]]);

    i := 0;
    randomize(i);
    ret := SmithMassager(A, 4);
    checkValidMassagerWrapper(A, 4, ret);
end proc;

checkSmithMassager := proc()
  checkOdd();
  checkEven();
  checkRandom();
  checkLarge();
  print("unfoundAns:", unfoundAns);
  print("total test cases:", tot);
  print("total unfound:", unfoundCount);
  print("Success %:", 1 - unfoundCount/tot);
end proc;

checkSmithMassager();
