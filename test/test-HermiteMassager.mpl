read "test/test-Helpers.mpl";
read "HermiteForm/HowellHelper.mpl";

checkValidHermiteWrapper := proc(A, n, H := false)
  global tot;
  ++tot;
  if (H <> false) then
    checkValidHermite(A, n, H);
  else
    recordMatrices(A, n);
  end if;

end proc;

checkValidHermite := proc(A, n, H)
  ASSERT(Equal(HermiteForm(A), H));
end proc;

checkRandom := proc()
  local A, n, roll, ret;
  roll := rand(2..30);
  for i from 1 to 100 do
    n := roll();
    A := randomNonsingular(n);
    #ret := SmithMassager(A, n);
    H := HermiteMassager(n, A);
    checkValidHermiteWrapper(A, n, H);
  end do;
end proc;

checkHermiteMassager := proc()
  checkRandom();
  print("unfoundAns:", unfoundAns);
  print("total test cases:", tot);
  print("total unfound:", unfoundCount);
  print("Success %:", 1 - unfoundCount/tot);
end proc;

checkHermiteMassager();
