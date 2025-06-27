interface(quiet=true);
read "HermiteForm/HowellHelper.mpl";

for i to 5 do
  read cat("SmithMassager/smithMassager", i, "-new.out");
  n, m := Dimensions(A);
  resetprofile();
  profile(HermiteDiagonals, compactHcol,
  SpecialHowellTransform, CompactScaledMatVecProd, compactRmodScalarMultiply, compactRmodAdd, rmodCompactMult,
  HermiteViaHowell, CompactScaledMatVecProd, cmodCompactHermiteMult,
  cmodCompactHermiteMult1
  );
  HermiteMassagerHelper(n, M, S);
  showprofile();
end do;
