-- Change input here:

nMin = 1
nMax = 5
gMin = 4
gMax = 5
prefix = "Hodge_";

-- End of input area.

loadPackage "HodgeIntegrals"
R = hodgeRing(max(1, gMax), nMax)

--
-- PREPARE PRINT FUNCTION
--

printList = (f, A) -> (
  f << "(";
  len = length(toList(A));
  for i from 1 to len do(
    f << A#(i-1);
    if i != len then (
      f << ", ";
    );
  );
  f << ")"
);

printPartition = (f, g, L) -> printList(f, for i from 1 to g list number(L, j -> j == i));

printRational = (f, r) -> (
  if instance(r, ZZ) then (
    f << r << " // 1";
  ) else if instance(r, QQ) then (
    f << numerator(r) << " // " << denominator(r);
  ) else if isConstant(r) then (
    cc = coefficient(1_R, r);
    f << numerator(cc) << " // " << denominator(cc);
  ) else (
    f << "ERROR: Output of integral(...) not constant: " << r;
  );
);

printHodge = (f, g, n, A, L, r) -> (
  f << g << ", " << n << ", ";
  printList(f, A);
  f << ", ";
  printPartition(f, g, L);
  f << "; ";
  printRational(f, r);
  f << endl;
);

--
-- START ITERATING
--

for g from gMin to gMax do (
  for n from nMin to nMax do (

    dimM = 3*g - 3 + n;
    -- If the following file name format is changed, we need to revise Hodge_integrals.jl, in
    -- particular the function _load_Hodge_integrals(...).
    f = concatenate(prefix, "g_", toString(g), " n_", toString(n), ".txt") << "";

    -- a1 + ... + an = cut while l1 + ... + g lambda_g = dimM - cut.
    for cut from 0 to dimM do (
    
      for A in unique(apply(compositions(n, cut), comp -> rsort(comp))) do (
      
        psiProd = 1_R;
        for i from 1 to n do (
          psiProd = psiProd * psi_i^(A#(i-1));
        );

        for L in partitions(dimM - cut, g) do (
        
          fullProd = psiProd;
          for l in L do (
            fullProd = fullProd * lambda_l;
          );

          r = integral(g, n, fullProd);
          printHodge(f, g, n, A, L, r);
        );



      );

    );

    f << close;

  );
);