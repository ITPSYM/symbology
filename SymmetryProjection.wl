(* ::Package:: *)

BeginPackage["SymmetryProjection`"];

GenerateSymmetryProjection::usage =
  "GenerateSymmetryProjection[repRule, letterRep, rootDef] generates the \
symmetry projection matrix from a variable replacement rule.\n\nrepRule: \
list of rules like {s12->s23, s23->s34, ...}\nletterRep: {W[1]->expr1, \
..., W[n]->exprn}\nrootDef: {eps5->Sqrt[Delta5]} or similar\n\nReturns a \
SparseArray matrix M where Log[W[i]_new] = Sum_j M[[i,j]]*Log[W[j]].";

Options[GenerateSymmetryProjection] = {deBug -> False, "specialkinematics" -> {}};

Begin["`Private`"];

(* ============================================================ *)
(* Helper: numeric point sampler (from user)                    *)
(* ============================================================ *)

Options[GenerateNumReal] = {"specialkinematics" -> {}};
GenerateNumReal[sq_, ovar_, OptionsPattern[]] := Module[{var, rep, tem, k = 1},
  var = ovar /. OptionValue["specialkinematics"] // Variables;
  rep = Thread@Rule[var, RandomChoice[{1, -1}]*RandomPrime[{10^5, 10^6}, Length[var]]];
  tem = sq /. {Power[a_, 1/2] :> a} /. OptionValue["specialkinematics"];
  While[Length[Values[rep] // Union] < Length[var] || ((tem /. rep) < 0 && k < 1000),
    rep = Thread@Rule[var, RandomChoice[{1, -1}]*RandomPrime[{10^5, 10^6}, Length[var]]];
    k = k + 1;
  ];
  Return[Thread@Rule[ovar, ovar /. OptionValue["specialkinematics"] /. rep]];
];

(* ============================================================ *)
(* Helper: find integer linear relation among Log[letters]      *)
(* ============================================================ *)

Options[FindLetterLinearRelation] = {deBug -> False, "specialkinematics" -> {}};
FindLetterLinearRelation[basis_, target_, OptionsPattern[]] := Module[{gr, var, sq, rep, rel},
  gr = Append[basis, target];
  var = gr /. {Log[z_] :> z} // Variables;
  sq = Cases[{gr[[-1]]}, Power[_, 1/2], Infinity] // DeleteDuplicates;
  If[OptionValue[deBug], Print["square root: ", Short[sq, 20]]];
  rep = Table[GenerateNumReal[sq[[1]], var, "specialkinematics" -> OptionValue["specialkinematics"]], {i, 1, 3}];
  rel = Commonest[Quiet[FindIntegerNullVector[ReplaceAll[gr /. {Log[z_] :> Log[Abs[z]]}, #], 30]] & /@ rep];
  Return[rel];
];

(* ============================================================ *)
(* Main: variable replacement -> symmetry projection matrix     *)
(* ============================================================ *)

GenerateSymmetryProjection[repRule_, letterRep_, rootDef_, OptionsPattern[]] := Module[
  {n, alphabet, letterExprs, logBasis, transformedExprs, result, i, j,
   exprNew, found, sqLetters, basisPairs, pair, pi, rel, cc, b1, b2,
   debug, deltaExpr, deltaNew},

  debug = OptionValue[deBug];
  n = Length[letterRep];
  alphabet = letterRep[[All, 1]];            (* {W[1], W[2], ..., W[n]} *)
  letterExprs = letterRep[[All, 2]] /. rootDef; (* expressions with eps5 -> Sqrt[Delta5] *)
  logBasis = Log /@ alphabet;                (* {Log[W[1]], ..., Log[W[n]]} *)

  (* Check Delta5 invariance *)
  deltaExpr = (eps5 /. rootDef)^2;           (* Delta5 *)
  deltaNew = deltaExpr /. repRule;
  If[!PossibleZeroQ[deltaNew - deltaExpr],
    Print["Warning: Delta5 is not invariant under the replacement rule."];
    Print["  Delta5     = ", Short[deltaExpr, 30]];
    Print["  Delta5_new = ", Short[deltaNew, 30]];
  ];

  (* Compute transformed expressions *)
  transformedExprs = Table[letterExprs[[i]] /. repRule, {i, n}];

  (* Identify square-root letters *)
  sqLetters = Select[Range[n], !FreeQ[letterExprs[[#]], Power[_, 1/2]] &];
  If[debug, Print["Square-root letters: ", sqLetters]];

  (* For each letter, find the transformation in Log space *)
  result = Table[0, {n}];

  For[i = 1, i <= n, i++,
    exprNew = transformedExprs[[i]];
    found = False;

    (* --- Step 1: Direct match (exprNew == ±letterExprs[[j]]) --- *)
    (* Log[±W[j]] = Log[W[j]] + const, so coefficient is +1 *)
    For[j = 1, j <= n && !found, j++,
      If[Simplify[exprNew - letterExprs[[j]]] === 0 ||
         Simplify[exprNew + letterExprs[[j]]] === 0,
        result[[i]] = logBasis[[j]];
        found = True;
        If[debug, Print["W[", i, "] -> +W[", j, "]"]];
      ];
    ];

    (* --- Step 2: Inverse match (exprNew == ±1/letterExprs[[j]]) --- *)
    (* Log[±1/W[j]] = -Log[W[j]] + const, so coefficient is -1 *)
    If[!found,
      For[j = 1, j <= n && !found, j++,
        If[Simplify[exprNew * letterExprs[[j]] - 1] === 0 ||
           Simplify[exprNew * letterExprs[[j]] + 1] === 0,
          result[[i]] = -logBasis[[j]];
          found = True;
          If[debug, Print["W[", i, "] -> -W[", j, "]"]];
        ];
      ];
    ];

    (* --- Step 3: FindLetterLinearRelation for square-root letters --- *)
    (* Try pairs of square-root letters as basis *)
    If[!found && !FreeQ[exprNew, Power[_, 1/2]],
      basisPairs = Subsets[sqLetters, {2}];
      For[pi = 1, pi <= Length[basisPairs] && !found, pi++,
        pair = basisPairs[[pi]];
        {b1, b2} = pair;
        rel = FindLetterLinearRelation[
          {Log[letterExprs[[b1]]], Log[letterExprs[[b2]]]},
          Log[exprNew],
          deBug -> debug,
          "specialkinematics" -> OptionValue["specialkinematics"]
        ];
        (* rel = {c1, c2, c3} with c1*Log[b1] + c2*Log[b2] + c3*Log[target] = 0 *)
        (* => Log[target] = -(c1/c3)*Log[b1] - (c2/c3)*Log[b2] *)
        If[rel =!= $Failed && ListQ[rel] && Length[rel] == 3 && rel[[3]] =!= 0,
          cc = -rel/rel[[3]];  (* {-(c1/c3), -(c2/c3), -1} -> take first two *)
          result[[i]] = cc[[1]] * logBasis[[b1]] + cc[[2]] * logBasis[[b2]];
          found = True;
          If[debug, Print["W[", i, "] -> ", cc[[1]], "*W[", b1, "] + ", cc[[2]], "*W[", b2, "]"]];
        ];
      ];
    ];

    (* --- Step 4: Extended search — try all pairs if still not found --- *)
    If[!found && !FreeQ[exprNew, Power[_, 1/2]],
      If[debug, Print["W[", i, "]: extended search over all pairs..."]];
      For[b1 = 1, b1 <= n && !found, b1++,
        For[b2 = b1 + 1, b2 <= n && !found, b2++,
          rel = FindLetterLinearRelation[
            {Log[letterExprs[[b1]]], Log[letterExprs[[b2]]]},
            Log[exprNew],
            deBug -> debug,
            "specialkinematics" -> OptionValue["specialkinematics"]
          ];
          If[rel =!= $Failed && ListQ[rel] && Length[rel] == 3 && rel[[3]] =!= 0,
            cc = -rel/rel[[3]];
            result[[i]] = cc[[1]] * logBasis[[b1]] + cc[[2]] * logBasis[[b2]];
            found = True;
            If[debug, Print["W[", i, "] -> ", cc[[1]], "*W[", b1, "] + ", cc[[2]], "*W[", b2, "] (extended)"]];
          ];
        ];
      ];
    ];

    If[!found,
      Print["Warning: no transformation found for W[", i, "], expr = ", Short[exprNew, 30]];
      result[[i]] = 0;
    ];
  ];

  (* Build projection matrix: M[[i,j]] = coeff of Log[W[j]] in result[[i]] *)
  projMat = CoefficientArrays[result, logBasis][[2]];
  Return[projMat];
];

End[];
EndPackage[];
