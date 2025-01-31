(* ::Package:: *)

(* ::Section:: *)
(*Quark-gluon plasma potential*)


(* ::Text:: *)
(*This is a workbook defines the partition function for quark-gluon plasma which existed in the early Universe's first few microseconds. Our temperature range of interest is \!\(TraditionalForm\`\(500\(\ \)\) \**)
(*StyleBox["MeV",*)
(*FontSlant->"Plain"] > \**)
(*SubscriptBox[*)
(*StyleBox["k", "TI"], *)
(*StyleBox["B", "TI"]] \**)
(*StyleBox["T", "TI"] > \(150\(\ \)\) \**)
(*StyleBox["MeV",*)
(*FontSlant->"Plain"]\).  Under these conditions we would like to analyze the magnetization of the quark-gluon plasma and its effect on the chemical potential of the species \!\(TraditionalForm\`\**)
(*StyleBox["i", "TI"] \[Element] \(\**)
(*StyleBox["u", "TI"], \**)
(*StyleBox["d", "TI"], \**)
(*StyleBox["s", "TI"], \**)
(*StyleBox["e", "TI"], \[Nu]\)\). This builds off our prior work in analyzing the magnetized Fermi electron-positron gas.*)


(* ::Subsection:: *)
(*Goals Blueprint*)


(* ::Text:: *)
(*This document sets out to accomplish the following goals:*)


(* ::Item:: *)
(*Implement the standard Fermi-Dirac distribution function.*)


(* ::Subitem:: *)
(*Define the partition function for a single species of fermion (particle and antiparticles) with a shared chemical potential.*)


(* ::Subitem:: *)
(*Define the grand potential for that fermion species.*)


(* ::Item:: *)
(*Evaluate the grand potential for a given energy.*)


(* ::Subitem:: *)
(*Evaluate the massless limit i.e. \!\(TraditionalForm\`\[Epsilon] \[GreaterGreater] \**)
(*StyleBox["m", "TI"]\) such that \!\(TraditionalForm\`\[Epsilon] -> \**)
(*StyleBox["p", "TI"]\).*)


(* ::Subitem:: *)
(*Evaluate the mass corrections which arise in the limit \!\(TraditionalForm\`\**)
(*StyleBox["m", "TI"] \[GreaterGreater] \[Mu]\).*)


(* ::Text:: *)
(*Combining the results from these limits allows us to obtain the grand potential within the regime of \!\(TraditionalForm\`\**)
(*StyleBox["T", "TI"] \[GreaterGreater] \**)
(*StyleBox["m", "TI"] \[GreaterGreater] \[Mu]\). This is specifically the regime where we have a high temperature Fermi gas (non-degenerate regime) but with a nearly matter-antimatter asymmetry. These are the conditions of the early Universe which are very different from the regime found in heavy-ion collisions of matter particles where \!\(TraditionalForm\`\[Mu] \[GreaterGreater] \**)
(*StyleBox["m", "TI"]\).*)


(* ::Item:: *)
(*Once we have the grand potential in terms of mass corrections in the limits appropriate for primordial quark-gluon plasma, we modify the mass to introduce the corrections due to magnetization.*)


(* ::Item:: *)
(*Evaluate and plot the magnetization of the quark-gluon plasma.*)


(* ::Section:: *)
(*Evaluating Fermi Integrals for various limits*)


(* ::Text:: *)
(*We're going to define a set of expressions which allow us to evaluate the Fermi integral for various important limits. To start, we define the Fermi Partition function which we need to integrate over to obtain the free energy. Below we write the fermion partition function and the coefficient for the free energy.*)


(* Define variables and conditions *)
chemicalPotentials={\[Mu],-\[Mu]};
thermodynamicVariables={T,m,\[Mu],p,\[Lambda]};
allThermodynamicVariablesPositive=Positive/@thermodynamicVariables;

(* Define particle energy and partition function *)
particleEnergy[p_,m_]:=Sqrt[m^2+p^2];
fermiPartition[p_,m_,\[Mu]_,T_]:=Log[1+Exp[\[Mu]/T]Exp[-particleEnergy[p,m]/T]];
freeEnergyCoeff[T_,V_,N_]:=-(4\[Pi] T V N)(2\[Pi])^(-3);


(* ::Subsection:: *)
(*Massless Fermi integrals*)


(* ::Text:: *)
(*The product of the two above lines when integrated over the momentum \!\(TraditionalForm\`\**)
(*StyleBox["p", "TI"]\) yields the free energy of the system. We can evaluate the massless limit by integrating with the mass set to zero \!\(TraditionalForm\`\**)
(*StyleBox["m", "TI"] == 0\) and summing over particle(antiparticle) states \!\(TraditionalForm\`\(\[PlusMinus]\[Mu]\)\).*)


(* Evaluate the massless limit with parallelization *)
masslessIntegrand[p_,\[Mu]_,T_]:=p^(2)*fermiPartition[p,0,\[Mu],T];
freeEnergyMassless=
  freeEnergyCoeff[T,V,n]*
   Total[Refine[
     ParallelMap[
       Integrate[masslessIntegrand[p,#,T],{p,0,\[Infinity]}]&,
       chemicalPotentials
     ],allThermodynamicVariablesPositive]
   ];
(* Display the series expansion of the free energy in the massless limit *)
Series[freeEnergyMassless/.\[Mu]->T*Log[\[Lambda]],{\[Lambda],0,4}]//Normal


(* ::Text:: *)
(*This is a check that our method of integration by parts agrees with the above integration.*)


(* Check method of integration by parts with parallelization *)
intByPartsMasslessIntegrand[p_,\[Mu]_,T_] := -(1/3)p^(3)*D[fermiPartition[p,0,\[Mu],T],p];
intByPartsCheck=
freeEnergyCoeff[T,V,n]*
  Total[Refine[
    ParallelMap[
      Integrate[intByPartsMasslessIntegrand[p,#,T],{p,0,\[Infinity]}]&,
      chemicalPotentials
    ],allThermodynamicVariablesPositive]
  ];
(* Display the series expansion of the free energy in the massless limit *)
Series[intByPartsCheck/.\[Mu]->T*Log[\[Lambda]],{\[Lambda],0,4}]//Normal
(* Equality check *)
intByPartsCheck==freeEnergyMassless


(* ::Text:: *)
(*As expected, the above two integration methods yield the same massless limit for the fermion free energy.*)


(* ::Subsection:: *)
(*Massive Fermi integrals (Nonrelativistic limit)*)


(* ::Text:: *)
(*Below is the nonrelativistic limit for \!\(TraditionalForm\`\**)
(*StyleBox["p", "TI"] \[LessLess] \**)
(*StyleBox["m", "TI"]\) for the energy eigenstate and the free energy. After evaluating the integral, we expand for \!\(TraditionalForm\`\[Lambda] \[TildeTilde] 1 \((\[Mu] \[TildeTilde] 0)\)\). These results agree with the results of Landau and Lifshitz, "Statistical Mechanics" (1970). pp. 164. Note that in their derivation \!\(TraditionalForm\`\**)
(*SubscriptBox["\[Mu]", *)
(*StyleBox["NR",*)
(*FontSlant->"Plain"]] \[TildeTilde] \[Mu] - \**)
(*StyleBox["m", "TI"] < 0\) such that \!\(TraditionalForm\`\**)
(*SuperscriptBox[*)
(*StyleBox["e", "TI"], *)
(*SubscriptBox["\[Mu]", *)
(*StyleBox["NR",*)
(*FontSlant->"Plain"]]] \[LessLess] 1\).*)


(* Define the nonrelativistic limit for the energy eigenstate *)
particleEnergySeries[p_,m_]:=
  Normal[Series[Sqrt[m^2+p^2],{p,0,2}]]//
    Refine[#,allThermodynamicVariablesPositive]&;

(* Evaluate the massive limit with integration by parts *)
massiveIntegrandSeries[p_,m_,\[Mu]_,T_]:=-(1/3)p^(3)*
  D[fermiPartition[p,m,\[Mu],T]/.
    particleEnergy[p,m]->particleEnergySeries[p,m],p];
freeEnergyMassiveSeries=
  freeEnergyCoeff[T,V,n]*
   Total[Refine[
     ParallelMap[
       Integrate[massiveIntegrandSeries[p,m,#,T],{p,0,\[Infinity]}]&,
       chemicalPotentials
     ],allThermodynamicVariablesPositive]
   ];

(* Display the series expansion of the free energy in the massive limit *)
freeEnergyMassiveSeries//
  Series[#,{m,0,2}]&//
  Series[#,{\[Mu],0,2}]&//Normal//FullSimplify
freeEnergyMassiveSeries/.\[Mu]->T*Log[\[Lambda]]//
  Series[#,{m,0,2}]&//
  Series[#,{\[Lambda],1,2}]&//Normal//FullSimplify


(* ::Text:: *)
(*We note that above expressions are for the relativistic definition for chemical potential. We have printed the free energy in terms of chemical potential \!\(TraditionalForm\`\[Mu]\) and fugacity \!\(TraditionalForm\`\[Lambda]\) which are not directly comparable to the nonrelativistic quantities. Before we continue, we want to note the distinction between the exponential of a Taylor series and the Taylor series of an exponential. Mathematically, the two should agree up to the order of expansion. However, they may differ if the evaluation of one yields is completed analytically while the other is only evaluated as a series. In other words for \!\(TraditionalForm\`\**)
(*StyleBox["F", "TI"](\**)
(*StyleBox["x", "TI"])\) and \!\(TraditionalForm\`\**)
(*StyleBox["G", "TI"] \((\**)
(*StyleBox["x", "TI"])\) == \**)
(*SuperscriptBox[*)
(*StyleBox["e", "TI"], *)
(*RowBox[{*)
(*StyleBox["F", "TI"], "(", *)
(*StyleBox["x", "TI"], ")"}]]\) then: \!\(TraditionalForm\`\[Integral]\**)
(*StyleBox["d", "TI"] \**)
(*StyleBox["x", "TI"] exp \((\**)
(*StyleBox["F", "TI"] \((\**)
(*StyleBox["a", "TI"])\) + \**)
(*SuperscriptBox[*)
(*StyleBox["F", "TI"], "\[Prime]"] \((\**)
(*StyleBox["a", "TI"])\) \((\**)
(*StyleBox["x", "TI"] - \**)
(*StyleBox["a", "TI"])\) + *)
(*\*FractionBox[\(1\), \(2!\)] \**)
(*SuperscriptBox[*)
(*StyleBox["F", "TI"], "\[Prime]\[Prime]"] \((\**)
(*StyleBox["a", "TI"])\) \**)
(*SuperscriptBox[*)
(*RowBox[{"(", *)
(*RowBox[{*)
(*StyleBox["x", "TI"], "-", *)
(*StyleBox["a", "TI"]}], ")"}], "2"] + \**)
(*StyleBox["O", "TI"] \**)
(*SuperscriptBox[*)
(*RowBox[{"(", *)
(*RowBox[{*)
(*StyleBox["x", "TI"], "-", *)
(*StyleBox["a", "TI"]}], ")"}], "3"])\) - \[Integral]\**)
(*StyleBox["d", "TI"] \**)
(*StyleBox["x", "TI"] \((\**)
(*StyleBox["G", "TI"] \((\**)
(*StyleBox["a", "TI"])\) + \**)
(*SuperscriptBox[*)
(*StyleBox["G", "TI"], "\[Prime]"] \((\**)
(*StyleBox["a", "TI"])\) \((\**)
(*StyleBox["x", "TI"] - \**)
(*StyleBox["a", "TI"])\) + *)
(*\*FractionBox[\(1\), \(2!\)] \**)
(*SuperscriptBox[*)
(*StyleBox["G", "TI"], "\[Prime]\[Prime]"] \((\**)
(*StyleBox["a", "TI"])\) \**)
(*SuperscriptBox[*)
(*RowBox[{"(", *)
(*RowBox[{*)
(*StyleBox["x", "TI"], "-", *)
(*StyleBox["a", "TI"]}], ")"}], "2"] + \**)
(*StyleBox["O", "TI"] \**)
(*SuperscriptBox[*)
(*RowBox[{"(", *)
(*RowBox[{*)
(*StyleBox["x", "TI"], "-", *)
(*StyleBox["a", "TI"]}], ")"}], "3"])\) == 0 + \[Integral]\**)
(*StyleBox["d", "TI"] \**)
(*StyleBox["x", "TI"] \**)
(*StyleBox["O", "TI"] \**)
(*SuperscriptBox[*)
(*RowBox[{"(", *)
(*RowBox[{*)
(*StyleBox["x", "TI"], "-", *)
(*StyleBox["a", "TI"]}], ")"}], "3"]\)*)
(*However, differences may result if the truncated integral is non-convergent and the higher order terms are singular. Care should be taken that the expansion is chosen such that the integration is convergent. As a practical example, the Fermi integral is NOT convergent for expansions in momentum up to \!\(TraditionalForm\`\**)
(*StyleBox["O", "TI"](\**)
(*SuperscriptBox[*)
(*StyleBox["p", "TI"], "4"])\) even though the Fermi integral is convergent for expansions of the microstate energies up to the same order.*)


(* ::Subsection:: *)
(*Mass corrections to relativistic Fermi integrals*)


(* ::Text:: *)
(*Now we want to evaluate the free energy in the form \!\(TraditionalForm\`\[CapitalOmega] \[TildeTilde] \[CapitalOmega]\**)
(*SubscriptBox["|", *)
(*RowBox[{*)
(*StyleBox["m", "TI"], "==", "0"}]]+\[CapitalDelta]\[CapitalOmega]\**)
(*SubscriptBox["|", *)
(*StyleBox["m", "TI"]]\) where \!\(TraditionalForm\`\[CapitalDelta]\[CapitalOmega]\**)
(*SubscriptBox["|", *)
(*StyleBox["m", "TI"]] \[Congruent] \[CapitalOmega]\**)
(*SubscriptBox["|", *)
(*RowBox[{*)
(*StyleBox["m", "TI"], "!=", "0"}]]-\[CapitalOmega]\**)
(*SubscriptBox["|", *)
(*RowBox[{*)
(*StyleBox["m", "TI"], "==", "0"}]]\). Therefore we need to define the massive correction to the free energy \!\(TraditionalForm\`\(\[CapitalDelta]\[CapitalOmega]\**)
(*SubscriptBox["|", *)
(*StyleBox["m", "TI"]]\)\).*)


(* Define the relativistic limit for the energy eigenstate *)
particleEnergyRelativisticSeries[p_,m_]:=
  Normal[Series[Sqrt[m^2+p^2],{m,0,2}]]//
    Refine[#,allThermodynamicVariablesPositive]&;

(* Define the massive correction with integration by parts *)
massiveIntegrandRelativisticSeries[p_,m_,\[Mu]_,T_]:=-(1/3)p^(3)*
  D[fermiPartition[p,m,\[Mu],T]/.
    particleEnergy[p,m]->particleEnergyRelativisticSeries[p,m],p];
massiveIntegrandCorrections[p_,m_,\[Mu]_,T_]:=massiveIntegrandRelativisticSeries[p,m,\[Mu],T]-intByPartsMasslessIntegrand[p,\[Mu],T];
massiveIntegrandCorrections[p,m,\[Mu],T]//
    Refine[#,allThermodynamicVariablesPositive]&;

(* Define the energy correction function *)
massiveFreeEnergyCorrectionsDifferential[p_,m_,\[Mu]_,T_]:=
	freeEnergyCoeff[T,1,1]*
	(-((E^(-(p/T)+\[Mu]/T) p^3)/(3 (1+E^(-(p/T)+\[Mu]/T)) T))+(E^(-((m^2/(2 p)+p)/T)+\[Mu]/T) (1-m^2/(2 p^2)) p^3)/(3 (1+E^(-((m^2/(2 p)+p)/T)+\[Mu]/T)) T));

(* Plot the integrand over the specified range *)
Plot[massiveFreeEnergyCorrectionsDifferential[p,1,0.01,200],{p,0.01,2000},
PlotRange->All,
AxesLabel->{"p [MeV]","\!\(\*TemplateBox[<|\"boxes\" -> FormBox[RowBox[{FractionBox[\"1\", StyleBox[\"V\", \"TI\"]], FractionBox[RowBox[{\"\[PartialD]\", \"\[CapitalOmega]\"}], RowBox[{\"\[PartialD]\", StyleBox[\"p\", \"TI\"]}]]}], TraditionalForm], \"errors\" -> {}, \"input\" -> \"\\\\frac{1}{V}\\\\frac{\\\\partial\\\\Omega}{\\\\partial p}\", \"state\" -> \"Boxes\"|>,\n\"TeXAssistantTemplate\"]\)"},
PlotStyle->Thick]
