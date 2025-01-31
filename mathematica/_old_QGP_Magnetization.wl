(* ::Package:: *)

(* ::Chapter:: *)
(*Quark-gluon plasma potential*)


(* ::Text:: *)
(*This is a workbook which defines the partition function for quark-gluon plasma which existed in the early Universe's first few microseconds. Our temperature range of interest is 500 MeV>Subscript[k, B]T>150 MeV. Under these conditions we would like to analyze the magnetization of the quark-gluon plasma and its effect on the chemical potential of the species i\[Element]u,d,s,e,\[Nu]. This builds off our prior work in analyzing the magnetized Fermi electron-positron gas.*)


(* ::Section:: *)
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
(*Evaluate the massless limit i.e. \[Epsilon]\[GreaterGreater]m such that \[Epsilon]\[RightArrow]p.*)


(* ::Subitem:: *)
(*Evaluate the mass corrections which arise in the limit m\[GreaterGreater]\[Mu].*)


(* ::Text:: *)
(*Combining the results from these limits allows us to obtain the grand potential within the regime of T\[GreaterGreater]m\[GreaterGreater]\[Mu] . This is specifically the regime where we have a high temperature Fermi gas (non-degenerate regime) but with a nearly matter-antimatter asymmetry. These are the conditions of the early Universe which are very different from the regime found in heavy-ion collisions of matter particles where \[Mu]\[GreaterGreater]m.*)


(* ::Item:: *)
(*Once we have the grand potential in terms of mass corrections in the limits appropriate for primordial quark-gluon plasma, we modify the mass to introduce the corrections due to magnetization.*)


(* ::Item:: *)
(*Evaluate and plot the magnetization of the quark-gluon plasma.*)


(* ::Section:: *)
(*Initial Setup*)


(* ::Subsection:: *)
(*Definitions*)


(* ::Text:: *)
(*We first define the Fermi distribution.*)


EnergyFull[p_,m_] := Sqrt[m^2+p^2];
BoltzFactor = (1/T)*(EnergyFull[p,m]-\[Sigma]*\[Mu]);
Fermi[BoltzFactor_] := 
	(Exp[BoltzFactor]+1)^(-1);


(* ::Text:: *)
(*The grand potential is  given by the integration over phase-space of the partition function defined by the Fermi distribution. Using (Fetter and Walecka, 1971) we write*)


Omega[V_,g_,IntFermi_] := 
	-(1/3)*((g*V)/(2 Pi^2))*IntFermi;


(* ::Text:: *)
(*The variable "IntFermi" is the integration of the Fermi distribution over momentum.*)


IntFermi[p_,BoltzFactor_] := 
	Sum[Integrate[p^(3)*Fermi[BoltzFactor],{p,0,\[Infinity]}],{\[Sigma],{-1,1}}];


(* ::Subsection:: *)
(*Massless limit*)


(* ::Text:: *)
(*We evaluate the Fermi integral in the limit \[Epsilon]\[RightArrow]p to obtain the massless limit.*)


BoltzFactorZeroMass = (1/T)*(EnergyFull[p,0]-\[Sigma]*\[Mu]);
IntegratedMassZero = 
	IntFermi[p,BoltzFactorZeroMass]//FullSimplify[#,T>0&&\[Mu]>0] &


(* ::Text:: *)
(*The above is the correct expression which is proportional to the massless grand potential.*)


(* ::Section:: *)
(*Old Scratch Notes on Fermi Integrals*)


(* ::Text:: *)
(*Now we turn our focus to obtaining mass corrections. These corrections are determined by considering the size of the mass relative to the chemical potential. Below, let us examine well - known limits and reproduce them first to ensure our understanding is correct. Let us take the limit of p\[RightArrow]0 and then integrate the non-relativistic result.*)


LimitP = Series[BoltzFactor,{p,0,2}]//Normal[#]&//Refine[#,p>0&&m>0]&
LimitM = Series[BoltzFactor,{m,0,2}]//Normal[#]&//Refine[#,p>0&&m>0]&
IntLimitP = 
	IntFermi[p,LimitP]//FullSimplify[#,T>0&&m>0&&\[Mu]>0] &


Series[p/EnergyFull[p,m],{m,0,2}]//Normal//Refine[#,p>0]&
p/(Series[EnergyFull[p,m],{m,0,2}]//Normal)//Refine[#,p>0]&


(* ::Text:: *)
(*Here we then expand around m\[TildeTilde]0. This is the same expansion as m/\[Mu]\[TildeTilde]0.*)


Series[IntLimitP,{m,0,2}]//Normal[#]&//FullSimplify[#,T>0&&m>0&&\[Mu]>0]&
Series[IntLimitP,{T,\[Infinity],1}]//Normal[#]&//FullSimplify[#,T>0&&m>0&&\[Mu]>0]&
Series[IntLimitP,{T,0,4}]//Normal[#]&//FullSimplify[#,T>0&&m>0&&\[Mu]>0]&


(* ::Text:: *)
(*This is the correct expansion which gives us the mass corrections to the grand potential. Let us now instead explore the opposite expansion where mass is large and chemical potential is small i.e. m\[RightArrow]\[Infinity] or \[Mu]/m\[TildeTilde]0.*)


Series[IntLimitP,{\[Mu],0,2}]//Normal[#]&//FullSimplify[#,T>0&&m>0&&\[Mu]>0]&


Series[Out[13],{m,0,3}]


FermiInt[T_,m_,\[Mu]_] := T (Integrate[p^2 Log[1 + Exp[\[Mu]/T] Exp[-Sqrt[m^2+p^2]/T]],{p,0,\[Infinity]}] 
+ Integrate[p^2 Log[1 + Exp[-\[Mu]/T] Exp[-Sqrt[m^2+p^2]/T]],{p,0,\[Infinity]}]);
FreeEnergy = -(4 \[Pi]) (2 \[Pi])^(-3) FermiInt


(* ::Section:: *)
(*Evaluating Fermi Integrals for various limits*)


(* ::Text:: *)
(*We're going to define a set of expressions which allow us to evaluate the Fermi integral for various important limits. To start, we define the Fermi Partition function which we need to integrate over to obtain the free energy.*)


FermiPartition[T_,m_,p_,\[Mu]_] := Log[1 + Exp[\[Mu]/T] Exp[-Sqrt[m^2+p^2]/T]];
FreeEnergyCoeff = -(4 \[Pi] T) (2 \[Pi])^(-3);


(* ::Text:: *)
(*We can evaluate the massless limit by integrating over the FermiPartition expression with the mass set to zero m=0 and summing over particle(antiparticle) states \[PlusMinus]\[Mu].*)


FreeEnergyCoeff Total[Refine[ParallelMap[Integrate[p^2 FermiPartition[T,0,p,#],{p,0,\[Infinity]}]&,{-\[Mu],\[Mu]}],T>0&&\[Mu]>0]]
