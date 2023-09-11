(* ::Package:: *)

(* ::Title:: *)
(*AdaptiveEstimator.m*)


(* ::Text:: *)
(*Ruvi Lecamwasam, September 2023*)
(*Wolfram Mathematica script*)
(*Evaluate performance of the CXI-adaptive estimator.*)
(**)
(*Warning: If you change NDim, you must re-start the Mathematica kernel.*)
(*The function \[Psi]PhotonHG is computationally demanding, but called multiple times with the same arguments. To speed up evaluation this function is "memoised", meaning it stores the results each time it is called with new arguments, and next time returns that result without computation. If you change the number of dimensions NDim you need to clear all those saved results, hence have to re-start the Mathematica kernel.*)


(* ::Section:: *)
(*Setup*)


SetDirectory@NotebookDirectory[];
NDim = 40; (* Truncation of the Hermite-Gauss basis *)
(* Re-start the kernel if you change NDim! *)

\[Theta]Range = 5.; (* Range of \[Theta] displacements to consider *)
n\[Theta] = 100; (* \[Theta] increment *)
\[Theta]Vals = Range[-\[Theta]Range,\[Theta]Range,2\[Theta]Range/n\[Theta]];

(* Max iterations when finding optimum measurement angle *)
NMinimizeMethod = "RandomSearch";
NMinimizeMaxIterations = Automatic;

(* Launch as many kernels as CPUs, as otherwise Mathematica sometimes uses less 
LaunchKernels[$ProcessorCount];
Print["Launched "<>ToString@$KernelCount<>" parallel kernels."];*)

PrintLog[string_] := Print[DateString["Time"]<>" "<>string];


(* ::Subsection:: *)
(*Prior distribution*)


\[Sigma]Prior = 0.25; (* Variance of the prior distribution *)
(* Continuous Gaussian as prior distribution *)
priorCts[\[Phi]_] :=  (1/Sqrt[2\[Pi] \[Sigma]Prior^2]) (Exp[-(\[Phi]-1)^2/(2 \[Sigma]Prior^2)]/2+Exp[-(\[Phi]+1)^2/(2 \[Sigma]Prior^2)]/2);

(* Sample the continuous prior to get discrete prior {{\[Phi]1,p1},{\[Phi]2,p2},...} *)
prior = Module[{samples, \[Phi]Vals, pVals, \[Phi]Range},
\[Phi]Range = 2; (* \[Phi] from -\[Phi]Range to +\[Phi]Range *)
samples = Table[{\[Phi],priorCts[\[Phi]]},{\[Phi],-\[Phi]Range,\[Phi]Range,(2\[Phi]Range)/50}]//Chop;
\[Phi]Vals = Transpose[samples][[1]]; (* Discrete parameter values *)
pVals = Transpose[samples][[2]]; (* Discrete parameter probabilities (unnormalised) *)

(* A list {{\[Phi]1,p1},{\[Phi]2,p2},...} *)
{\[Phi]Vals,pVals/Total[pVals]}//Transpose
];

\[Phi]List = Transpose[prior][[1]]; (* {\[Phi]1,\[Phi]2,...} *)
priorp = Transpose[prior][[2]]; (* {p1,p2,...} *)


(* ::Subsection:: *)
(*State in Hermite-Gauss basis*)


(* Wavefunction of a photon centred at \[Phi] *)
\[Psi]Photon[\[Phi]_]:=1/(2\[Pi])^(1/4) Exp[-((x-\[Phi])^2/4)];
(* Wavefunction of qth HG mode, as a function of position x *)
\[Psi]HGq[x_,q_]:=1/(2\[Pi] \[Sigma]hg^2)^(1/4) 1/Sqrt[2^q Factorial[q]] HermiteH[q,x/(Sqrt[2]\[Sigma]hg)]Exp[-(x^2/(4\[Sigma]hg^2))]/.\[Sigma]hg->2;

(* State of photon centered at \[Phi], in the \[Theta]-displaced HG basis *)
(*
We want to compute the overlap:
\[Psi]PhotonHGOld[\[Phi]_,\[Theta]_] := Table[
Quiet@NIntegrate[Conjugate[\[Psi]HGq[x-\[Theta],q]]\[Psi]Photon[\[Phi]],{x,-\[Infinity],\[Infinity]},Method->{Automatic,"SymbolicProcessing"->0}]//Chop
,{q,0,NDim-1}]

This integral has an analytic form which we can find using:
$Assumptions = {\[Theta]\[Element]Reals,\[Phi]\[Element]Reals,n>=0};
overlapIntegralFunction = Integrate[x^n\[ExponentialE]^(-(1/16) (x-\[Theta])^2-1/4 (x-\[Phi])^2),{x,-\[Infinity],\[Infinity]}]
*)

overlapIntegralFunction[n_,\[Theta]_,\[Phi]_] := 5^(-1-n/2) E^(1/16 (-\[Theta]^2-4 \[Phi]^2)) (2^(1+2 n) Sqrt[5] (1+(-1)^n) Gamma[(1+n)/2] Hypergeometric1F1[(1+n)/2,1/2,1/80 (\[Theta]+4 \[Phi])^2]-((-4)^n-4^n) (\[Theta]+4 \[Phi]) Gamma[1+n/2] Hypergeometric1F1[(2+n)/2,3/2,1/80 (\[Theta]+4 \[Phi])^2]);

Overlap[q_,\[Phi]_,\[Theta]_] := Module[{xnCoeffs},
	xnCoeffs = HermiteH[q,(x-\[Theta])/(2Sqrt[2])]//CoefficientList[#,x]&;
	1/(2.Sqrt[\[Pi]] Sqrt[2^q Factorial[q]]) Plus@@Table[xnCoeffs[[n+1]]overlapIntegralFunction[n,\[Theta],\[Phi]],{n,0,q}]
]

(* State of a photon centred at \[Theta], measured in a \[Phi]-displaced HG basis *)
(* Memoized, since \[Psi]PhotonHG is intensive *)
\[Psi]PhotonHG[\[Phi]_,\[Theta]_]:=\[Psi]PhotonHG[\[Phi],\[Theta]]=Table[Overlap[q,\[Phi],\[Theta]],{q,0,NDim-1}]


(* Density matrix of photon centred at \[Phi], in the \[Theta]-dispaced HG basis *)
\[Rho]PhotonHG[\[Phi]_,\[Theta]_] := Module[{\[Psi]hg},
\[Psi]hg=\[Psi]PhotonHG[\[Phi],\[Theta]];
KroneckerProduct[\[Psi]hg,Conjugate[\[Psi]hg]]
]

(* HG basis density matrix averaged over prior distribution for \[Phi] *)
\[Rho]\[Theta]HG[\[Theta]_]:=Sum[priorp[[i]]\[Rho]PhotonHG[\[Phi]List[[i]],\[Theta]],{i,1,Length[\[Phi]List]}]//Chop;


(* ::Subsection:: *)
(*Measurement and parameter distributions*)


(* Parameter distribution p\[CapitalPhi][\[Phi]i]=pi *)
p\[CapitalPhi][\[Phi]_] := Nearest[\[Phi]List->priorp,\[Phi]]//First;

(* Conditional probability of single measurement outcome given \[Phi] *)
pMl\[CapitalPhi][n_,\[Phi]_,\[Theta]_] := \[Rho]PhotonHG[\[Phi],\[Theta]][[n+1,n+1]];

(* Conditional probability of multiple measurement outcomes given \[Phi] *)
pMl\[CapitalPhi][nList_List,\[Phi]_,\[Theta]List_List] := Module[{measurements},
	measurements = Thread[{nList,\[Theta]List}];
	Times@@Table[pMl\[CapitalPhi][m[[1]],\[Phi],m[[2]]],{m,measurements}]
]

(* Joint probability distribution *)
pM\[CapitalPhi][nList_List,\[Phi]_,\[Theta]List_List] := pMl\[CapitalPhi][nList,\[Phi],\[Theta]List]p\[CapitalPhi][\[Phi]]

(* Measurement probability distribution *)
pM[nList_List,\[Theta]List_List] := Sum[pM\[CapitalPhi][nList,\[Phi],\[Theta]List],{\[Phi],\[Phi]List}]

(* Posterior distribution *)
p\[CapitalPhi]lM[\[Phi]_,nList_List,\[Theta]List_] := pM\[CapitalPhi][nList,\[Phi],\[Theta]List]/pM[nList,\[Theta]List];

(* Convert a list of probabilities to distribution function *) 
DistributionFromPlist[probs_]:=Module[{pDist},
	pDist[\[Phi]_] := Nearest[\[Phi]List->probs,\[Phi]]//First;
	pDist
]


(* ::Subsection:: *)
(*Coherence*)


(* Von Neumann entropy *)
S[\[Rho]_] := Module[{evals},
	(* Select non-zero eigenvalues, this will be empty for a pure state *)
	evals = Select[Chop@Eigenvalues[\[Rho]],#>0&];
	(* Return zero entropy for a pure state *)
	If[Length@evals>0,
		-# Log[#]&/@evals//Total,
	0]
]

(* Decohere a matrix in the HG basis *)
\[CapitalDelta]HG[\[Rho]_] := DiagonalMatrix@Diagonal@\[Rho];

(* Find the coherence in the HG basis *)
CHG[\[Rho]_] := S[\[CapitalDelta]HG[\[Rho]]]-S[\[Rho]];
(* Memoized version *)
CHG[\[Phi]_,\[Theta]_] := CHG[\[Phi],\[Theta]] = CHG[\[Rho]PhotonHG[\[Phi],\[Theta]]];

(* Coherence of encoding with respect to a probabilty distribution *)
Ce[\[Theta]_,probs_] := Module[{Ci,Cf,\[Rho]pDist,pDist},
	pDist = DistributionFromPlist[probs];
	(* Initial coherence *)
	Ci = Sum[CHG[\[Phi],\[Theta]] pDist[\[Phi]],{\[Phi],\[Phi]List}] // Chop;
	(* Final coherence *)
	\[Rho]pDist = Sum[\[Rho]PhotonHG[\[Phi],\[Theta]]pDist[\[Phi]],{\[Phi],\[Phi]List}] // Chop;
	Cf = CHG[\[Rho]pDist];
	(* Ce is the difference *)
	Ci-Cf
]


(* ::Subsection:: *)
(*Measurement functions*)


(* Initial distribution *)
(* {{probability distribution},{measurement sequence}} *)
zeroMeasurements = {priorp,{}};

(* Measure at constant displacement *)
MeasureAtConst\[Theta][probs_,measurements_,\[Theta]0_] := Module[
	(* XPre \[Rule] X prior to measurement *)
	{p\[CapitalPhi]Pre,pM\[CapitalPhi]Pre,pMPre,p\[CapitalPhi]lMPre,q,newProbs,newMeasurementList},

	(* Turn list of probabilities into a probability distribution function *)
	p\[CapitalPhi]Pre = DistributionFromPlist[probs];

	(* Construct other probability distributions *)
	pM\[CapitalPhi]Pre[n_,\[Phi]_] := pMl\[CapitalPhi][{n},\[Phi],{\[Theta]0}]p\[CapitalPhi]Pre[\[Phi]];
	pMPre[n_] := Sum[pM\[CapitalPhi]Pre[n,\[Phi]],{\[Phi],\[Phi]List}];
	p\[CapitalPhi]lMPre[\[Phi]_,n_] := pM\[CapitalPhi]Pre[n,\[Phi]]/pMPre[n];
	
	(* Randomly simulate a measurement q with probabilities pMPre[q] *)
	q = RandomChoice[(pMPre/@Range[0,NDim-1])->Range[0,NDim-1]];
	
	newProbs = p\[CapitalPhi]lMPre[#,q]&/@\[Phi]List; (* p\[CapitalPhi]Post[\[Phi]] = p\[CapitalPhi]lMPre[\[Phi],q] *)
	newMeasurementList = Append[measurements,{q,\[Theta]0}];
	
	{newProbs,newMeasurementList}
]


(* Use coherence to find the best measurement angle *)
Options[FindBestMeasurement]={parallel->False,maxIterations->NMinimizeMaxIterations,method->NMinimizeMethod};
FindBestMeasurement[probs_,OptionsPattern[]] := Module[{ceList,ceFun,table},
	
	table = If[OptionValue@parallel,ParallelTable,Table];
	
	ceList = table[{\[Theta],Ce[\[Theta],probs]},{\[Theta],\[Theta]Vals}];
	ceFun = Interpolation@ceList;
	\[Theta]/.Last@NMinimize[{ceFun[\[Theta]],-\[Theta]Range<\[Theta]<\[Theta]Range},\[Theta],MaxIterations->OptionValue@maxIterations,Method->OptionValue@method]
]

MeasureAdaptively[probs_,measurements_] := Module[
	{\[Theta]Best,p\[CapitalPhi]Pre,pM\[CapitalPhi]Pre,pMPre,p\[CapitalPhi]lMPre,q,newProbs,newMeasurementList},
	(* Turn list of probabilities into a probability distribution function *)
	p\[CapitalPhi]Pre = DistributionFromPlist[probs];

	(* Find best measurement angle *)
	\[Theta]Best = FindBestMeasurement[probs];

	(* Construct other probability distributions*)
	pM\[CapitalPhi]Pre[n_,\[Phi]_] := pMl\[CapitalPhi][{n},\[Phi],{\[Theta]Best}]p\[CapitalPhi]Pre[\[Phi]];
	pMPre[n_] := Sum[pM\[CapitalPhi]Pre[n,\[Phi]],{\[Phi],\[Phi]List}];
	p\[CapitalPhi]lMPre[\[Phi]_,n_] := pM\[CapitalPhi]Pre[n,\[Phi]]/pMPre[n];
	
	(* Randomly simulate a measurement q with probabilities pMPre[q] *)
	q = RandomChoice[(pMPre/@Range[0,NDim-1])->Range[0,NDim-1]];

	newProbs = p\[CapitalPhi]lMPre[#,q]&/@\[Phi]List; (* p\[CapitalPhi]Post[\[Phi]] = p\[CapitalPhi]lMPre[\[Phi],q] *)
	newMeasurementList = Append[measurements,{q,\[Theta]Best}];
	
	{newProbs,newMeasurementList}
]


(* ::Section:: *)
(*Simulate measurements*)


nMeas = 10; (* Number of measurements to make *)
nCores = $KernelCount; (* Number of cores to simulate *)
Print["Number of cores: "<>ToString@nCores];
nTrajPerCore = 10; (* Number of trajectories to simulate per core *)
Print["Trajectories per core: "<>ToString@nTrajPerCore];

Print["Simulating "<>ToString@nMeas<>" measurements, "<>ToString[nCores*nTrajPerCore]<>" trajectories."];


\[Theta]opt = FindBestMeasurement[zeroMeasurements[[1]],parallel->True];


(* Check truncation NDim is appropriate by summing pM[q] *)
Module[{probsum},
		probsum = ParallelTable[pM[{q},{0}],{q,0,NDim-1}]//Flatten//Total;
	Print["For NDim = "<>ToString@NDim<>", sum of measurement proabilities pM[q] is "<>ToString@probsum]
]


(* ::Subsection:: *)
(*Constant \[Phi]*)


MeasureAtZero[probs_,measurements_] := MeasureAtConst\[Theta][probs,measurements,0];
MeasureAtZeroPointFive[probs_,measurements_] := MeasureAtConst\[Theta][probs,measurements,0.5];
MeasureAtOpt[probs_,measurements_]:=MeasureAtConst\[Theta][probs,measurements,\[Theta]opt];


Print["Measuring \[Theta]=0"];
resultsZero = ParallelTable[
		Table[NestList[MeasureAtZero@@#&,zeroMeasurements,nMeas],{i,1,nTrajPerCore}],{j,1,nCores}
]//Flatten[#,1]&//Transpose;
	
Print["Measuring \[Theta]=0.5"];
resultsZeroPointFive = ParallelTable[
		Table[NestList[MeasureAtZeroPointFive@@#&,zeroMeasurements,nMeas],{i,1,nTrajPerCore}],{j,1,nCores}
]//Flatten[#,1]&//Transpose;

Print["Measuring \[Theta]=\[Theta]opt"]
resultsOpt = ParallelTable[
		Table[NestList[MeasureAtOpt@@#&,zeroMeasurements,nMeas],{i,1,nTrajPerCore}],{j,1,nCores}
]//Flatten[#,1]&//Transpose;


(* ::Subsection:: *)
(*Adaptive*)


Print["Measuring adaptively"];
resultsAdaptive = ParallelTable[
	Table[NestList[MeasureAdaptively@@#&,zeroMeasurements,nMeas],{i,1,nTrajPerCore}],{j,1,nCores}
]//Flatten[#,1]&//Transpose;


(* ::Section:: *)
(*Analyse performance*)


(* ::Subsection:: *)
(*Calculate error*)


(* Squared error as our error function *)
Error[\[CapitalDelta]_] := \[CapitalDelta]^2;

(* Our estimate is the mean of the distribution *)
EstimateFromDist[pDist_] := Sum[\[Phi] pDist[\[Phi]],{\[Phi],\[Phi]List}]//Chop;
EstimateFromPList[probs_] := EstimateFromDist[DistributionFromPlist[probs]]

(* Error in one trajectory *)
ErrorInOneTrajectory[results_] := Module[
	{probs,trajectory\[CapitalPhi]lM,measurements,\[Phi]estimate,performance},
	probs = results[[1]];
	measurements = results[[2]];
	(* Distribution for \[CapitalPhi] conditioned on measurement results M *)
	trajectory\[CapitalPhi]lM = DistributionFromPlist[probs];
	
	(* Estimate of \[Phi] from measurement results *)
	\[Phi]estimate = EstimateFromPList[probs];
	(* \!\(
\*SubscriptBox[\(\[Sum]\), \(\[Phi]\)]\(
\(\*SubscriptBox[\(p\), \(\[CapitalPhi] | M\)]\)[\[Phi]]
\*SuperscriptBox[\((\[Phi] - 
\*OverscriptBox[\(\[Phi]\), \(_\)])\), \(2\)]\)\) *)
	(* In MeasurementPerformance[] we take the mean of this, which corresponds
		to dividing by the number of measurements *)
	performance = Sum[trajectory\[CapitalPhi]lM[\[Phi]] Error[\[Phi]-\[Phi]estimate],{\[Phi],\[Phi]List}];
	
	{performance,measurements}
]

(* Map ErrorInOneTrajectory[] over all trajectories *)
ErrorInAllTrajectories[results_] := ParallelMap[
	(* Divide all results into nCores blocks, one block for each CPU *)
	(* Gives each core one 'big' job, rather than repeatedly giving it small ones *)
	ErrorInOneTrajectory/@#&,Partition[Flatten[results,1],{nCores}]
	]//Flatten[#,1]&//GatherBy[#,Length@Last@#&]&;
	(* Flatten \[Rule] Combine the nCores blocks black into one list *)
	(* GatherBy \[Rule] Group the results by the number of measurements *)
	
(* {{n,Error after n measurements}} *)
(* Mean[] corresponds to dividing by number of measurement results *)
MeasurementPerformance[results_] := results//ErrorInAllTrajectories\
	//Table[{i-1,Mean[First/@#[[i]]]},{i,1,Length@#}]&;
	
(* Variance *)
MeasurementPerformanceVariance[results_] := results//ErrorInAllTrajectories\
	//Table[{i-1,Variance[First/@#[[i]]]},{i,1,Length@#}]&;


Print["Analysing results"];
performance=<||>;

performance["ntraj"] = nCores * nTrajPerCore;
performance["nDim"] = NDim;
performance["\[Theta]Range"] = {\[Theta]Range,\[CapitalDelta]\[Theta]};

performance["0"] = resultsZero//MeasurementPerformance;
performance["0.5"] = resultsZeroPointFive//MeasurementPerformance;
performance["Opt"] = resultsOpt//MeasurementPerformance;
performance["Adaptive"] = resultsAdaptive//MeasurementPerformance;

performance["0Error"] = resultsZero//MeasurementPerformanceVariance;
performance["0.5Error"] = resultsZeroPointFive//MeasurementPerformanceVariance;
performance["OptError"] = resultsOpt//MeasurementPerformanceVariance;
performance["AdaptiveError"] = resultsAdaptive//MeasurementPerformanceVariance;

(* Print the results *)
Table[Print["Performance for "<>key<>":"];Print[performance[key]],{key,{"0","0.5","Opt","Adaptive"}}];


Print["Exporting results"];
Export["performance.wdx",performance];
