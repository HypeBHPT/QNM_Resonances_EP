#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"]
(*SetDirectory[NotebookDirectory[]];*)
Print["Starting Notebook with ",ToString[$ProcessorCount], " threads"];
LaunchKernels[$ProcessorCount];


(* ::Section:: *)
(*Routine To Call Notebook with Parameters*)


runWith[aValue_,\[Epsilon]Value_,spinValue_,lValue_, N$dom0Value_,N$dom1Value_, PrecValue_]:=Block[
{ao,\[Epsilon],spin,l, Nz$dom0,Nz$dom1,Prec},
ao=aValue;
\[Epsilon]=\[Epsilon]Value;
spin=spinValue;
l=lValue; 
Nz$dom0=N$dom0Value;
Nz$dom1=N$dom1Value; 
Prec=PrecValue;
<<241002_ElephantFlea_BaseCode.m;
]


(* ::Subsection:: *)
(*Set Parameters*)


ao$min=5;
ao$max=10;
N$a0=5;
ao[i_]:=ao$min + (ao$max-ao$min)*i/N$a0

\[Epsilon]=10^(-2);

spin=2;
l=2;

Nz$dom0=30;
Nz$dom1=30;

Prec=10*MachinePrecision;


(* ::Subsection:: *)
(*Check if differentiation matrices exits*)


fn$D\[Chi]="OperatorMatrix/Dchi_dom0_"<>"N_"<>ToString[Nz$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";
fn$D2\[Chi]="OperatorMatrix/D2chi_dom0_"<>"N_"<>ToString[Nz$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";

check$Dx=FileExistsQ[fn$D\[Chi]];

If[!check$Dx ,
(*Print["Solving for a0 = "<>ToString[N[ao[0]]]];
Time=AbsoluteTiming[runWith[ao[0],\[Epsilon],spin,l, Nz$dom0,Nz$dom1, Prec]]; 
Print["Done in "<>ToString[N@Time[[1]]]<>"s"];*)
runWith[ao[0],\[Epsilon],spin,l, Nz$dom0,Nz$dom1, Prec];
i0=1, 
i0=0;
]


(* ::Subsection:: *)
(*Run over parameter space ao*)


ParallelDo[

(*Print["Solving for a0 = "<>ToString[N[ao[i]]]<>" in thread "<>ToString[$KernelID]<>"("<>ToString[$ProcessorCount]<>")"];
Time=AbsoluteTiming[runWith[ao[i],\[Epsilon],spin,l, Nz$dom0,Nz$dom1, Prec]];
Print["Done in "<>ToString[N@Time[[1]]]<>"s"],*)

runWith[ao[i],\[Epsilon],spin,l, Nz$dom0,Nz$dom1, Prec],

{i,i0,N$a0}
]//AbsoluteTiming


CloseKernels[];
