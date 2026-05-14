#!/usr/bin/env wolframscript
(* ::Package:: *)

Clear["Global`*"]
SetDirectory[NotebookDirectory[]];
Print["Starting Notebook with ",ToString[$ProcessorCount], " threads"];
LaunchKernels[$ProcessorCount];


(* ::Section:: *)
(*Routine To Call Notebook with Parameters*)


runWith[aValue_,\[Epsilon]Value_,spinValue_,lValue_, N$dom0Value_,N$dom1Value_,  NRef$dom0Value_,NRef$dom1Value_, \[Epsilon]Ref$Value_, PrecValue_]:=Block[
{ao,\[Epsilon],spin,l, Nz$dom0,Nz$dom1, NzRef$dom0,NzRef$dom1, \[Epsilon]QNM, Prec},
ao=aValue;
\[Epsilon]=\[Epsilon]Value;
spin=spinValue;
l=lValue; 
Nz$dom0=N$dom0Value;
Nz$dom1=N$dom1Value; 
NzRef$dom0=NRef$dom0Value;
NzRef$dom1=NRef$dom1Value; 
\[Epsilon]QNM=\[Epsilon]Ref$Value;
Prec=PrecValue;
<<"250319_Cluster_ElephantFlea_QNM_EV_EF.m";
]


(* ::Subsection:: *)
(*Set Parameters*)


ao$min=1154/100;
ao$max=116/10;
N$a0=100;
ao[i_]:=ao$min + (ao$max-ao$min)*i/N$a0

\[Epsilon]=204/100000;

spin=2;
l=2;

Nz$dom0=200;
Nz$dom1=200;

NzRef$dom0=150;
NzRef$dom1=150;

\[Epsilon]Ref=10^(-4);

Prec=10*MachinePrecision;


(* ::Subsection:: *)
(*Check if differentiation matrices exits*)


fn$D\[Chi]="OperatorMatrix/Dchi_dom0_"<>"N_"<>ToString[Nz$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";
fn$D2\[Chi]="OperatorMatrix/D2chi_dom0_"<>"N_"<>ToString[Nz$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";

check$Dx=FileExistsQ[fn$D\[Chi]];

fn$D\[Chi]Ref="OperatorMatrix/Dchi_dom0_"<>"N_"<>ToString[NzRef$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";
fn$D2\[Chi]Ref="OperatorMatrix/D2chi_dom0_"<>"N_"<>ToString[NzRef$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";

check$DxRef=FileExistsQ[fn$D\[Chi]Ref];

If[(!check$Dx || !check$DxRef) ,
(*Print["Solving for a0 = "<>ToString[N[ao[0]]]];
Time=AbsoluteTiming[runWith[ao[0],\[Epsilon],spin,l, Nz$dom0,Nz$dom1, Prec]]; 
Print["Done in "<>ToString[N@Time[[1]]]<>"s"];*)
runWith[ao[0],\[Epsilon],spin,l, Nz$dom0,Nz$dom1, NzRef$dom0,NzRef$dom1, \[Epsilon]Ref,Prec];
i0=1, 
i0=0;
]


(* ::Subsection:: *)
(*Run over parameter space ao*)


ParallelDo[

(*Print["Solving for a0 = "<>ToString[N[ao[i]]]<>" in thread "<>ToString[$KernelID]<>"("<>ToString[$ProcessorCount]<>")"];
Time=AbsoluteTiming[runWith[ao[i],\[Epsilon],spin,l, Nz$dom0,Nz$dom1, Prec]];
Print["Done in "<>ToString[N@Time[[1]]]<>"s"],*)
runWith[ao[i],\[Epsilon],spin,l, Nz$dom0,Nz$dom1, NzRef$dom0,NzRef$dom1, \[Epsilon]Ref,Prec],

{i,i0,N$a0}
]//AbsoluteTiming


CloseKernels[];
