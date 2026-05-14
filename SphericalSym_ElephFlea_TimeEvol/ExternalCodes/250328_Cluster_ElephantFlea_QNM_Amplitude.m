(* ::Package:: *)

(* ::Title:: *)
(*QNM Amplitude: Time and Frequency domain analysis*)


(*Clear["Global`*"]
<<MaTeX`
SetDirectory[NotebookDirectory[]];
$Assumptions={M>0}
DirName=NotebookDirectory[];*)


(* ::Chapter::Closed:: *)
(* Parameters*)


(* ::Section:: *)
(*Physical Parameters*)


(* ::Subsection:: *)
(*Wave equation Parameters*)


Mbh=rh/2;
\[Lambda]=rh;
(*spin=2;
l=2;*)


(* ::Subsection:: *)
(*Bump Parameters*)


 ModelName="BumpPT";
(* ao=1;(* 115722/10000;*) (*120000/10000;*)
 \[Epsilon]=204/100000;*)


(* ::Section:: *)
(*Numerical Parameters*)


(*Prec=10*MachinePrecision;*)


(*Nz$dom0=75;
Nz$dom1=75;*)


(*Nz$QNMload$dom0=100;
Nz$QNMload$dom1=100;*)


(* ::Subsection:: *)
(*QNM Selection *)


(*N$QNM$max=3;*)


(* ::Chapter::Closed:: *)
(*Schwarzschild Spacetime*)


(* ::Section:: *)
(*Metric Function*)


a=1-rh/r;
b=a;
f=FullSimplify[Sqrt[a*b], Assumptions->{rh>0, r>rh}];
\[Zeta]$of$r=FullSimplify[Sqrt[a/b], Assumptions->{rh>0, r>rh}];


m=r/2*(1-b);

fAsymp=Normal@Series[f,{r,Infinity,1}]//Expand;
M$ADM=-Coefficient[fAsymp,r,-1]/2//FullSimplify;
SubsM={M->M$ADM}


(* ::Chapter::Closed:: *)
(*Hyperboloidal Formulation*)


(* ::Section::Closed:: *)
(*Radial Compactfication*)


(* ::Subsection:: *)
(*Areal Radius*)


\[Rho]1=0;(*rh/\[Lambda] 1/(1-\[Sigma]0)*)
\[Rho]0=rh/\[Lambda]-\[Rho]1;

\[Rho][\[Sigma]_]:=\[Rho]0+\[Rho]1*\[Sigma]
\[Beta]=\[Rho][\[Sigma]]-\[Sigma]*D[\[Rho][\[Sigma]],\[Sigma]]//FullSimplify;


(* ::Subsection:: *)
(*Radial Mapping*)


r$of$\[Sigma][\[Sigma]_]:=\[Lambda] \[Rho][\[Sigma]]/\[Sigma];
dr$d\[Sigma]=D[r$of$\[Sigma][\[Sigma]],\[Sigma]];

\[Sigma]$of$r=\[Sigma]/.Solve[r$of$\[Sigma][\[Sigma]]==r,\[Sigma]][[1]];


\[Sigma]h=\[Sigma]hrz/.Solve[r$of$\[Sigma][\[Sigma]hrz]==rh,\[Sigma]hrz][[1]];


(* ::Subsection:: *)
(*Metric function and deviation function in compact coordinate*)


F=Simplify[f/.r->r$of$\[Sigma][\[Sigma]]];
\[Zeta]=Simplify[\[Zeta]$of$r/.r->r$of$\[Sigma][\[Sigma]]];


(* ::Subsubsection:: *)
(*Find Horizons*)


Subs\[Sigma]Hrz=Solve[{F==0, \[Sigma]>0},\[Sigma],Reals];
\[Sigma]Hrz=\[Sigma]/.Subs\[Sigma]Hrz;
nHrz=Length@\[Sigma]Hrz;
rHrz=r$of$\[Sigma][\[Sigma]Hrz];


(* ::Subsubsection:: *)
(*Representation as in eq. 2.13 - Event Horizon*)


root$rh=(1-rHrz/r)/.r->r$of$\[Sigma][\[Sigma]]//FullSimplify;
Kh=Expand[F/root$rh]//FullSimplify;


(* ::Section::Closed:: *)
(*Dimensionless Tortoise Coordinate*)


(* ::Subsection:: *)
(*Differential Equation (eq. 2.11)*)


dx$d\[Sigma]=FullSimplify[-\[Beta]/(\[Sigma]^2F)];


(* ::Subsection:: *)
(*Contribution from Infinity (eq. 4.7) *)


x$0=\[Rho]0/\[Sigma]-2M/\[Lambda] Log[\[Sigma]];
dx0$d\[Sigma]=FullSimplify[D[x$0,\[Sigma]]];


(* ::Subsection:: *)
(*Contribution from Event Horizon (eq. 4.5)*)


KBHHrz=Table[Limit[Kh[[i]],Subs\[Sigma]Hrz[[i]]],{i,1,nHrz}];
x$BHHrz=Table[rHrz[[i]]/(\[Lambda]*KBHHrz[[i]])*Log[\[Sigma]Hrz[[i]]-\[Sigma]],{i,1,nHrz}];

dx$BHHrz$d\[Sigma]=Table[D[x$BHHrz[[i]],\[Sigma]],{i,1,nHrz}];


(* ::Subsection:: *)
(*Regular Contribution*)


dx$d\[Sigma]$Reg=dx$d\[Sigma]-dx0$d\[Sigma] - Sum[dx$BHHrz$d\[Sigma][[i]],{i,1,nHrz}]/.SubsM;
Subs$dxd\[Sigma]$Reg=x$Reg'[\[Sigma]]->dx$d\[Sigma]$Reg;


(* ::Subsubsection:: *)
(*Check Tortoise*)


x$Reg[\[Sigma]_]:=0;
x=x$0+ Sum[x$BHHrz[[i]],{i,1,nHrz}]+x$Reg[\[Sigma]]/.SubsM//Simplify;
dx$d\[Sigma]$Paper=D[x,\[Sigma]]/.Subs$dxd\[Sigma]$Reg;
FullSimplify[dx$d\[Sigma]$Paper==dx$d\[Sigma]]


(* ::Section:: *)
(*Height Function*)


(* ::Subsection:: *)
(*Out-In Strategy*)


H=-x$0 +Sum[x$BHHrz[[i]],{i,1,nHrz}]-x$Reg[\[Sigma]]/.SubsM//Simplify;
dH$d\[Sigma]=D[H,\[Sigma]]/.Subs$dxd\[Sigma]$Reg;


(* ::Subsection:: *)
(*Metric Functions*)


p=FullSimplify[-1/dx$d\[Sigma]]
(*Plot[p/.\[Lambda]->rh//Simplify,{\[Sigma],0,\[Sigma]Hrz[[1]]}, WorkingPrecision->100, PlotRange->All]*)


\[Gamma]=dH$d\[Sigma]*p//FullSimplify
(*Plot[\[Gamma]/.\[Lambda]->rh//Simplify,{\[Sigma],0,\[Sigma]Hrz[[1]]}, WorkingPrecision->100, PlotRange->All]*)


(*Series[\[Gamma],{\[Sigma],\[Sigma]Hrz[[1]],0}]*)
w=(1-\[Gamma]^2)/p//FullSimplify
(*Plot[w/.\[Lambda]->rh//Simplify,{\[Sigma],0,\[Sigma]Hrz[[1]]}, WorkingPrecision->100, PlotRange->{0,All}]*)


(* ::Chapter::Closed:: *)
(*Wave Equation*)


(* ::Section:: *)
(*Regge Wheeler Potential*)


If[spin==0,VRW=a/r^2*(l*(l+1) + r*D[a*b,r]/(2*a))];
If[Abs@spin==2, VRW=a/r^2* ( l*(l+1) - 6m/r + D[m,r])];

qRW=\[Lambda]^2/p * VRW/.r->r$of$\[Sigma][\[Sigma]]//FullSimplify


(* ::Subsection:: *)
(*Bump*)


Vbump =  \[Lambda]^(-2) Sech[x-ao]^2/.r->r$of$\[Sigma][\[Sigma]]//FullSimplify;
(* \[Lambda]^(-2) Exp[ -(x-ao)^2/(2 \[Delta])  ]*)

qBump=\[Lambda]^2/p * Vbump;


If[\[Epsilon]==0,
\[Sigma]Bump=1/2,
\[Sigma]Bump=\[Sigma]/.N[Solve[x-ao==0,\[Sigma]],50*MachinePrecision][[1]]
];


(* ::Subsection:: *)
(*Final Potential*)


V = VRW + \[Epsilon] Vbump/.r->r$of$\[Sigma][\[Sigma]];
q = \[Lambda]^2/p * V//FullSimplify;


(* ::Section::Closed:: *)
(*Functions in Wave Equation*)


W=w;

\[Alpha]2=p;
\[Alpha]1=D[p,\[Sigma]];
\[Alpha]0=-q;

\[Beta]1=2 \[Gamma];
\[Beta]0=D[\[Gamma],\[Sigma]];


(* ::Chapter:: *)
(*QNM Eigenfunctions*)


(* ::Section:: *)
(*Discrete Functions and Matrices*)


(* ::Subsection:: *)
(*AnMR*)


x$of$\[Chi][\[Chi]_,\[Kappa]_,xb_]:= If[ \[Kappa]==0,
\[Chi]
,
xb(1- 2 Sinh[ \[Kappa] (1-xb \[Chi])]/Sinh[2\[Kappa]]    )
]
\[Chi]$of$x[X_,\[Kappa]_,xb_]:= If[ \[Kappa]==0,
X
,
xb*(1-ArcSinh[(1 - X*xb)*Sinh[2\[Kappa]]/2]/\[Kappa])
]


(* ::Subsection:: *)
(*Spectral Routines*)


(* ::Subsubsection::Closed:: *)
(*Chebyshev Coefficients*)


ChebGaussCoefficients[f_, Prec_]:=Module[
{nf, Nf,c,x},
nf=Length@f;
Nf=nf-1;

x[j_]:=N[ Cos[Pi*(j+1/2)/(Nf+1)] ,Prec];

c=Table[  2 Sum[f[[j+1]]*ChebyshevT[m,x[j]],{j,0,Nf}]  /nf        ,{m,0,Nf}]
];

ChebLobattoCoefficients[f_, Prec_]:=Module[
{nf, Nf,c,x},
nf=Length@f;
Nf=nf-1;

x[j_]:=N[Cos[j \[Pi]/Nf],Prec];

c=Table[  ( 2-KroneckerDelta[m,Nf]) /(2Nf) (  f[[1]]+(-1)^m f[[nf]] + 2 Sum[f[[j+1]]*ChebyshevT[m,x[j]],{j,1,Nf-1}]  )        ,{m,0,Nf}]
];


(* ::Subsubsection::Closed:: *)
(*Spectral Differentiation matrices*)


DiffMatrix[a_,b_,Nz_,Prec_]:=Module[
{k,i,j,dx, \[CapitalDelta],Dy},

k[i_]:=If[i*(i-Nz)==0, 2,1];
dx[i_,j_]:=If[i!=j,
k[i]*(-1)^(i-j)/(k[j]*(x[i,Nz]-x[j,Nz])),
0
];
\[CapitalDelta]=b-a;
Dy=N[2/\[CapitalDelta] * Table[dx[i,j], {i,0,Nz}, {j,0,Nz}],Prec ];

For[i=0, i<=Nz, i++,
Dy[[i+1,i+1]]= -Sum[ Dy[[i+1,j+1]], {j,0,Nz} ];
];
Dy
]


(* ::Subsubsection:: *)
(*Chebyshev Interpolation*)


ChebInterpolation[c_,a_,b_,x_, \[Kappa]_, xb_, Prec_]:=Module[
{nc, Nc, \[Chi], y},
nc=Length@c;
Nc=nc-1;

y=(2 x-(a+b))/(b-a);
\[Chi]=\[Chi]$of$x[y,\[Kappa],xb];

c[[1]]/2+Sum[c[[i+1]]*ChebyshevT[i,\[Chi]],{i,1,Nc}]
];

ChebyLobattoInterpolationVector[xl_,Nl_,x\[Alpha]_,Prec_]:=Module[
{\[Alpha],l},
Table[
If[l==0,
1/(2Nl)+Sum[(2-KroneckerDelta[j,Nl])ChebyshevT[j,x\[Alpha]]/(2Nl),{j,1,Nl}],
If[l==Nl,
1/(2Nl)+Sum[(2-KroneckerDelta[j,Nl])(-1)^j ChebyshevT[j,x\[Alpha]]/(2Nl),{j,1,Nl}],
1/(Nl)+Sum[(2-KroneckerDelta[j,Nl]) ChebyshevT[j,xl[[l+1]]]ChebyshevT[j,x\[Alpha]]/(Nl),{j,1,Nl}]
]
],
{l,0,Nl}
]

];

ChebyLobattoInterpolationMatrix[xl_,Nl_,x\[Alpha]_,N\[Alpha]_,Prec_]:=Module[
{\[Alpha],l},
Table[
If[l==0,
1/(2Nl)+Sum[(2-KroneckerDelta[j,Nl])ChebyshevT[j,x\[Alpha][[\[Alpha]+1]]]/(2Nl),{j,1,Nl}],
If[l==Nl,
1/(2Nl)+Sum[(2-KroneckerDelta[j,Nl])(-1)^j ChebyshevT[j,x\[Alpha][[\[Alpha]+1]]]/(2Nl),{j,1,Nl}],
1/(Nl)+Sum[(2-KroneckerDelta[j,Nl]) ChebyshevT[j,xl[[l+1]]]ChebyshevT[j,x\[Alpha][[\[Alpha]+1]]]/(Nl),{j,1,Nl}]
]
],
{\[Alpha],0,N\[Alpha]},{l,0,Nl}
]

]


(* ::Subsubsection::Closed:: *)
(*Chebyshev Integration*)


DefiniteIntegral[c_,a_,b_]:=Module[
{kf,nc,Nc, Int},
nc=Length@c;
Nc=nc-1;
kf=Floor[Nc/2];
Int=c[[1]]- 2Sum[c[[2 k +1]]/(4 k ^2 -1),{k,1,kf}];
(b-a)*Int/2
];

ChebyLobattoHilbertIntegration[z_,Nz_,Prec_]:=Module[
{i,j},
Table[
If[i!=j,
0,
If[i==0 || i==Nz,
(1/Nz - Sum[(2-KroneckerDelta[2k,Nz])/(Nz(4k^2-1)),{k,1,Floor[Nz/2]}]),
(2/Nz - 2Sum[(2-KroneckerDelta[2k,Nz])*ChebyshevT[2k,z[[i+1]]]/(Nz(4k^2-1)),{k,1,Floor[Nz/2]}])
]
],
 {i,0,Nz}, {j,0,Nz}
]
]


(* ::Subsection:: *)
(*Numerical parameters*)


x=.

zi=0;
zp=\[Sigma]Bump;
zf=1;

\[CapitalDelta]z$dom0=zp-zi;
\[CapitalDelta]z$dom1=zf-zp;


nz$dom0=Nz$dom0+1; (*Total number of points*)
nz$dom1=Nz$dom1+1; (*Total number of points*)

x[i_,Nz_]:=Cos[(\[Pi] i)/Nz];

\[Chi]$dom0=N[Table[x[i,Nz$dom0],{i,0,Nz$dom0}],Prec];
\[Chi]$dom1=N[Table[x[i,Nz$dom1],{i,0,Nz$dom1}],Prec];

ZeroArray$dom0=ConstantArray[0,nz$dom0];
ZeroArray$dom1=ConstantArray[0,nz$dom1];


(* ::Subsubsection:: *)
(*AnMR Parameter*)


(* ::Text:: *)
(*DOMAIN 0 \[Sigma] \[Element] [0, \[Sigma]bump]*)


xb$dom0=1;

\[Kappa]$dom0=If[ao<10, 0,
If[10<=ao<=35,77/100 + 14/100 Exp[6*ao/100],2]
];

If[\[Kappa]$dom0==0,
X$dom0=\[Chi]$dom0;
XRef$dom0=\[Chi]Ref$dom0,

X$dom0=N[xb$dom0*(1-2*Sinh[\[Kappa]$dom0*(1-xb$dom0*\[Chi]$dom0)]/Sinh[2*\[Kappa]$dom0]),Prec];
XRef$dom0=N[xb$dom0*(1-2*Sinh[\[Kappa]$dom0*(1-xb$dom0*\[Chi]Ref$dom0)]/Sinh[2*\[Kappa]$dom0]),Prec];
];
z$dom0=N[zi+1/2 \[CapitalDelta]z$dom0*(1+X$dom0),Prec];
zRef$dom0=N[zi+1/2 \[CapitalDelta]z$dom0*(1+XRef$dom0),Prec];


(* ::Text:: *)
(*DOMAIN 1 \[Sigma] \[Element] [\[Sigma]bump,1]*)


xb$dom1=-1;
\[Kappa]$dom1=125/100 * ao^(3/10);
If[\[Kappa]$dom1==0,
X$dom1=\[Chi]$dom1;
XRef$dom1=\[Chi]Ref$dom1
,
X$dom1=N[xb$dom1*(1-2*Sinh[\[Kappa]$dom1*(1-xb$dom1*\[Chi]$dom1)]/Sinh[2*\[Kappa]$dom1]),Prec];
XRef$dom1=N[xb$dom1*(1-2*Sinh[\[Kappa]$dom1*(1-xb$dom1*\[Chi]Ref$dom1)]/Sinh[2*\[Kappa]$dom1]),Prec];
];
z$dom1=N[zp+1/2 \[CapitalDelta]z$dom1*(1+X$dom1),Prec];
zRef$dom1=N[zp+1/2 \[CapitalDelta]z$dom1*(1+XRef$dom1),Prec];


(* ::Subsection::Closed:: *)
(*Differentiation Matrices*)


(* ::Subsubsection:: *)
(*Create directory to store differentiation matrices*)


DirName$Matrices="OperatorMatrix";
If[!DirectoryQ[DirName$Matrices],CreateDirectory[DirName$Matrices]]


(* ::Subsubsection:: *)
(*Create or Load Matrices*)


(* ::Text:: *)
(*Domain 0*)


fn$D\[Chi]="OperatorMatrix/Dchi_dom0_"<>"N_"<>ToString[Nz$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";
fn$D2\[Chi]="OperatorMatrix/D2chi_dom0_"<>"N_"<>ToString[Nz$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";

check$Dx=FileExistsQ[fn$D\[Chi]];
check$D2x=FileExistsQ[fn$D2\[Chi]];

If[!check$Dx || !check$D2x, 
(*Chebyshev Differeniation Matrix*)
Print["\t\tCalculating Differentiation Matrices"];
D\[Chi]$dom0=DiffMatrix[-1,1,Nz$dom0,Prec]; 
D2\[Chi]$dom0=D\[Chi]$dom0 . D\[Chi]$dom0;
Export[fn$D\[Chi],IntegerPart[10^(Prec+10)*D\[Chi]$dom0],"Table"];
Export[fn$D2\[Chi],IntegerPart[10^(Prec+10)*D2\[Chi]$dom0],"Table"];,

Print["\t\tLoading Differentiation Matrices"];
D\[Chi]$dom0=N[Import[fn$D\[Chi],"Table"]/10^(Prec+10),Prec];
D2\[Chi]$dom0=N[Import[fn$D2\[Chi],"Table"]/10^(Prec+10),Prec];
]

Id$dom0=IdentityMatrix[nz$dom0];
Zero$dom0=ConstantArray[0,{nz$dom0,nz$dom0}];
Zero$dom0dom1=ConstantArray[0,{nz$dom0,nz$dom1}];


(* ::Text:: *)
(*Domain 1*)


fn$D\[Chi]="OperatorMatrix/Dchi_dom1_"<>"N_"<>ToString[Nz$dom1]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";
fn$D2\[Chi]="OperatorMatrix/D2chi_dom1_"<>"N_"<>ToString[Nz$dom1]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";

check$Dx=FileExistsQ[fn$D\[Chi]];
check$D2x=FileExistsQ[fn$D2\[Chi]];

If[!check$Dx || !check$D2x, 
(*Chebyshev Differeniation Matrix*)
Print["\t\tCalculating Differentiation Matrices"];
D\[Chi]$dom1=DiffMatrix[-1,1,Nz$dom1,Prec]; 
D2\[Chi]$dom1=D\[Chi]$dom1 . D\[Chi]$dom1;
Export[fn$D\[Chi],IntegerPart[10^(Prec+10)*D\[Chi]$dom1],"Table"];
Export[fn$D2\[Chi],IntegerPart[10^(Prec+10)*D2\[Chi]$dom1],"Table"];,

Print["\t\tLoading Differentiation Matrices"];
D\[Chi]$dom1=N[Import[fn$D\[Chi],"Table"]/10^(Prec+10),Prec];
D2\[Chi]$dom1=N[Import[fn$D2\[Chi],"Table"]/10^(Prec+10),Prec];
]

Id$dom1=IdentityMatrix[nz$dom1];
Zero$dom1=ConstantArray[0,{nz$dom1,nz$dom1}];
Zero$dom1dom0=ConstantArray[0,{nz$dom1,nz$dom0}];


(* ::Subsubsection:: *)
(*Map derivatives \[Delta]\[Chi] into \[Delta]x*)


(* ::Text:: *)
(*Domain 0*)


If[\[Kappa]$dom0==0,
Dx$dom0=D\[Chi]$dom0;
D2x$dom0=D2\[Chi]$dom0;

DxRef$dom0=D\[Chi]Ref$dom0;
D2xRef$dom0=D2\[Chi]Ref$dom0;

,

J =N[ 2*\[Kappa]$dom0*Cosh[\[Kappa]$dom0*(1-\[Chi]$dom0*xb$dom0)]/Sinh[2*\[Kappa]$dom0],Prec];
J2 =N[ -2*\[Kappa]$dom0^2*xb$dom0*Sinh[\[Kappa]$dom0*(1-\[Chi]$dom0*xb$dom0)]/Sinh[2*\[Kappa]$dom0],Prec];

Dx$dom0 = N[D\[Chi]$dom0/J,Prec];
D2x$dom0 = N[(D2\[Chi]$dom0 - J2*Dx$dom0)/J^2,Prec];

JRef =N[ 2*\[Kappa]$dom0*Cosh[\[Kappa]$dom0*(1-\[Chi]Ref$dom0*xb$dom0)]/Sinh[2*\[Kappa]$dom0],Prec];
J2Ref =N[ -2*\[Kappa]$dom0^2*xb$dom0*Sinh[\[Kappa]$dom0*(1-\[Chi]Ref$dom0*xb$dom0)]/Sinh[2*\[Kappa]$dom0],Prec];

DxRef$dom0 = N[D\[Chi]Ref$dom0/JRef,Prec];
D2xRef$dom0 = N[(D2\[Chi]Ref$dom0 - J2Ref*DxRef$dom0)/JRef^2,Prec];
];


(* ::Text:: *)
(*Domain 1*)


If[\[Kappa]$dom1==0,
Dx$dom1=D\[Chi]$dom1;
D2x$dom1=D2\[Chi]$dom1;

DxRef$dom1=D\[Chi]Ref$dom1;
D2xRef$dom1=D2\[Chi]Ref$dom1
,

J =N[ 2*\[Kappa]$dom1*Cosh[\[Kappa]$dom1*(1-\[Chi]$dom1*xb$dom1)]/Sinh[2*\[Kappa]$dom1],Prec];
J2 =N[ -2*\[Kappa]$dom1^2*xb$dom1*Sinh[\[Kappa]$dom1*(1-\[Chi]$dom1*xb$dom1)]/Sinh[2*\[Kappa]$dom1],Prec];

Dx$dom1 = N[D\[Chi]$dom1/J,Prec];
D2x$dom1 = N[(D2\[Chi]$dom1 - J2*Dx$dom1)/J^2,Prec];

JRef =N[ 2*\[Kappa]$dom1*Cosh[\[Kappa]$dom1*(1-\[Chi]Ref$dom1*xb$dom1)]/Sinh[2*\[Kappa]$dom1],Prec];
J2Ref =N[ -2*\[Kappa]$dom1^2*xb$dom1*Sinh[\[Kappa]$dom1*(1-\[Chi]Ref$dom1*xb$dom1)]/Sinh[2*\[Kappa]$dom1],Prec];

DxRef$dom1 = N[D\[Chi]Ref$dom1/JRef,Prec];
D2xRef$dom1 = N[(D2\[Chi]Ref$dom1 - J2Ref*DxRef$dom1)/JRef^2,Prec];
];


(* ::Subsubsection:: *)
(*Re-scale Derivative Operators*)


D\[Sigma]$dom0=2/\[CapitalDelta]z$dom0 * Dx$dom0;
D2\[Sigma]$dom0=(2/\[CapitalDelta]z$dom0)^2 * D2x$dom0;

D\[Sigma]Ref$dom0=2/\[CapitalDelta]z$dom0 * DxRef$dom0;
D2\[Sigma]Ref$dom0=(2/\[CapitalDelta]z$dom0)^2 * D2xRef$dom0;


D\[Sigma]$dom1=2/\[CapitalDelta]z$dom1 * Dx$dom1;
D2\[Sigma]$dom1=(2/\[CapitalDelta]z$dom1)^2 * D2x$dom1;

D\[Sigma]Ref$dom1=2/\[CapitalDelta]z$dom1 * DxRef$dom1;
D2\[Sigma]Ref$dom1=(2/\[CapitalDelta]z$dom1)^2 * D2xRef$dom1;


(* ::Subsubsection:: *)
(*TEST OPERATORS*)


(* ::Text:: *)
(*(SET TO NON-EVALUATED: CHANGE IF NEED TO TEST)*)


(*fTest[\[Sigma]_]:=Cos[2*\[Pi]*\[Sigma]]*Exp[-\[Sigma]^2]*)
(*fVec$dom0=fTest[z$dom0];*)
(*fVec$dom1=fTest[z$dom1];*)
(**)
(*Error$d\[Sigma]$dom0=Total[D\[Sigma]$dom0 . fVec$dom0-fTest'[z$dom0]]*)
(*Error$d\[Sigma]$dom1=Total[D\[Sigma]$dom1 . fVec$dom1-fTest'[z$dom1]]*)
(**)
(**)
(*Error$d\[Sigma]d\[Sigma]$dom0=Total[D\[Sigma]$dom0 . D\[Sigma]$dom0 . fVec$dom0-fTest''[z$dom0]]*)
(*Error$d\[Sigma]d\[Sigma]$dom1=Total[D\[Sigma]$dom1 . D\[Sigma]$dom1 . fVec$dom1-fTest''[z$dom1]]*)
(**)
(**)
(*Error$d2\[Sigma]$dom0=Total[D2\[Sigma]$dom0 . fVec$dom0-fTest''[z$dom0]]*)
(*Error$d2\[Sigma]$dom1=Total[D2\[Sigma]$dom1 . fVec$dom1-fTest''[z$dom1]]*)


(* ::Subsection::Closed:: *)
(*Discrete Functions*)


\[Alpha]2$dom0=\[Alpha]2/.\[Sigma]->z$dom0;
\[Alpha]1$dom0=\[Alpha]1/.\[Sigma]->z$dom0;
\[Alpha]0$dom0=Limit[\[Alpha]0, \[Sigma]->z$dom0, Direction->"FromAbove"];

\[Beta]1$dom0=\[Beta]1/.\[Sigma]->z$dom0;
\[Beta]0$dom0=\[Beta]0/.\[Sigma]->z$dom0;
W$dom0=W/.\[Sigma]->z$dom0;


\[Alpha]2$dom1=\[Alpha]2/.\[Sigma]->z$dom1;
\[Alpha]1$dom1=\[Alpha]1/.\[Sigma]->z$dom1;
\[Alpha]0$dom1=Limit[\[Alpha]0, \[Sigma]->z$dom1, Direction->"FromBelow"];

\[Beta]1$dom1=\[Beta]1/.\[Sigma]->z$dom1;
\[Beta]0$dom1=\[Beta]0/.\[Sigma]->z$dom1;

W$dom1=W/.\[Sigma]->z$dom1;


(* ::Section:: *)
(*QNM Eigenvalue Problem (Load Data)*)


(* ::Subsection:: *)
(*Discrete ODE operators *)
(*( - NOTE: TIME DOMAIN WAVE EQUATION IS NOT DIVIDED BY w(\[Sigma]);*)
(*  - L matrix is not needed as we wil import QNM data )*)


L1$dom0=(\[Alpha]2$dom0*D2\[Sigma]$dom0 + \[Alpha]1$dom0*D\[Sigma]$dom0 + \[Alpha]0$dom0*Id$dom0);
L2$dom0=( \[Beta]1$dom0*D\[Sigma]$dom0 + \[Beta]0$dom0*Id$dom0);

L1$dom1=(\[Alpha]2$dom1*D2\[Sigma]$dom1 + \[Alpha]1$dom1*D\[Sigma]$dom1 + \[Alpha]0$dom1*Id$dom1);
L2$dom1=( \[Beta]1$dom1*D\[Sigma]$dom1 + \[Beta]0$dom1*Id$dom1);


(* ::Subsection:: *)
(*Load Eigenvalues and Eigenfunctions generated with BH_QNM_Eigenvalue.nb*)


n$QNM$max=N$QNM$max+1;
n$QNM=Range[0,N$QNM$max];

DirName$InputData="QNM_Eigenfunctions/eps_"<>ToString[N[\[Epsilon]],CForm]<>"_a0_Over_rh_"<>ToString[N[ao],CForm];


Nz=Nz$QNMload$dom0;
fn=Table[DirName$InputData<>"/Cheb_QNM_dom0_N_"<>ToString[Nz]<>"_n"<>ToString[n$QNM[[iQNM]]]<>"_l"<>ToString[l]<>"_spin"<>ToString[spin]<>"_Prec_"<>ToString[N[Prec]]<>".dat",{iQNM,1,n$QNM$max}]

QNM$EigensytemData$dom0=Table[Import[fn[[iQNM]],"Table"],{iQNM,1,n$QNM$max}];
QNM$EigenFunctionsData$dom0=Table[Delete[QNM$EigensytemData$dom0[[iQNM]],{{1},{2},{3},{4}}],{iQNM,1,n$QNM$max}];

Print["QNM s"]
s$QNM = Table[QNM$EigensytemData$dom0[[iQNM,1,2]]+I*QNM$EigensytemData$dom0[[iQNM,1,3]],{iQNM,1,n$QNM$max}];
N[s$QNM]
s$QNM$Conj=Conjugate[s$QNM];

N$QNM$dom0=Table[QNM$EigensytemData$dom0[[iQNM,2,2]],{iQNM,1,n$QNM$max}];

Print["Chebyschev Coefficients QNM Eigenfunctions"]
cRE\[Phi]$dom0=Table[QNM$EigenFunctionsData$dom0[[iQNM,All,2]],{iQNM,1,n$QNM$max}];
cIM\[Phi]$dom0=Table[QNM$EigenFunctionsData$dom0[[iQNM,All,3]],{iQNM,1,n$QNM$max}];

Print["QNM Eigenfunctions"]
\[Phi]QNM$Func$dom0[\[Sigma]_, n_]:=ChebInterpolation[cRE\[Phi]$dom0[[n+1]],zi,zp,\[Sigma], \[Kappa]$dom0, xb$dom0, Prec]+I*ChebInterpolation[cIM\[Phi]$dom0[[n+1]],zi,zp,\[Sigma], \[Kappa]$dom0, xb$dom0,Prec]

Plot$\[Phi]QNM$Func$dom0=Table[Plot[Re@\[Phi]QNM$Func$dom0[\[Sigma], n],{\[Sigma],zi,zp}, PlotRange->{{zi,zf},All}, PlotStyle->ColorData[n+1, "ColorFunction"]], {n,0,N$QNM$max}];

\[Phi]QNM$dom0=Table[\[Phi]QNM$Func$dom0[z$dom0,n],{n,0,N$QNM$max}];
\[Phi]QNM$Conj$dom0=Conjugate[\[Phi]QNM$dom0];


Nz=Nz$QNMload$dom1;

fn=Table[DirName$InputData<>"/Cheb_QNM_dom1_N_"<>ToString[Nz]<>"_n"<>ToString[n$QNM[[iQNM]]]<>"_l"<>ToString[l]<>"_spin"<>ToString[spin]<>"_Prec_"<>ToString[N[Prec]]<>".dat",{iQNM,1,n$QNM$max}]

QNM$EigensytemData$dom1=Table[Import[fn[[iQNM]],"Table"],{iQNM,1,n$QNM$max}];
QNM$EigenFunctionsData$dom1=Table[Delete[QNM$EigensytemData$dom1[[iQNM]],{{1},{2},{3},{4}}],{iQNM,1,n$QNM$max}];

Print["QNM s"]
s$QNM = Table[QNM$EigensytemData$dom1[[iQNM,1,2]]+I*QNM$EigensytemData$dom1[[iQNM,1,3]],{iQNM,1,n$QNM$max}];
N[s$QNM]
s$QNM$Conj=Conjugate[s$QNM];

N$QNM$dom1=Table[QNM$EigensytemData$dom1[[iQNM,2,2]],{iQNM,1,n$QNM$max}];

Print["Chebyschev Coefficients QNM Eigenfunctions"]
cRE\[Phi]$dom1=Table[QNM$EigenFunctionsData$dom1[[iQNM,All,2]],{iQNM,1,n$QNM$max}];
cIM\[Phi]$dom1=Table[QNM$EigenFunctionsData$dom1[[iQNM,All,3]],{iQNM,1,n$QNM$max}];

Print["QNM Eigenfunctions"]
\[Phi]QNM$Func$dom1[\[Sigma]_, n_]:=ChebInterpolation[cRE\[Phi]$dom1[[n+1]],zp,zf,\[Sigma], \[Kappa]$dom1, xb$dom1,Prec]+I*ChebInterpolation[cIM\[Phi]$dom1[[n+1]],zp,zf,\[Sigma], \[Kappa]$dom1, xb$dom1,Prec]

\[Phi]QNM$dom1=Table[\[Phi]QNM$Func$dom1[z$dom1,n],{n,0,N$QNM$max}];
\[Phi]QNM$Conj$dom1=Conjugate[\[Phi]QNM$dom1];


(* ::Chapter:: *)
(*QNM Amplitudes*)


(* ::Section:: *)
(*QNM Excitation Operators*)


(* ::Subsection:: *)
(*Operators*)


A$dom0[s_]:=- W$dom0*s^2*Id$dom0 + L1$dom0 + s*L2$dom0;
c$dom0[s_, sQNM_]:=L2$dom0- W$dom0*(s+sQNM)*Id$dom0;
B$dom0[s_,V0_,W0_]:=L2$dom0 . V0 - W$dom0*(s*V0 + W0);

A$dom1[s_]:=- W$dom1*s^2*Id$dom1 + L1$dom1 + s*L2$dom1;
c$dom1[s_, sQNM_]:=L2$dom1- W$dom1*(s+sQNM)*Id$dom1;
B$dom1[s_,V0_,W0_]:=L2$dom1 . V0 - W$dom1*(s*V0 + W0);


(* ::Subsubsection:: *)
(*Choose QNM*)


An$dom0=Table[A$dom0[s$QNM[[iQNM]]],{iQNM,1,n$QNM$max}];
Cn$dom0=Table[c$dom0[s$QNM[[iQNM]],s$QNM[[iQNM]]] . \[Phi]QNM$dom0[[iQNM]],{iQNM,1,n$QNM$max}];

An$dom1=Table[A$dom1[s$QNM[[iQNM]]],{iQNM,1,n$QNM$max}];
Cn$dom1=Table[c$dom1[s$QNM[[iQNM]],s$QNM[[iQNM]]] . \[Phi]QNM$dom1[[iQNM]],{iQNM,1,n$QNM$max}];

(*An$Conj=Table[A[s$QNM$Conj[[iQNM]]],{iQNM,1,n$QNM$max}];
Cn$Conj=Table[c[s$QNM$Conj[[iQNM]],s$QNM$Conj[[iQNM]]].\[Phi]QNM$Conj$Vec[[iQNM]],{iQNM,1,n$QNM$max}];*)


(* ::Subsubsection:: *)
(*Multidomain Matrix*)


An=Table[ N[ArrayFlatten[{
{   An$dom0[[iQNM]]    ,       Zero$dom0dom1   } , 
{  Zero$dom1dom0,   An$dom1[[iQNM]] } 
}],
Prec],
{iQNM,1,n$QNM$max}
];




Cn=Table[ Join[Cn$dom0[[iQNM]],Cn$dom1[[iQNM]]],{iQNM,1,n$QNM$max}];


(* ::Subsubsection:: *)
(*Continuity of Auxiliary function at Boundary Domain*)


An[[All,1,All]]=0;
An[[All,1,1]]=1;
An[[All,1,-1]]=-1;

Cn[[All,1]]=0;


(* ::Subsubsection:: *)
(*Continuity of Auxiliary function's derivative at Boundary Domain*)


An[[All,nz$dom0 + nz$dom1,All]]=0;
An[[All,nz$dom0 + nz$dom1, ;;nz$dom0]]=D\[Sigma]$dom0[[1, All]];
An[[All,nz$dom0 + nz$dom1, nz$dom0+1;;nz$dom0+nz$dom1 ]]=-D\[Sigma]$dom1[[-1, All]];

Cn[[All,-1]]=0;


(* ::Subsection:: *)
(*Amplitude Matrix*)


m1t =Table[Join[Transpose@An[[iQNM]], {Cn[[iQNM]]}],{iQNM,1,n$QNM$max}];
m1 =Table[Transpose[m1t[[iQNM]]],{iQNM,1,n$QNM$max}];

\[Delta]=Table[KroneckerDelta[i,nz$dom0+1],{i,1,nz$dom0+nz$dom1+1}];

Mamp=Table[Join[m1[[iQNM]], {\[Delta]}],{iQNM,1,n$QNM$max}];
(*MampInv=ParallelTable[Inverse@Mamp[[iQNM]],{iQNM,1,n$QNM$max}];*)


(*m1t$Conj=Table[Join[Transpose@An$Conj[[iQNM]], {Cn$Conj[[iQNM]]}],{iQNM,1,n$QNM$max}];
m1$Conj =Table[Transpose[m1t$Conj[[iQNM]]],{iQNM,1,n$QNM$max}];

Mamp$Conj=Table[Join[m1$Conj[[iQNM]], {\[Delta]}],{iQNM,1,n$QNM$max}];
Mamp$ConjInv=ParallelTable[Inverse@Mamp$Conj[[iQNM]],{iQNM,1,n$QNM$max}];*)


(* ::Section:: *)
(*Calculate QNM Amplitude *)


(* ::Subsection:: *)
(*Initial Data*)


\[Phi]0$dom0=ConstantArray[1,nz$dom0]; (*2*Re[\[Phi]QNM$dom0];*) 
\[Psi]0$dom0=ConstantArray[0,nz$dom0];(*2*Re[\[Phi]QNM$dom0*s$QNM];*) 

\[Phi]0$dom1=ConstantArray[1,nz$dom1]; (*2*Re[\[Phi]QNM$dom1];*) 
\[Psi]0$dom1= ConstantArray[0,nz$dom1];(*2*Re[\[Phi]QNM$dom1*s$QNM];*) 


(* ::Subsubsection:: *)
(*Source Function ID*)


Bn$dom0=Table[B$dom0[s$QNM[[iQNM]],\[Phi]0$dom0,\[Psi]0$dom0],{iQNM,1,n$QNM$max}];
Bn$dom1=Table[B$dom1[s$QNM[[iQNM]],\[Phi]0$dom1,\[Psi]0$dom1],{iQNM,1,n$QNM$max}];

Bn=Table[Join[Bn$dom0[[iQNM]],Bn$dom1[[iQNM]]],{iQNM,1,n$QNM$max}];
(*Bn$Conj=Table[B[s$QNM$Conj[[iQNM]],\[Phi]0,\[Psi]0],{iQNM,1,n$QNM$max}];*)

Bn[[All,1]]=0;
Bn[[All,-1]]=0;


(* ::Subsection:: *)
(*Total Source*)


g0=1;

S=Table[Append[Bn[[iQNM]] , g0],{iQNM,1,n$QNM$max}];
(*S$Conj=Table[Append[Bn$Conj[[iQNM]] +S2ndOrder$Vec[s$QNM$Conj[[iQNM]]], g0],{iQNM,1,n$QNM$max}];*)


(* ::Subsubsection:: *)
(*Amplitude*)


\[Eta]Sol=Table[LinearSolve[Mamp[[iQNM]] , S[[iQNM]]],{iQNM,1,n$QNM$max}];
(*Table[MampInv[[iQNM]] . S[[iQNM]],{iQNM,1,n$QNM$max}];*)

\[Eta]=\[Eta]Sol[[All,nz$dom0+nz$dom1+1]];

\[Phi]Scri=\[Phi]QNM$dom0[[All,nz$dom0]];
\[Phi]Hrz=\[Phi]QNM$dom1[[All,1]];

\[Eta]\[Phi]Scri=Table[ \[Eta][[iQNM]] * \[Phi]Scri[[iQNM]],{iQNM,1,n$QNM$max}];
\[Eta]\[Phi]Hrz =Table[ \[Eta][[iQNM]] * \[Phi]Hrz[[iQNM]],{iQNM,1,n$QNM$max}];

\[Eta]\[Phi]Scri;
Abs[\[Eta]\[Phi]Hrz];


(* ::Chapter:: *)
(*Export Data*)


(* ::Section:: *)
(*Map QNM into regular QNM*)


sData=Table[Null,{i,1,n$QNM$max}];
AScriData=Table[Null,{i,1,n$QNM$max}];
AHrzData=Table[Null,{i,1,n$QNM$max}];

For[i=1, i<=n$QNM$max, i++,
If[Im[s$QNM[[i]]]<0, 

sData[[i]]=s$QNM[[i]];
AScriData[[i]]=\[Eta]\[Phi]Scri[[i]];
AHrzData[[i]]=\[Eta]\[Phi]Hrz[[i]];

, 
sData[[i]]=Conjugate[s$QNM[[i]]];
AScriData[[i]]=Conjugate[\[Eta]\[Phi]Scri[[i]]];
AHrzData[[i]]=Conjugate[\[Eta]\[Phi]Hrz[[i]]];
]
]


(* ::Section:: *)
(*Export Data*)


Print["Exporting Data..."]

DirName$InputData="QNM_Amplitudes/eps_"<>ToString[N[\[Epsilon]],CForm]<>"_a0_Over_rh_"<>ToString[N[ao],CForm];
If[!DirectoryQ[DirName$InputData],CreateDirectory[DirName$InputData]]


fn=DirName$InputData<>"/QNM_Amplitude_Scri_N_"<>ToString[Nz]<>"_nQNM"<>ToString[n$QNM$max]<>"_l"<>ToString[l]<>"_spin"<>ToString[spin]<>"_Prec_"<>ToString[N[Prec]]<>".dat";

fp=OpenWrite[fn];
WriteString[fp,"#1:Real(s)\t #2:Imag(s)\t #3:Re(A)\t #4:Im(A)\n"]
For[i=1, i<=n$QNM$max, i++,
WriteString[fp,
ToString[Re[sData[[i]]],CForm],"\t", ToString[Im[sData[[i]]],CForm],"\t",ToString[Re[AScriData[[i]]],CForm],"\t", ToString[Im[AScriData[[i]]],CForm],"\n"
]

] ;
Close[fp];



fn=DirName$InputData<>"/QNM_Amplitude_Hrzn_N_"<>ToString[Nz]<>"_nQNM"<>ToString[n$QNM$max]<>"_l"<>ToString[l]<>"_spin"<>ToString[spin]<>"_Prec_"<>ToString[N[Prec]]<>".dat";

fp=OpenWrite[fn];
WriteString[fp,"#1:Real(s)\t #2:Imag(s)\t #3:Re(A)\t #4:Im(A)\n"]
For[i=1, i<=n$QNM$max, i++,
WriteString[fp,
ToString[Re[sData[[i]]],CForm],"\t", ToString[Im[sData[[i]]],CForm],"\t",ToString[Re[AHrzData[[i]]],CForm],"\t", ToString[Im[AHrzData[[i]]],CForm],"\n"
]

] ;
Close[fp];
