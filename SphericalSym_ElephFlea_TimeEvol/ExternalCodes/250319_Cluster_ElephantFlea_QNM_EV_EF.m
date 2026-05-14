(* ::Package:: *)

(* ::Chapter:: *)
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
(* ao=10; (*120000/10000;*)(*115722/10000;*)
 \[Epsilon]=204/100000;*)


(* ::Section:: *)
(*Numerical Parameters*)


(*Prec=10*MachinePrecision;

Nz$dom0=75;
Nz$dom1=75;

NzRef$dom0=70;
NzRef$dom1=70;
\[Epsilon]QNM=10^(-4); (*Relative error between QNMs from different resolutions*)*)


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


(* ::Section:: *)
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


(* ::Section:: *)
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


(* ::Chapter:: *)
(*Eigenvalue Solver*)


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


(* ::Subsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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


nzRef$dom0=NzRef$dom0+1; (*Total number of points*)
nzRef$dom1=NzRef$dom1+1; (*Total number of points*)


\[Chi]Ref$dom0=N[Table[x[i,NzRef$dom0],{i,0,NzRef$dom0}],Prec];
\[Chi]Ref$dom1=N[Table[x[i,NzRef$dom1],{i,0,NzRef$dom1}],Prec];

ZeroArrayRef$dom0=ConstantArray[0,nzRef$dom0];
ZeroArrayRef$dom1=ConstantArray[0,nzRef$dom1];


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


(* ::Subsection:: *)
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


fnRef$D\[Chi]="OperatorMatrix/Dchi_dom0_"<>"N_"<>ToString[NzRef$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";
fnRef$D2\[Chi]="OperatorMatrix/D2chi_dom0_"<>"N_"<>ToString[NzRef$dom0]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";

check$Dx=FileExistsQ[fnRef$D\[Chi]];
check$D2x=FileExistsQ[fnRef$D2\[Chi]];

If[!check$Dx || !check$D2x, 
(*Chebyshev Differeniation Matrix*)
Print["\t\tCalculating Differentiation Matrices"];
D\[Chi]Ref$dom0=DiffMatrix[-1,1,NzRef$dom0,Prec]; 
D2\[Chi]Ref$dom0=D\[Chi]Ref$dom0 . D\[Chi]Ref$dom0;
Export[fnRef$D\[Chi],IntegerPart[10^(Prec+10)*D\[Chi]Ref$dom0],"Table"];
Export[fnRef$D2\[Chi],IntegerPart[10^(Prec+10)*D2\[Chi]Ref$dom0],"Table"];,

Print["\t\tLoading Differentiation Matrices"];
D\[Chi]Ref$dom0=N[Import[fnRef$D\[Chi],"Table"]/10^(Prec+10),Prec];
D2\[Chi]Ref$dom0=N[Import[fnRef$D2\[Chi],"Table"]/10^(Prec+10),Prec];
]

IdRef$dom0=IdentityMatrix[nzRef$dom0];
ZeroRef$dom0=ConstantArray[0,{nzRef$dom0,nzRef$dom0}];
ZeroRef$dom0dom1=ConstantArray[0,{nzRef$dom0,nzRef$dom1}];


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

fnRef$D\[Chi]="OperatorMatrix/Dchi_dom1_"<>"N_"<>ToString[NzRef$dom1]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";
fnRef$D2\[Chi]="OperatorMatrix/D2chi_dom1_"<>"N_"<>ToString[NzRef$dom1]<>"_Prec_"<>ToString[Floor[Prec]]<>".dat";

check$Dx=FileExistsQ[fnRef$D\[Chi]];
check$D2x=FileExistsQ[fnRef$D2\[Chi]];

If[!check$Dx || !check$D2x, 
(*Chebyshev Differeniation Matrix*)
Print["\t\tCalculating Differentiation Matrices"];
D\[Chi]Ref$dom1=DiffMatrix[-1,1,NzRef$dom1,Prec]; 
D2\[Chi]Ref$dom1=D\[Chi]Ref$dom1 . D\[Chi]Ref$dom1;
Export[fnRef$D\[Chi],IntegerPart[10^(Prec+10)*D\[Chi]Ref$dom1],"Table"];
Export[fnRef$D2\[Chi],IntegerPart[10^(Prec+10)*D2\[Chi]Ref$dom1],"Table"];,

Print["\t\tLoading Differentiation Matrices"];
D\[Chi]Ref$dom1=N[Import[fnRef$D\[Chi],"Table"]/10^(Prec+10),Prec];
D2\[Chi]Ref$dom1=N[Import[fnRef$D2\[Chi],"Table"]/10^(Prec+10),Prec];
]

IdRef$dom1=IdentityMatrix[nzRef$dom1];
ZeroRef$dom1=ConstantArray[0,{nzRef$dom1,nzRef$dom1}];
ZeroRef$dom1dom0=ConstantArray[0,{nzRef$dom1,nzRef$dom0}];


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


p$dom0=p/.\[Sigma]->z$dom0;
p$dom1=p/.\[Sigma]->z$dom1;

dp$dom0=D[p,\[Sigma]]/.\[Sigma]->z$dom0;
dp$dom1=D[p,\[Sigma]]/.\[Sigma]->z$dom1;

\[Gamma]$dom0=\[Gamma]/.\[Sigma]->z$dom0;
\[Gamma]$dom1=\[Gamma]/.\[Sigma]->z$dom1;

d\[Gamma]$dom0=D[\[Gamma],\[Sigma]]/.\[Sigma]->z$dom0;
d\[Gamma]$dom1=D[\[Gamma],\[Sigma]]/.\[Sigma]->z$dom1;

w$dom0=w/.\[Sigma]->z$dom0;
w$dom1=w/.\[Sigma]->z$dom1;

q$dom0=Limit[q,\[Sigma]->z$dom0, Direction->"FromAbove"];
q$dom1=Limit[q,\[Sigma]->z$dom1, Direction->"FromBelow"];


pRef$dom0=p/.\[Sigma]->zRef$dom0;
pRef$dom1=p/.\[Sigma]->zRef$dom1;

dpRef$dom0=D[p,\[Sigma]]/.\[Sigma]->zRef$dom0;
dpRef$dom1=D[p,\[Sigma]]/.\[Sigma]->zRef$dom1;

\[Gamma]Ref$dom0=\[Gamma]/.\[Sigma]->zRef$dom0;
\[Gamma]Ref$dom1=\[Gamma]/.\[Sigma]->zRef$dom1;

d\[Gamma]Ref$dom0=D[\[Gamma],\[Sigma]]/.\[Sigma]->zRef$dom0;
d\[Gamma]Ref$dom1=D[\[Gamma],\[Sigma]]/.\[Sigma]->zRef$dom1;

wRef$dom0=w/.\[Sigma]->zRef$dom0;
wRef$dom1=w/.\[Sigma]->zRef$dom1;

qRef$dom0=Limit[q,\[Sigma]->zRef$dom0, Direction->"FromAbove"];
qRef$dom1=Limit[q,\[Sigma]->zRef$dom1, Direction->"FromBelow"];

(*cq$dom0=ChebLobattoCoefficients[q$dom0,Prec];
cq$dom1=ChebLobattoCoefficients[q$dom1,Prec];
ListLogPlot[{Abs@cq$dom0,Abs@cq$dom1},PlotLabel->"Cheb Coeff func Pot", PlotStyle->{Red,Purple}]*)


(*Plot$p=Plot[p,{\[Sigma],0,1}];
Plot$pVec$dom0=ListPlot[Transpose@{z$dom0,p$dom0}, PlotStyle->Red];
Plot$pVec$dom1=ListPlot[Transpose@{z$dom1,p$dom1}, PlotStyle->Purple];
Show[Plot$p,Plot$pVec$dom0,Plot$pVec$dom1]

Plot$\[Gamma]=Plot[\[Gamma],{\[Sigma],0,1}];
Plot$\[Gamma]Vec$dom0=ListPlot[Transpose@{z$dom0,\[Gamma]$dom0}, PlotStyle->Red];
Plot$\[Gamma]Vec$dom1=ListPlot[Transpose@{z$dom1,\[Gamma]$dom1}, PlotStyle->Purple];
Show[Plot$\[Gamma],Plot$\[Gamma]Vec$dom0,Plot$\[Gamma]Vec$dom1]

Plot$w=Plot[w//Simplify,{\[Sigma],0,1}, WorkingPrecision->100];
Plot$wVec$dom0=ListPlot[Transpose@{z$dom0,w$dom0}, PlotStyle->Red];
Plot$wVec$dom1=ListPlot[Transpose@{z$dom1,w$dom1}, PlotStyle->Purple];
Show[Plot$w,Plot$wVec$dom0,Plot$wVec$dom1]

Plot$q=Plot[q//Simplify,{\[Sigma],0,1}, WorkingPrecision->100];
Plot$qVec$dom0=ListPlot[Transpose@{z$dom0,q$dom0}, PlotStyle->Red];
Plot$qVec$dom1=ListPlot[Transpose@{z$dom1,q$dom1}, PlotStyle->Purple];
Show[Plot$q,Plot$qVec$dom0,Plot$qVec$dom1]*)


(* ::Subsubsection:: *)
(*Chebyshev Coefficients*)


(*cp$dom0=ChebLobattoCoefficients[p$dom0,Prec];
cp$dom1=ChebLobattoCoefficients[p$dom1,Prec];
ListLogPlot[{Abs@cp$dom0,Abs@cp$dom1}, PlotLabel->"Cheb Coeff func p", PlotStyle->{Red,Purple}]

c\[Gamma]$dom0=ChebLobattoCoefficients[\[Gamma]$dom0,Prec];
c\[Gamma]$dom1=ChebLobattoCoefficients[\[Gamma]$dom1,Prec];
ListLogPlot[{Abs@c\[Gamma]$dom0, Abs@c\[Gamma]$dom1},PlotLabel->"Cheb Coeff func \[Gamma]", PlotStyle->{Red,Purple}]

cw$dom0=ChebLobattoCoefficients[w$dom0,Prec];
cw$dom1=ChebLobattoCoefficients[w$dom1,Prec];
ListLogPlot[{Abs@cw$dom0,Abs@cw$dom1},PlotLabel->"Cheb Coeff func w", PlotStyle->{Red,Purple}]*)

(*cq$dom0=ChebLobattoCoefficients[q$dom0,Prec];
cq$dom1=ChebLobattoCoefficients[q$dom1,Prec];
ListLogPlot[{Abs@cq$dom0,Abs@cq$dom1},PlotLabel->"Cheb Coeff func Pot", PlotStyle->{Red,Purple}]*)


(* ::Section::Closed:: *)
(*QNM Eigenvalue Problem*)


(* ::Subsection:: *)
(*Discrete ODE operators*)


L1$dom0=p$dom0*D2\[Sigma]$dom0 + dp$dom0*D\[Sigma]$dom0 - q$dom0*Id$dom0;
L1$dom1=p$dom1*D2\[Sigma]$dom1 + dp$dom1*D\[Sigma]$dom1 - q$dom1*Id$dom1;

L2$dom0=2*\[Gamma]$dom0*D\[Sigma]$dom0 + d\[Gamma]$dom0*Id$dom0 ;
L2$dom1=2*\[Gamma]$dom1*D\[Sigma]$dom1 + d\[Gamma]$dom1*Id$dom1 ;

L=N[ArrayFlatten[{
{     Zero$dom0  ,       Id$dom0   ,  Zero$dom0dom1,     Zero$dom0dom1} , 
{  L1$dom0/w$dom0,   L2$dom0/w$dom0,  Zero$dom0dom1,     Zero$dom0dom1} ,
{     Zero$dom1dom0  ,       Zero$dom1dom0   ,  Zero$dom1,     Id$dom1},
{     Zero$dom1dom0  ,       Zero$dom1dom0   ,  L1$dom1/w$dom1,     L2$dom1/w$dom1}
}],
Prec];

B=IdentityMatrix[2 (nz$dom0+nz$dom1)];

L1Ref$dom0=pRef$dom0*D2\[Sigma]Ref$dom0 + dpRef$dom0*D\[Sigma]Ref$dom0 - qRef$dom0*IdRef$dom0;
L1Ref$dom1=pRef$dom1*D2\[Sigma]Ref$dom1 + dpRef$dom1*D\[Sigma]Ref$dom1 - qRef$dom1*IdRef$dom1;

L2Ref$dom0=2*\[Gamma]Ref$dom0*D\[Sigma]Ref$dom0 + d\[Gamma]Ref$dom0*IdRef$dom0 ;
L2Ref$dom1=2*\[Gamma]Ref$dom1*D\[Sigma]Ref$dom1 + d\[Gamma]Ref$dom1*IdRef$dom1 ;

LRef=N[ArrayFlatten[{
{     ZeroRef$dom0  ,       IdRef$dom0   ,  ZeroRef$dom0dom1,     ZeroRef$dom0dom1} , 
{  L1Ref$dom0/wRef$dom0,   L2Ref$dom0/wRef$dom0,  ZeroRef$dom0dom1,     ZeroRef$dom0dom1} ,
{     ZeroRef$dom1dom0  ,       ZeroRef$dom1dom0   ,  ZeroRef$dom1,     IdRef$dom1},
{     ZeroRef$dom1dom0  ,       ZeroRef$dom1dom0   ,  L1Ref$dom1/wRef$dom1,     L2Ref$dom1/wRef$dom1}
}],
Prec];

BRef=IdentityMatrix[2 (nzRef$dom0+nzRef$dom1)];


(* ::Subsubsection:: *)
(*Continuity of the function at the domain boundary*)


B[[nz$dom0+1,All]]=0;
L[[nz$dom0+1,All]]=0;
L[[nz$dom0+1,1]]=1;
L[[nz$dom0+1,2*nz$dom0+nz$dom1]]=-1;

BRef[[nzRef$dom0+1,All]]=0;
LRef[[nzRef$dom0+1,All]]=0;
LRef[[nzRef$dom0+1,1]]=1;
LRef[[nzRef$dom0+1,2*nzRef$dom0+nzRef$dom1]]=-1;


(* ::Subsubsection:: *)
(*Continuity of the function's Derivative at the domain boundary*)


B[[2*nz$dom0+nz$dom1+1,All]]=0;
L[[2*nz$dom0+nz$dom1+1,All]]=0;
L[[2*nz$dom0+nz$dom1+1, ;;nz$dom0]]=D\[Sigma]$dom0[[1, All]];
L[[2*nz$dom0+nz$dom1+1, 2*nz$dom0+1;;2*nz$dom0+nz$dom1]]=-D\[Sigma]$dom1[[-1, All]];

BRef[[2*nzRef$dom0+nzRef$dom1+1,All]]=0;
LRef[[2*nzRef$dom0+nzRef$dom1+1,All]]=0;
LRef[[2*nzRef$dom0+nzRef$dom1+1, ;;nzRef$dom0]]=D\[Sigma]Ref$dom0[[1, All]];
LRef[[2*nzRef$dom0+nzRef$dom1+1, 2*nzRef$dom0+1;;2*nzRef$dom0+nzRef$dom1]]=-D\[Sigma]Ref$dom1[[-1, All]];


(* ::Subsection:: *)
(*Solve Eigenvalue*)


Print["\t\tCalculating QNM values and Functions"];
EigenSol=Eigensystem[{L,B}];
Print["\t\tDone"];

Print["\t\tCalculating QNM References"];
EigenSolRef=Eigenvalues[{LRef,BRef}];
Print["\t\tDone"];


(* ::Subsubsection:: *)
(*Re-scale Laplace to Fourier*)


Print["\t Getting QNMs"];

\[Omega]=-I*EigenSol[[1]];
u=EigenSol[[2]];

\[Omega]Ref=-I*EigenSolRef;


Print["\t Getting QNMs Functions"];
\[Phi]QNM$dom0=Table[u[[All,i  ]],{i,1,nz$dom0}   ];
\[Phi]QNM$dom1=Table[u[[All,i  ]],{i,2*nz$dom0+1,nz$dom1}   ];
\[Phi]QNM=Join[\[Phi]QNM$dom0,\[Phi]QNM$dom1];


(* ::Subsubsection:: *)
(*Filter discrete from continuous spectra*)


QNM$filter[x_,y_,z_]:=Block[{list={}},tab=y;tab1=z;Do[Do[If[Abs[(tab[[-j]]-tab1[[-i]])/tab[[-j]]]<x
(*If[Abs[(tab[[-j]][[1]]-tab1[[-i]][[1]])/tab[[-j]][[1]]]<x\[And]Abs[(tab[[-j]][[2]]-tab1[[-i]][[2]])/tab[[-j]][[2]]]<x*),AppendTo[list,tab1[[-i]]]],{j,1,Length[tab]}],{i,1,Length[tab1]}];list]


Print["\t Filter QNMs Functions \[Epsilon] =", \[Epsilon]QNM];

PosQNMLong=Position[EigenSol[[1]],x_/; Im@x!=0 && Re@x < 0]//Flatten;
nQNMLong=Length@PosQNMLong;
\[Omega]$QNMLong=Table[\[Omega][[PosQNMLong[[i]]]],{i,1,nQNMLong}];
u$QNMLong=Table[u[[PosQNMLong[[i]]]],{i,1,nQNMLong}];

PosBranch=Position[EigenSol[[1]],x_/;Im@x==0]//Flatten;
nBranch=Length@PosBranch;
\[Omega]$BranchCut=Table[\[Omega][[PosBranch[[i]]]],{i,1,nBranch}];

BranchCutData=Table[{Re@\[Omega]$BranchCut[[iqnm]],Im@\[Omega]$BranchCut[[iqnm]]},{iqnm,1,nBranch}];


PosQNMRef=Position[EigenSolRef,x_/; Im@x!=0 && Re@x < 0]//Flatten;
nRefQNM=Length@PosQNMRef;
\[Omega]Ref$QNM=Table[\[Omega]Ref[[PosQNMRef[[i]]]],{i,1,nRefQNM}];

\[Omega]$QNM=QNM$filter[\[Epsilon]QNM,\[Omega]Ref$QNM,\[Omega]$QNMLong];
PosQNM=Position[\[Omega]$QNMLong, x_ /; MemberQ[\[Omega]$QNM, x]] // Flatten//Reverse

nQNM=Length@PosQNM;
QNMData=Table[{Re@\[Omega]$QNM[[iqnm]],Im@\[Omega]$QNM[[iqnm]]},{iqnm,1,nQNM}];
u$QNM=u$QNMLong[[PosQNM]];

\[Phi]QNM$Norm=Table[u$QNM[[iQNM,2*nz$dom0+1]],{iQNM,1,nQNM}];

\[Phi]QNM$dom0=Table[u$QNM[[iQNM,i]]/\[Phi]QNM$Norm[[iQNM]],{iQNM,1,nQNM},{i,1,nz$dom0}];
\[Phi]QNM$dom1=Table[u$QNM[[iQNM,i]]/\[Phi]QNM$Norm[[iQNM]],{iQNM,1,nQNM},{i,2*nz$dom0+1,2*nz$dom0+nz$dom1}];

\[Phi]QNM$dom0$data=Table[Transpose@{z$dom0,\[Phi]QNM$dom0[[iQNM]]},{iQNM,1,nQNM}];
\[Phi]QNM$dom1$data=Table[Transpose@{z$dom1,\[Phi]QNM$dom1[[iQNM]]},{iQNM,1,nQNM}];



ListPlot[QNMData , PlotRange->{{0.7,0.75}, {0,1}}]
PlotQNMFunc=ListPlot[{Re@\[Phi]QNM$dom0$data[[1]],Re@\[Phi]QNM$dom1$data[[1]]}, (*Joined->True,*) PlotRange->All]


(* ::Subsection:: *)
(*QNM Data for post-processing*)


Print["\t Post prodessing QNM"];
s$n0=I*(QNMData[[All,1]] + I*QNMData[[All,2]]);

cheb$\[Phi]QNM$dom0=Table[ChebLobattoCoefficients[\[Phi]QNM$dom0[[iQNM]],Prec],{iQNM,1,nQNM}];
cheb$\[Phi]QNM$dom1=Table[ChebLobattoCoefficients[\[Phi]QNM$dom1[[iQNM]],Prec],{iQNM,1,nQNM}];

ListLogPlot[{Abs@cheb$\[Phi]QNM$dom0[[1]],Abs@cheb$\[Phi]QNM$dom1[[1]] }]


(* ::Section:: *)
(*Export Data*)


(* ::Subsection:: *)
(*Initial Data Data C code Pure QNM Evolution*)


Print["\t Export Data"];
DirName$InputData="../input_data/eps_"<>ToString[N[\[Epsilon]],CForm]<>"_a0_Over_rh_"<>ToString[N[ao],CForm];
If[!DirectoryQ[DirName$InputData],CreateDirectory[DirName$InputData]]

For[iQNM=0, iQNM<nQNM/2, iQNM++,
Nz=Nz$dom0;

fn=DirName$InputData<>"/Cheb_QNM_dom0_N_"<>ToString[Nz]<>"_n"<>ToString[iQNM]<>"_l"<>ToString[l]<>"_spin"<>ToString[spin]<>"_Prec_"<>ToString[N[Prec]]<>".dat";
NumericalData={
"s\t"<>ToString[Re[s$n0[[2*iQNM+1]]]]<>"\t"<>ToString[Im[s$n0[[2*iQNM+1]]]]<>"\n",
"N\t"<>ToString[Nz]<>"\n\n"
};

fp=OpenWrite[fn];
Table[ WriteString[fp,NumericalData[[i]]], {i,1,Length@NumericalData}];

WriteString[fp,"#1:i\t #2:Real(Cheb_phi_QNM)\t #3:Imag(Cheb_phi_QNM)\n"]
For[i=0, i<=Nz, i++,
WriteString[fp,ToString[i],"\t",
ToString[Re[cheb$\[Phi]QNM$dom0[[2*iQNM+1,i+1]]],CForm],"\t", ToString[Im[cheb$\[Phi]QNM$dom0[[2*iQNM+1,i+1]]],CForm],"\n"]] ;
Close[fp];
fn=.;



Nz=Nz$dom1;
fn=DirName$InputData<>"/Cheb_QNM_dom1_N_"<>ToString[Nz]<>"_n"<>ToString[iQNM]<>"_l"<>ToString[l]<>"_spin"<>ToString[spin]<>"_Prec_"<>ToString[N[Prec]]<>".dat";
NumericalData={
"s\t"<>ToString[Re[s$n0[[2*iQNM+1]]]]<>"\t"<>ToString[Im[s$n0[[2*iQNM+1]]]]<>"\n",
"N\t"<>ToString[Nz]<>"\n\n"
};

fp=OpenWrite[fn];
Table[ WriteString[fp,NumericalData[[i]]], {i,1,Length@NumericalData}];

WriteString[fp,"#1:i\t #2:Real(Cheb_phi_QNM)\t #3:Imag(Cheb_phi_QNM)\n"]
For[i=0, i<=Nz, i++,
WriteString[fp,ToString[i],"\t",
ToString[Re[cheb$\[Phi]QNM$dom1[[2*iQNM+1,i+1]]],CForm],"\t", ToString[Im[cheb$\[Phi]QNM$dom1[[2*iQNM+1,i+1]]],CForm],"\n"]] 
Close[fp];
fn=.
]


(* ::Subsection:: *)
(*Amplitude ExtractionEvolution*)


DirName$InputData="QNM_Eigenfunctions/eps_"<>ToString[N[\[Epsilon]],CForm]<>"_a0_Over_rh_"<>ToString[N[ao],CForm];
If[!DirectoryQ[DirName$InputData],CreateDirectory[DirName$InputData]]

For[iQNM=0, iQNM<nQNM/2, iQNM++,

Nz=Nz$dom0;

fn$dom0=DirName$InputData<>"/Cheb_QNM_dom0_N_"<>ToString[Nz]<>"_n"<>ToString[iQNM]<>"_l"<>ToString[l]<>"_spin"<>ToString[spin]<>"_Prec_"<>ToString[N[Prec]]<>".dat";


NumericalData={
"s\t"<>ToString[Re[s$n0[[2*iQNM+1]]]]<>"\t"<>ToString[Im[s$n0[[2*iQNM+1]]]]<>"\n",
"N\t"<>ToString[Nz]<>"\n\n"
};

fp=OpenWrite[fn$dom0];
Table[ WriteString[fp,NumericalData[[i]]], {i,1,Length@NumericalData}];

WriteString[fp,"#1:i\t #2:Real(Cheb_phi_QNM)\t #3:Imag(Cheb_phi_QNM)\n"]
For[i=0, i<=Nz, i++,
WriteString[fp,ToString[i],"\t",
ToString[Re[cheb$\[Phi]QNM$dom0[[2*iQNM+1,i+1]]],CForm],"\t", ToString[Im[cheb$\[Phi]QNM$dom0[[2*iQNM+1,i+1]]],CForm],"\n"]] ;
Close[fp];


Nz=Nz$dom1;
fn$dom1=DirName$InputData<>"/Cheb_QNM_dom1_N_"<>ToString[Nz]<>"_n"<>ToString[iQNM]<>"_l"<>ToString[l]<>"_spin"<>ToString[spin]<>"_Prec_"<>ToString[N[Prec]]<>".dat";


NumericalData={
"s\t"<>ToString[Re[s$n0[[2*iQNM+1]]]]<>"\t"<>ToString[Im[s$n0[[2*iQNM+1]]]]<>"\n",
"N\t"<>ToString[Nz]<>"\n\n"
};


fp=OpenWrite[fn$dom1];
Table[ WriteString[fp,NumericalData[[i]]], {i,1,Length@NumericalData}];

WriteString[fp,"#1:i\t #2:Real(Cheb_phi_QNM)\t #3:Imag(Cheb_phi_QNM)\n"]
For[i=0, i<=Nz, i++,
WriteString[fp,ToString[i],"\t",
ToString[Re[cheb$\[Phi]QNM$dom1[[2*iQNM+1,i+1]]],CForm],"\t", ToString[Im[cheb$\[Phi]QNM$dom1[[2*iQNM+1,i+1]]],CForm],"\n"]] 
Close[fp];

]

