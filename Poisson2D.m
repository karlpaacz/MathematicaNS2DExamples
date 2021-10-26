(* ::Package:: *)

(* ::Section:: *)
(*Solution of a Poisson equation*)


(* ::Input:: *)
(*f[x_,y_]:=Exp[-x^2-10y^4-Sin[x y]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*ff=Plot3D[f[x,y],{x,-1,1},{y,-1,1},PlotRange->All,Mesh->None,ColorFunction->"Rainbow",PlotLegends->Automatic,AxesStyle->Directive[Black,Thick],AxesLabel->{"x","y","z"},Boxed->False]*)


(* ::Input:: *)
(*Export["G:\\Unidades compartidas\\Doctorado\\Semestre1\\Proyect-ComputationalAplicationEngeeneringEducation\\Imagenes\\fpoisson.pdf",ff]*)


(* ::Input:: *)
(*uval=NDSolveValue[{D[u[x,y],x,x]+D[u[x,y],y,y]==f[x,y]+NeumannValue[Sin[y],x==-1]+NeumannValue[0,x==1],u[x,-1]==0,u[x,1]==0},u,{x,-1,1},{y,-1,1}]*)


(* ::Input:: *)
(*sNDSolve=Plot3D[uval[x,y],{x,-1,1},{y,-1,1},PlotRange->All,PlotRange->All,Mesh->None,ColorFunction->"Rainbow",PlotLegends->Automatic,AxesStyle->Directive[Black,Thick],AxesLabel->{"x","y","z"},Boxed->False];*)


(* ::Input:: *)
(*Export["G:\\Unidades compartidas\\Doctorado\\Semestre1\\Proyect-ComputationalAplicationEngeeneringEducation\\Imagenes\\sNDSolve.pdf",sNDSolve]*)


(* ::Text:: *)
(*El problema es:*)


(* ::Text:: *)
(*Laplacian[u[x,y],{x,y}]==f[x,y]*)


(* ::Text:: *)
(*Condiciones de frontera :*)


(* ::Text:: *)
(*u[x,-1]==u[x,1]==0;*)
(*\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(x\)]\((u[\(-1\), y])\)\)==Sin[y]*)
(*\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(x\)]\((u[1, y])\)\)==0*)


(* ::Section:: *)
(*Soluci\[OAcute]n Cheb*)


(* ::Subsection:: *)
(*packages*)


(* ::Input:: *)
(*GeneralCollocationPoints[nx_,r0_,rl_]:=Module[{a,b},a=(rl+r0)*0.5;b=(rl-r0)*0.5;*)
(*a+b*Cos[(\[Pi] Range[0,nx])/nx]];*)


(* ::Input:: *)
(*CollocationMatrix[x_]:=-NDSolve`FiniteDifferenceDerivative[1,-x,DifferenceOrder->"Pseudospectral"]["DifferentiationMatrix"];*)
(**)


(* ::Input:: *)
(*fbcRR[d2x_,d1x_,{a0_,b0_,an_,bn_}]:=Module[{i0=1, in=Length[d1x],nx=Length[d1x]-1,c00,c0n,cnn,cn0,ex,bj0,bjn,bcLx,bf0,bfn},*)
(*c00=a0+b0*d1x[[i0,i0]];*)
(*cnn=an+bn*d1x[[in,in]];*)
(*c0n=b0*d1x[[i0,in]];*)
(*cn0=bn*d1x[[in,i0]];*)
(**)
(*ex=c00*cnn-c0n*cn0;*)
(**)
(*bj0=-b0*cnn*d1x[[i0,2;;nx]]+bn*c0n*d1x[[in,2;;nx]];*)
(*bjn=b0*cn0*d1x[[i0,2;;nx]]-bn*c00*d1x[[in,2;;nx]];*)
(**)
(**)
(*bcLx=(KroneckerProduct[d2x[[2;;nx,i0]],bj0]+KroneckerProduct[d2x[[2;;nx,in]],bjn])/ex;*)
(**)
(*bf0=-(cnn*d2x[[2;;nx,i0]]-cn0*d2x[[2;;nx,in]])/ex;*)
(*bfn=-(-c0n*d2x[[2;;nx,i0]]+c00*d2x[[2;;nx,in]])/ex;*)
(*{bcLx,{bf0,bfn},{bj0,bjn,c00,cnn,c0n,cn0}/ex}*)
(*]*)


(* ::Input:: *)
(*fbch[{bf0_,bfn_},{g0_,gn_}]:=TensorProduct[bf0,g0]+TensorProduct[bfn,gn];*)
(*Diagonalization[d2x_,bcux_,d2y_,bcuy_,sgu_]:=Module[{Lx,Ly,p,pi,q,qi,dg,dgx,dgy},*)
(*Lx=d2x+bcux;*)
(*Ly=d2y+bcuy;*)
(*p=Transpose[Eigenvectors[Lx]];*)
(*pi=Inverse[p];*)
(*dgx=Eigenvalues[Lx];*)
(*q=Eigenvectors[Ly];*)
(*qi=Inverse[q];*)
(*dgy=Eigenvalues[Ly];*)
(*dg=Outer[Plus[#1,#2]-sgu&,dgx,dgy];*)
(*{p,pi,q,qi,dg}*)
(*];*)
(*DiagSolve[{p_,pi_,q_,qi_,dg_},h_]:=p . (((pi . h) . qi/dg) . q);*)
(**)
(**)


(* ::Subsection:: *)
(*inputs*)


(* ::Input:: *)
(*nx=30; ny=32;*)
(*x0=-1; xn=1;*)
(*y0=-1; yn=1;*)
(*f[x_,y_]:=Exp[-x^2-10y^4-Sin[x y]]*)


(* ::Input:: *)
(*fx0[x_]:=0;*)
(*fxn[x_]:=Sin[x];*)
(*fy0[x_]:=0;*)
(*fyn[x_]:=0;*)


(* ::Subsection:: *)
(*definitions*)


(* ::Input:: *)
(**)
(*x=GeneralCollocationPoints[nx,x0,xn];*)
(*y=GeneralCollocationPoints[ny,y0,yn];*)
(*gx0=fx0/@y;*)
(*gxn=fxn/@y;*)
(*gy0=fy0/@x;*)
(*gyn=fyn/@x;*)


(* ::Input:: *)
(*u0=Outer[0.#&,x,y];*)
(*u0[[1,All]]=gx0;*)
(*u0[[nx+1,All]]=gxn;*)
(*u0[[All,1]]=gy0;*)
(*u0[[All,ny+1]]=gyn;*)
(*d1x=CollocationMatrix[x];*)
(*d1y=CollocationMatrix[y];*)
(*d2x=d1x . d1x;*)
(*d2y=d1y . d1y;*)
(*idx=IdentityMatrix[nx-1];*)
(*idy=IdentityMatrix[ny-1];*)


(* ::Input:: *)
(**)
(*{ubcLx,{ubfx0,ubfxn},{ubj0x,ubjnx,uc00x,ucnnx,uc0nx,ucn0x}}=fbcRR[d2x,d1x,{0,1,0,1}];*)
(*ubchx=fbch[{ubfx0,ubfxn},{gx0[[2;;ny]],gxn[[2;;ny]]}];*)
(*{ubcLy,{ubfy0,ubfyn},{ubj0y,ubjny,uc00y,ucnny,uc0ny,ucn0y}}=fbcRR[d2y,d1y,{1,0,1,0}];*)
(*ubchy=Transpose[fbch[{ubfy0,ubfyn},{gy0[[2;;nx]],gyn[[2;;nx]]}]];*)


(* ::Input:: *)
(*fi=Outer[f[#1,#2]&,x,y];*)


(* ::Subsection:: *)
(*solve*)


(* ::Input:: *)
(*{p,pi,q,qi,dg}=Diagonalization[d2x[[2;;nx,2;;nx]],ubcLx,d2y[[2;;ny,2;;ny]],ubcLy,0];*)
(*h=fi[[2;;nx,2;;ny]]+ubchx+ubchy;*)


(* ::Input:: *)
(*usol=DiagSolve[{p,pi,q,qi,dg},h];*)


(* ::Input:: *)
(*u0[[2;;nx,2;;ny]]=usol;*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*data=Flatten[Table[{x[[i]],y[[j]],u0[[i,j]]},{i,2,nx},{j,2,ny}],1];*)


(* ::Input:: *)
(*ss=ListPlot3D[data,Mesh->False,ColorFunction->Hue]*)


(* ::Input:: *)
(*Export["G:\\Unidades compartidas\\Doctorado\\Semestre1\\Proyect-ComputationalAplicationEngeeneringEducation\\Imagenes\\solChev.pdf",ss]*)


(* ::Input:: *)
(*s2=Plot3D[uval[x,y],{x,-1,1},{y,-1,1},PlotRange->All,PlotRange->All,Mesh->None,PlotStyle->Blue];*)


(* ::Input:: *)
(*comp=Show[ss,s2]*)


(* ::Input:: *)
(*Export["G:\\Unidades compartidas\\Doctorado\\Semestre1\\Proyect-ComputationalAplicationEngeeneringEducation\\Imagenes\\comp.pdf",comp]*)
