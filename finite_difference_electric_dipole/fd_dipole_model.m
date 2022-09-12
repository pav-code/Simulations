clc
close all; clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finite Difference Solution of Laplace Equation

%%    4BC3 - Modelling of Biomedical Systems
%%             Pavel Gueorguiev                   %%
%%            McMaster University                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Modifiable variables                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nR = 200;
nTheta = 8;

dR = 1;
dTheta = 2*pi/nTheta;
dThOff = pi/8;

q = 1*10^-9;
eps = 8.854*10^-12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Matrix Dimensions/Definitions         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nPoints = (nR-2)*nTheta;
z = zeros(nPoints,nPoints);
b = zeros(nPoints,1);
V = zeros(nTheta,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Calculating Boundary Condition         %%
%%     Voltage at a close distance(inner ring)    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = dR / 100;
pQ = [0 d];
nQ = [0 -d];
for (i = 1:nTheta)
  theta = i*dTheta + dThOff;
  [x,y] = pol2cart(theta,dR);
  r1 = ((pQ(1) - x)^2 + (pQ(2) - y)^2)^0.5;
  r2 = ((nQ(1) - x)^2 + (nQ(2) - y)^2)^0.5;
  V(i) = q/(4*pi*eps)*(1/r1 - 1/r2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Matrix set-up loop               %%
%%         Solution of Laplace Equation           %%
%%              using Finite Difference           %%
%%    Dirichlet BC : Inner ring (see above)       %%
%%    Neumann BC   : Outermost ring equals inner  %%
%%                   dV/dR = 0                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eqCount = 1;
for i = 2:nR-1
  for j = 1:nTheta
    A = (1/(dR^2) + 1/(i*(dR^2)));
    B = (-2/(dR^2) - 2/((i*dR)^2*(dTheta)^2));
    C = (-1/(i*(dR^2)) + 1/(dR^2));
    D = (1/((i*dR)^2*(dTheta)^2) + cot((j*dTheta)+dThOff)/(2*(dTheta)*(i*dR)^2));
    E = (1/((i*dR)^2*(dTheta)^2) - cot((j*dTheta)+dThOff)/(2*(dTheta)*(i*dR)^2));
	
	if (i == 2)			                        %BC 1) Dirichlet Boundary	
      b(eqCount) = -1*C*V(j);
	end
	
    if (j == 1)				                    %Case: V(#,1)->V(#,0) = V(#,N)
      %C->(theta-2)-> 0BD ->(theta-2)-> A  *Caution* E (shift +Theta) [C 0BD EA]
      z(eqCount,eqCount) = B;
      z(eqCount,eqCount+1) = D;
      z(eqCount,(eqCount-1)+nTheta) = E;
	  if (i != 2)                               %BC 1) Dirichlet Boundary
        z(eqCount,eqCount-nTheta) = C;
	  end
	  if (i != (nR-1))
        z(eqCount,eqCount+nTheta) = A;
      else
        %A = (1/(dR^2) + 1/((i+1)*(dR^2)));	  
        z(eqCount,eqCount) += A;                %BC 2) Neumann Boundary 
      end	  
	elseif (j == nTheta)			            %Case: V(#,N)->V(#,N+1) = V(#,1)
	  %C->(theta-2)-> EB0 ->(theta-2)-> A *Caution* D (shift -Theta) [CD EB0 A]
	    z(eqCount,eqCount) = B;
	    z(eqCount,(eqCount+1)-nTheta) = D;
	    z(eqCount,eqCount-1) = E;
		if (i != 2)		                        %BC 1) Dirichlet Boundary
	      z(eqCount,eqCount-nTheta) = C;
		end
		if (i != (nR-1))                        %BC 2) Neumann Boundary 
          z(eqCount,eqCount+nTheta) = A;
        else
		  %A = (1/(dR^2) + 1/((i+1)*(dR^2)));
		  z(eqCount,eqCount) += A;
        end		
	else
	  %C->(theta-2)-> EBD ->(theta-2)-> A [C EBD A]
	  z(eqCount,eqCount) = B;
	  z(eqCount,eqCount+1) = D;
	  z(eqCount,eqCount-1) = E;
	  if (i != 2)		                        %BC 1) Dirichlet Boundary
	    z(eqCount,eqCount-nTheta) = C;
	  end
	  if (i != (nR-1))
        z(eqCount,eqCount+nTheta) = A;          %BC 2) Neumann Boundary 
      else
	    %A = (1/(dR^2) + 1/((i+1)*(dR^2)));
	    z(eqCount,eqCount) += A;
      end		
	end
	eqCount = eqCount + 1;
  end
end

Vmap = inv(z)*b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Remapping of Results             %%
%%       Adding initial boundary voltage       %%
%% Convert polar to cylindrical and prepare    %%
%% the matrix for display                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vmap = [V;Vmap];
count=1; Rcount=1; 
Vcontour = zeros((nR-1),nTheta);
theta = zeros((nR-1)*nTheta,1);
r = zeros((nR-1)*nTheta,1);

for i = 1:nR-1
  for j = 1:nTheta
    theta(count) = dTheta*count;
	r(count) = dR*Rcount;
    Vcontour(i,j) = Vmap(count);	
	count=count+1;
  end
  Rcount=Rcount+1;
end

theta = theta + dThOff;
[x, y] = pol2cart (theta, r);
count = 1;
for i = 1:nR-1
  for j = 1:nTheta
    xMap(i,j) = x(count);
	yMap(i,j) = y(count);
	count+=1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Graph - Using Trimap             %%
%%           See Trisuft Reference             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

center = zeros(nTheta,1)';
Vcontour = [V'; Vcontour];
xMap = [center; xMap];
yMap = [center; yMap];

clf
tri = delaunay (xMap, yMap);
h = trisurf (tri, xMap, yMap, Vcontour, "facecolor", "interp");
%axis tight
axis([-10 10 -10 10])
zlim auto

title (sprintf("Numerical Solution of Laplace Equation. #Radius = %i, #Theta = %i, dR = %i, dTheta = %i",nR,nTheta,dR,dTheta))
xlabel('Distance in X');
ylabel('Distance in Y');
zlabel('Voltage');
 
%Trisurf Code Reference: http://octave.sourceforge.net/octave/function/trisurf.html
