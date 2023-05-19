clc
clear
close all

fd = @(p) drectangle(p,-1,1,-1,1);
[p,t] = distmesh2d(fd,@huniform,0.04,[-2,-2;2,2],[]);
title("Triangulation")
axis on

coordinates = p;

%%
%load elements3.dat; elements3(:,1)=[];
elements3 = t;
%% Dirichlet
%load dirichlet.dat; dirichlet(:,1) = [];
eps = 1e-5;
ind_l = find(abs(p(:,1)+1)<eps);
lc1 = length(ind_l);
ind_r = find(abs(p(:,1)-1)<eps);
lc2 = length(ind_r);
ind_u = find(abs(p(:,2)+1)<eps);
lc3 = length(ind_u);
ind_d = find(abs(p(:,2)-1)<eps);
lc4 = length(ind_d);


dirichlet = [ind_r;ind_l;ind_u;ind_d];
%%
% figure
% scatter(coordinates(unique(dirichlet),1),coordinates(unique(dirichlet),2))

%%
FreeNodes=setdiff(1:size(coordinates,1),unique(dirichlet));

% Initial value
%U = -ones(size(coordinates,1),1);
U = sign(coordinates(:,1));
%U = 4*rand(size(coordinates,1),1);
%U = U*0;
U(unique(dirichlet)) = u_d(coordinates(unique(dirichlet),:));
%U(unique(dirichlet)) = u_d2(coordinates(unique(dirichlet),:));
show(elements3,[],coordinates,full(U));
title("random initial guess")
%%
% Newton-Raphson iteration
for i=1:1000
  
  % Assembly of DJ(U)
  A = sparse(size(coordinates,1),size(coordinates,1));
  for j = 1:size(elements3,1)
    A(elements3(j,:),elements3(j,:)) = A(elements3(j,:),elements3(j,:)) ...
	+ localdj(coordinates(elements3(j,:),:),U(elements3(j,:)));
  end
  
  % Assembly of J(U)
  b = sparse(size(coordinates,1),1);
  for j = 1:size(elements3,1)
    b(elements3(j,:)) = b(elements3(j,:)) ...
   	+ localj(coordinates(elements3(j,:),:),U(elements3(j,:)));
  end
  
  % Volume Forces
  for j = 1:size(elements3,1)
    b(elements3(j,:)) = b(elements3(j,:)) + ...
	det([1 1 1; coordinates(elements3(j,:),:)']) * ...
	f(sum(coordinates(elements3(j,:),:))/3)/6;
  end
  

  
  % Dirichlet conditions
  W = zeros(size(coordinates,1),1);
  W(unique(dirichlet)) = 0;
  
  % Solving one Newton step
  W(FreeNodes) = A(FreeNodes,FreeNodes)\b(FreeNodes);
  U = U - W;
  if norm(W) < 10^(-10)
    disp(i)
    break
  end
end

% graphic representation
show(elements3,[],coordinates,full(U));

%%
%figure
%plot(U(unique(dirichlet)))
%%
function DirichletBoundaryValue = u_d2(x)
DirichletBoundaryValue =  zeros(size(x,1),1);
L = size(x,1);
eps = 1e-5;
for i = 1:L
    if (abs(x(i,1)-1)<eps)| (abs(x(i,1)+1)<eps)
        DirichletBoundaryValue(i) = 1;
    end
    if (abs(x(i,2)-1)<eps)| (abs(x(i,2)+1)<eps)
        DirichletBoundaryValue(i) = -1;
    end
end
end