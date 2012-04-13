clear

%% Magnetostatic potential for a uniformly magnetised square with M = [-1,0,0].

max_i = 20;
max_j = max_i;

x_list = linspace(0.01,0.99,max_i);
y_list = linspace(0,1,max_j);


for i=1:max_i
  for j=1:max_j
    x = x_list(i);
    y = y_list(j);

    fun1 = @(v) 1/sqrt(x^2 + (y - v)^2);
    fun2 = @(v) 1/sqrt((x - 1)^2 + (y - v)^2);

    potential(j,i) = (1/(4*pi))*(quad(fun1,0,1) - quad(fun2,0,1));
    end
end

surf(x_list,y_list,potential)


%% Plot from oomph-lib data

load "messing/oomph"
x_oomph=oomph(:,1);
y_oomph=oomph(:,2);
z_oomph=oomph(:,3);

dx=0.05;
dy=0.05;

% Do something so that we can plot the data using surf (copied from the internet)
x_edge=[floor(min(x_oomph)):dx:ceil(max(x_oomph))];
y_edge=[floor(min(y_oomph)):dy:ceil(max(y_oomph))];
[X,Y]=meshgrid(x_edge,y_edge);
Z=griddata(x_oomph,y_oomph,z_oomph,X,Y);

figure
surf(X,Y,Z)


%% get a sort of difference between the two

a = interp2(x_list,y_list,potential,X,Y);
diff = abs(a - Z);

figure
surf(X,Y,diff)