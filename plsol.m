close all;
x = load('x.txt');
y = load('y.txt');
[X,Y] = meshgrid(x,y);
U = zeros(length(x));

mask = zeros(length(x));
a = dir("Node_Info*");
b = dir("mask*");
utmp = load('out.txt');
n_nodes = length(a);

for i = 1:n_nodes
	data = load(strcat(a(i).folder,'/',a(i).name))
	n_proc = size(data,1);
	mask_loc = load(strcat(b(i).folder,'/',b(i).name));
	ofst = 0;
	for j = 1:n_proc
		nx = data(j,1);
		ny = data(j,2);
		len = nx*ny;
		xstart = data(j,3);
		ystart = data(j,4);
		mask(xstart:xstart+nx-1,ystart:ystart+ny-1) = reshape(mask_loc(ofst+1:ofst+len),ny,nx)';
		U(xstart:xstart+nx-1,ystart:ystart+ny-1) = reshape(utmp(ofst+1:ofst+len),ny,nx)';
		ofst = ofst + len;
	end
end

figure(1)
% plot(X,Y,'k*');
% hold on
I = find(mask == 1);
plot(X(I),Y(I),'b*');
hold on
I = find(mask == -1);
plot(X(I),Y(I),'r*');
I = find(mask == 0);
plot(X(I),Y(I),'k*');
I = find(mask == -2);
plot(X(I),Y(I),'c*');
I = find(mask == -3);
plot(X(I),Y(I),'g*');

legend('Interior', 'Int. Ghost', 'Exterior','Ext. Ghost','Boundary','Location','NorthEastOutside')
% theta = linspace(0,2*pi,1000);
% plot(0.8*cos(theta),0.8*sin(theta),'m--')

% Plot the solution
figure(4)
surf(X,Y,((mask)==1).*U)

xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'fontsize',21)