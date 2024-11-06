L = 10;     %lenght
N = 400;    %segments for discretization
T=20;       %max time
dt=0.01;  % time step for disc.

xmesh = linspace(0,L,N);
dx = xmesh(2)-xmesh(1);

% System properties
pde.damping=1e-2;

x_poly = [-0.25 -0.5 -1.0 -4.0 0 0.25 0.5 1.0 4.0];
y_poly = [-0.4 -0.7 -1.1 -8.0 0 0.4 0.7 1.1 8.0];
p_poly = polyfit(x_poly,y_poly,4);
ip_poly=[p_poly./[length(p_poly):-1:1], 0];
clear x_poly y_poly
pde.dEdq= @(p,q) polyval(p_poly,q);      %stress/strain
pde.Eq= @(q) polyval(ip_poly,q);

pde.dEdq= @(p,q) 1./(1+exp(-q*2))-0.5;
pde.Eq= @(q) 0.5*log(exp(2*q)+1)-0.5*q;

pde.dEdp= @(p,q) p/1;            % momentum / velocity
pde.Ep= @(p) 0.5*p.^2;   %H=Ep+Eq    

figure(123)
clf
plot(linspace(-3,3,100),pde.dEdq(0,linspace(-3,3,100)))
hold on
plot(linspace(0,3,100),1*linspace(0,3,100),'--')
xlabel('strain')
ylabel('stress')
%% Conditions
% Init conditions
%xt0 = zeros(length(xmesh),2);
diffgauss=@(x) -4*(x-L/2).*exp(-(x-L/2).^2);
xt0 = [zeros(N,1),diffgauss(xmesh)'];
% Dirichet Boundary conditions

lbf =@(t) 0;
rbf =@(t) 0;%-0.5*(sign(t-2*pi)-1).*sin(t);

%%
[t,x_temp]=ode45(@(t,y) odefun(t,y,N,dx,pde,lbf,rbf),0:dt:T,xt0(:));
x=zeros(N,length(t),2);
x(:,:,1)=x_temp(:,1:N)';
x(:,:,2)=x_temp(:,N+1:end)';
disp('done')
%% over time
pos=cumsum(x(:,:,2))*dx;
% Plot the initial condition.
figure(1)
clf
handle_line = plot(xmesh,pos(:,1),'LineWidth',2);
hold on;
%hold on;
%handle_dot = plot(xmesh(end/2),pos(end/2,1),'o',MarkerSize=10);
axis([0,L,-2,2]);
xlabel('x'); ylabel('u');
title('Wave equation');

for ii=1:5:length(t)
    handle_line.YData = pos(:,ii);
    title(['Wave equation. Time step = ' num2str(ii) 's']);
    drawnow;
end
%% Energy
energy=sum(pde.Ep(x(:,:,1))+pde.Eq(x(:,:,2)))*dx;
figure(2)
plot(t,energy)
%% MODEL
% Get data
dxdt=zeros(2*N,length(t));
for i=1:length(t)
    dxdt(:,i)=odefun(t(i),x_temp(i,:)',N,dx,pde,lbf,rbf);
end
p_data=x(:,:,1);
q_data=x(:,:,2);
dp_data=dxdt(1:N,:);
dq_data=dxdt(N+1:end,:);

%% Integrate and add constant
dp_int = cumtrapz(dp_data)*dx;
dq_int = cumtrapz(dq_data)*dx;
for i=1:length(t)
    [~,closestIndex] = min(abs(q_data(:,i)));
    dp_int(:,i)=dp_int(:,i)-dp_int(closestIndex,i);
end
%% Plot actual dEdq / dEdp and data points
[Z_mesh,T_mesh]=meshgrid(1:20:N,1:10:length(t));
z_vec=Z_mesh(:);
t_vec=T_mesh(:);
real_dEdp=zeros(length(z_vec),1);
model_dEdp=zeros(length(z_vec),1);
real_dEdq=zeros(length(z_vec),1);
model_dEdq=zeros(length(z_vec),1);
p_in=zeros(length(z_vec),1);
q_in=zeros(length(z_vec),1);
%
for i=1:length(z_vec)
    p_in(i)=p_data(z_vec(i),t_vec(i));
    q_in(i)=q_data(z_vec(i),t_vec(i));
    %
    real_dEdp(i)=pde.dEdp(p_in(i),q_in(i));
    model_dEdp(i)=dq_int(z_vec(i),t_vec(i));
    %
    real_dEdq(i)=pde.dEdq(p_in(i),q_in(i));
    model_dEdq(i)=dp_int(z_vec(i),t_vec(i));
end
%
[Z_mesh2, T_mesh2] = meshgrid (1:5:size(q_data, 1), 1:200:size(q_data, 2));
z_vec2 = Z_mesh2(:);
t_vec2 = T_mesh2(:);
q_values = arrayfun(@(z,t) q_data(z, t), z_vec2, t_vec2);
[Z_mesh1, T_mesh1] = meshgrid (1:20:size(q_data, 1), 1:200:size(q_data, 2));
z_vec1 = Z_mesh1(:);
t_vec1 = T_mesh1(:);
downsampled_q_values = arrayfun(@(z,t) q_data(z, t), z_vec1, t_vec1);
p_values = arrayfun(@(z,t) p_data(z, t), z_vec, t_vec);
downsampled_p_values = arrayfun(@(z,t) p_data(z, t), z_vec1, t_vec1);

scalingFactor = 10/400;
figure(10);
clf
time_derivative = diff(downsampled_q_values, 1, 2);
[S, T] = meshgrid(1:size(q_data, 2), 1:size(q_data, 1)); % S for spatial, T for time
% Plotting the original q_data
T = T*scalingFactor;
surf(S, T, q_data,'FaceColor','flat','FaceAlpha',0.7);
shading interp; % This makes the color transitions smooth
hold on; % Keep the current plot
%quiver3(t_vec1, z_vec1, downsampled_q_values, zeros(size(time_derivative(1,:))), time_derivative(1,:), zeros(size(time_derivative(1,:))), 'k');
scatter3(t_vec2, z_vec2*scalingFactor, q_values, 'r'); % 36 is the marker size, 'r' for red color
scatter3(t_vec1, z_vec1*scalingFactor, downsampled_q_values, 50, 'b *');
xlabel('time step t', 'FontSize', 12);
ylabel('spatial z', 'FontSize', 12);
zlabel('state x', 'FontSize', 12);
%title('3D Plot of q Observation and Upsampled Data', 'FontSize', 12);
legend('State Surface','Upsampled x', 'Observed x', 'FontSize', 12);
hold off; % Release the current plot
view(3); % Sets the view to 3D perspective
view(-20, 30);

figure(11);
clf
[S, T] = meshgrid(1:size(p_data, 2), 1:size(p_data, 1)); % S for spatial, T for time
% Plotting the original q_data
surf(S, T, p_data,'FaceColor','r','FaceAlpha',0.6);
shading interp; % This makes the color transitions smooth
hold on; % Keep the current plot
scatter3(t_vec, z_vec, p_values, 'g v'); % 36 is the marker size, 'r' for red color
scatter3(t_vec1, z_vec1, downsampled_p_values, 50, 'k');
xlabel('Time', 'FontSize', 12);
ylabel('Spatial Dimension', 'FontSize', 12);
zlabel('p Value', 'FontSize', 12);
title('3D Plot of p Observation and Upsampled Data', 'FontSize', 12);
legend('p Surface','Upsampled p', 'Observed p', 'FontSize', 12);
hold off; % Release the current plot
view(3); % Sets the view to 3D perspective
view(-20, 30);

figure(1234);
clf
surf(reshape(p_in,size(Z_mesh,1),size(Z_mesh,2)),reshape(q_in,size(Z_mesh,1),size(Z_mesh,2)),reshape(real_dEdp,size(Z_mesh,1),size(Z_mesh,2)));
hold on;
plot3(p_in,q_in,model_dEdp,'g+');
title('dEdp');
figure(1235);
clf
surf(reshape(p_in,size(Z_mesh,1),size(Z_mesh,2)),reshape(q_in,size(Z_mesh,1),size(Z_mesh,2)),reshape(real_dEdq,size(Z_mesh,1),size(Z_mesh,2)));
hold on;
plot3(p_in,q_in,model_dEdq,'g+');
title('dEdq');
%% Train GP
downscale=50;
gpphs.dEdp_model=fitrgp([p_in(1:downscale:end),q_in(1:downscale:end)],model_dEdp(1:downscale:end),'KernelFunction','squaredExponential','FitMethod', ...
    'sr','PredictMethod','sd');
compact_dEdp = compact(gpphs.dEdp_model);
gpphs.dEdq_model=fitrgp([p_in(1:downscale:end),q_in(1:downscale:end)],model_dEdq(1:downscale:end),'KernelFunction','squaredExponential','FitMethod', ...
    'sr','PredictMethod','sd');
compact_dEdq = compact(gpphs.dEdq_model);

gpphs.dEdp = @(p,q) predict(gpphs.dEdp_model,[p,q]);
gpphs.dEdq = @(p,q) predict(gpphs.dEdq_model,[p,q]);


%% Plot model vs actual
[P_mesh,Q_mesh]=meshgrid(-1:0.1:1,-1:0.1:1);
P_vec=P_mesh(:);
Q_vec=Q_mesh(:);
gpdEdp_pred = gpphs.dEdp(P_vec,Q_vec);
gpdEdq_pred = gpphs.dEdq(P_vec,Q_vec);
[dppred,~,dpci] = predict(compact_dEdp,[P_vec,Q_vec]);
[dqpred,~,dqci] = predict(compact_dEdq,[P_vec,Q_vec]);



alphaList = [0.01, 0.05, 0.10]; % Corresponds to 99%, 95%, 90% confidence intervals

% Initialize variables to store the maximum losses
maxLoss_dEdp = 0;
maxLoss_dEdq = 0;

% Loop over each alpha value
% for i = 1:length(alphaList)
%     alpha = alphaList(i);
% 
%     dpdtGP_model=fitrgp([z_vec,t_vec],p_in,'KernelFunction','squaredExponential','FitMethod', 'sr','PredictMethod','sd');
%     compact_dpdtGP = compact(dpdtGP_model);
%     dqdtGP_model=fitrgp([z_vec,t_vec],q_in,'KernelFunction','squaredExponential','FitMethod','sr','PredictMethod','sd');
%     compact_dqdtGP= compact(dqdtGP_model);
% 
% end

% [U_T_mesh, U_Z_mesh] = meshgrid(1:N,1:2001);
% U_T_vec = U_T_mesh(:);
% U_Z_vec = U_Z_mesh(:);
% alpha = 0.05;
% [dpdtpred_sample, ~, dpci_1] = predict(compact_dpdtGP, [U_Z_vec, U_T_vec], 'Alpha', alpha);
% [dqdtpred_sample, ~, dpci_2] = predict(compact_dqdtGP, [U_Z_vec, U_T_vec], 'Alpha', alpha);
% dpdtpred_sample = reshape(dpdtpred_sample, 400, 2001);
% dpci_1_up = reshape(dpci_1(:,1), 400, 2001);
% dpci_1_down = reshape(dpci_1(:,2), 400, 2001);
% dqdtpred_sample = reshape(dqdtpred_sample, 400, 2001);
% dpci_2_up = reshape(dpci_2(:,1), 400, 2001);
% dpci_2_down = reshape(dpci_2(:,2), 400, 2001);


figure(1001)
clf
surf(P_mesh,Q_mesh,reshape(gpdEdp_pred,size(P_mesh)),'FaceColor','r','FaceAlpha',0.5);
hold on
surf(P_mesh,Q_mesh,reshape(dpci(:,2),size(P_mesh)),'FaceColor','g','FaceAlpha',0.5);
hold on
plot3(p_in,q_in,model_dEdp,'g+');
surf(P_mesh,Q_mesh,reshape(pde.dEdp(P_vec,Q_vec),size(P_mesh)),'FaceColor','b','FaceAlpha',0.7);
legend(["GP-DPHS model","GP-DPHS model", "Observation data", "True system"]);
% colorbar; % Adds a color bar to indicate height levels
zlabel('dHdp');
title('3D Slices of dHdp Visualization');
view(3); % Sets the view to 3D perspective
view(70, 20);
hold off; % Release the plot hold

figure(1002)
clf
surf(P_mesh,Q_mesh,reshape(gpdEdq_pred,size(P_mesh)),'FaceColor','r','FaceAlpha',0.3);
hold on
surf(P_mesh,Q_mesh,reshape(dqci(:,1),size(P_mesh)),'FaceColor','g','FaceAlpha',0.5);
hold on
plot3(p_in,q_in,model_dEdq,'g+');
surf(P_mesh,Q_mesh,reshape(pde.dEdq(P_vec,Q_vec),size(P_mesh)),'FaceColor','b','FaceAlpha',0.7);
legend(["GP-DPHS model","GP-DPHS model","Observation data", "True system"]);
% colorbar; % Adds a color bar to indicate height levels
xlabel('state p', 'FontSize', 12);
ylabel('state q', 'FontSize', 12);
zlabel('derivative of Hamiltonian \delta_q H', 'FontSize', 14);
% title('3D  Derivative of Hamiltonian to Q Visualization', 'FontSize', 12);
view(3); % Sets the view to 3D perspective
view(70, 20);
hold off; % Release the plot hold

figure(1003)
clf
surf(P_mesh,Q_mesh,reshape(gpdEdq_pred,size(P_mesh)),'FaceColor','r','FaceAlpha',0.8);
hold on
surf(P_mesh,Q_mesh,reshape(dqci(:,1),size(P_mesh)),'FaceColor','m','FaceAlpha',0.2);
hold on
surf(P_mesh,Q_mesh,reshape(dqci(:,2),size(P_mesh)),'FaceColor','m','FaceAlpha',0.2);
hold on
surf(P_mesh,Q_mesh,reshape(pde.dEdq(P_vec,Q_vec),size(P_mesh)),'FaceColor','b','FaceAlpha',0.6);
hold on
plot3(p_in,q_in,model_dEdq,'g+');
legend(["Mean Prediction","Upper Bound by GP-DPHS", "Lower Bound by GP-DPHS" , "True System" ,"Observation data"]);
% colorbar; % Adds a color bar to indicate height levels
xlabel('state p', 'FontSize', 12);
ylabel('state q', 'FontSize', 12);
zlabel('derivative of Hamiltonian \delta_q H', 'FontSize', 14);
% title('3D  Derivative of Hamiltonian to Q Visualization', 'FontSize', 12);
view(3); % Sets the view to 3D perspective
view(70, 20);
hold off; % Release the plot hold
%% Run GP model vs real
T=40; 
gpphs.damping=pde.damping;
xt0 = [zeros(N,1),4*pi/L*cos(pi/L*xmesh')];
% Dirichet Boundary conditions
lbf =@(t) 0;
rbf =@(t) 0;


[t,x_temp]=ode45(@(t,y) odefun(t,y,N,dx,pde,lbf,rbf),0:dt:T,xt0(:));
x=zeros(N,length(t),2);
x(:,:,1)=x_temp(:,1:N)';
x(:,:,2)=x_temp(:,N+1:end)';
disp('done real')

pos=cumsum(x(:,:,2))*dx;

[t,x_temp]=ode45(@(t,y) odefun(t,y,N,dx,gpphs,lbf,rbf),0:dt:T,xt0(:));
x=zeros(N,length(t),2);
x_gp(:,:,1)=x_temp(:,1:N)';
x_gp(:,:,2)=x_temp(:,N+1:end)';
disp('done gp')

pos_gp=cumsum(x_gp(:,:,2))*dx;

%% Plot 
figure(1)
clf
handle_line = plot(xmesh,pos(:,1),'-','LineWidth',2);
hold on;
handle_line_gp = plot(xmesh,pos_gp(:,1),'LineWidth',2);
legend({'true','model'});
%hold on;
%handle_dot = plot(xmesh(end/2),pos(end/2,1),'o',MarkerSize=10);
axis([0,L,-5,5]);
xlabel('x'); ylabel('u');
title('Wave equation');

for ii=1:10:length(t)
    handle_line.YData = pos(:,ii);
    handle_line_gp.YData = pos_gp(:,ii);
    title(['Wave equation. Time step = ' num2str(ii) 's']);
    drawnow;
end

figure(999);
clf
[Ys, Xs] =  meshgrid(1:400, 1:400);
surf(Ys*5, Xs*scalingFactor, pos(:,1:10:end-1));
% colormap(jet); % You can choose different colormaps like 'jet', 'hsv', etc.
shading interp; % This makes the color transitions smooth
colorbar; % Adds a color bar to indicate height levels
xlabel('time t', 'FontSize', 14);
ylabel('spatial z', 'FontSize', 14);
zlabel('deflection x', 'FontSize', 14);
% title('Wave Landscape (Unified)', 'FontSize', 14);
view(3);
view(30, 20);

figure(997);
clf
xSlices = [120];
% xSlices = [2, 40, 80, 120, 160, 200, 280, 320, 360, 400];
[Ys, Xs] = meshgrid(1:400, 1:400);
hold on; % Keep the plot held so all slices are plotted on the same figure

colors = winter(length(xSlices));
colors_gp = winter(length(xSlices));

% Plot the landscape slices
for i = 1:length(xSlices)
    sliceX = xSlices(i);
    sliceY = Ys(sliceX, :);
    aggpos = pos(:,1:10:end-1);
    aggpos_gp = pos_gp(:,1:10:end-1);

    newpos = aggpos(:, sliceX);
    newpos_gp = squeeze(aggpos_gp(:, sliceX));

    plot(sliceY*scalingFactor, newpos_gp, 'Color', 'r', 'LineWidth', 2);%colors(i,:)
    % plot(sliceY*scalingFactor, newpos_gp_u,'--', 'Color', 'b', 'LineWidth', 2);
    % plot(sliceY*scalingFactor, newpos_gp_u85,'--', 'Color', 'b', 'LineWidth', 2);
    % plot(sliceY*scalingFactor, newpos_gp_l,'--', 'Color', 'b', 'LineWidth', 2);
    % plot(sliceY*scalingFactor, newpos_gp_l85,'--', 'Color', 'b', 'LineWidth', 2);
    % plot(sliceY*scalingFactor, newpos_gp_l60,'--', 'Color', 'b', 'LineWidth', 2);
end
colormap(colors); % Apply the colormap
legend(["Mean Prediction","Sample 1","Sample 2","Sample 3","Sample 4","Sample 5"], 'FontSize', 14)
% colorbar; % Adds a color bar to indicate height levels
xlabel('spatial z', 'FontSize', 14);
ylabel('deflection x', 'FontSize', 14);
% title('3D Slices of wave Visualization', 'FontSize', 14);
hold off; % Release the plot hold

figure(998);
clf
% xSlices = [200];
xSlices = [2, 40, 80, 120, 160, 200, 280, 320, 360, 400];
[Ys, Xs] = meshgrid(1:400, 1:400);
hold on; % Keep the plot held so all slices are plotted on the same figure

colors = winter(length(xSlices));
colors_gp = winter(length(xSlices));

% Plot the landscape slices
for i = 1:length(xSlices)
    sliceX = xSlices(i);
    sliceY = Ys(sliceX, :);
    aggpos = pos(:,1:10:end-1);
    aggpos_gp = pos_gp(:,1:10:end-1);

%     aggpos_gp_u = pos_gp_u(:,1:10:end-1);
%     aggpos_gp_l = pos_gp_l(:,1:10:end-1);

    newpos = aggpos(:, sliceX);
    newpos_gp = aggpos_gp(:, sliceX);

%     newpos_gp_u = aggpos_gp_u(:, sliceX);
%     newpos_gp_l = aggpos_gp_l(:, sliceX);

    plot3(repmat(sliceX, size(sliceY))*5, sliceY*scalingFactor, newpos_gp,'--', 'Color', 'g', 'LineWidth', 2);%colors(i,:)
    plot3(repmat(sliceX, size(sliceY))*5, sliceY*scalingFactor, newpos, 'Color', 'b', 'LineWidth', 2);
%     plot3(repmat(sliceX, size(sliceY))*5, sliceY*scalingFactor, newpos_gp_l,'--', 'Color', 'b', 'LineWidth', 2);
end
colormap(colors); % Apply the colormap
legend(["GP-dPHS model","True System"], 'FontSize', 14)
% colorbar; % Adds a color bar to indicate height levels
xlabel('time step t', 'FontSize', 14);
ylabel('spatial z', 'FontSize', 14);
zlabel('deflection x', 'FontSize', 14);
% title('3D Slices of wave Visualization', 'FontSize', 14);
view(3); % Sets the view to 3D perspective
view(30, 20);
hold off; % Release the plot hold
%% Compare
downscale=200;
gpphs1.dEdp_model=fitrgp([p_in(1:downscale:end),q_in(1:downscale:end)],model_dEdp(1:downscale:end));
size(model_dEdp(1:downscale:end))
gpphs1.dEdq_model=fitrgp([p_in(1:downscale:end),q_in(1:downscale:end)],model_dEdq(1:downscale:end));
gpphs1.dEdp = @(p,q) predict(gpphs1.dEdp_model,[p,q]);
gpphs1.dEdq = @(p,q) predict(gpphs1.dEdq_model,[p,q]);
gpdEdp_pred_l = gpphs1.dEdp(P_vec,Q_vec);
gpdEdq_pred_l = gpphs1.dEdq(P_vec,Q_vec);
gpphs1.damping=pde.damping;
[t,x_temp]=ode45(@(t,y) odefun(t,y,N,dx,gpphs1,lbf,rbf),0:dt:T,xt0(:));
x=zeros(N,length(t),2);
x_gp1(:,:,1)=x_temp(:,1:N)';
x_gp1(:,:,2)=x_temp(:,N+1:end)';
pos_gp1=cumsum(x_gp1(:,:,2))*dx;

downscale=20;
gpphs2.dEdp_model=fitrgp([p_in(1:downscale:end),q_in(1:downscale:end)],model_dEdp(1:downscale:end));
gpphs2.dEdq_model=fitrgp([p_in(1:downscale:end),q_in(1:downscale:end)],model_dEdq(1:downscale:end));
gpphs2.dEdp = @(p,q) predict(gpphs2.dEdp_model,[p,q]);
gpphs2.dEdq = @(p,q) predict(gpphs2.dEdq_model,[p,q]);
gpdEdp_pred_2 = gpphs2.dEdp(P_vec,Q_vec);
gpdEdq_pred_2 = gpphs2.dEdq(P_vec,Q_vec);
gpphs2.damping=pde.damping;
[t,x_temp]=ode45(@(t,y) odefun(t,y,N,dx,gpphs2,lbf,rbf),0:dt:T,xt0(:));
x=zeros(N,length(t),2);
x_gp2(:,:,1)=x_temp(:,1:N)';
x_gp2(:,:,2)=x_temp(:,N+1:end)';
pos_gp2=cumsum(x_gp2(:,:,2))*dx;

downscale=10;
gpphs3.dEdp_model=fitrgp([p_in(1:downscale:end),q_in(1:downscale:end)],model_dEdp(1:downscale:end));
gpphs3.dEdq_model=fitrgp([p_in(1:downscale:end),q_in(1:downscale:end)],model_dEdq(1:downscale:end));
gpphs3.dEdp = @(p,q) predict(gpphs3.dEdp_model,[p,q]);
gpphs3.dEdq = @(p,q) predict(gpphs3.dEdq_model,[p,q]);
gpdEdp_pred_3 = gpphs3.dEdp(P_vec,Q_vec);
gpdEdq_pred_3 = gpphs3.dEdq(P_vec,Q_vec);
gpphs3.damping=pde.damping;
[t,x_temp]=ode45(@(t,y) odefun(t,y,N,dx,gpphs3,lbf,rbf),0:dt:T,xt0(:));
x=zeros(N,length(t),2);
x_gp3(:,:,1)=x_temp(:,1:N)';
x_gp3(:,:,2)=x_temp(:,N+1:end)';
pos_gp3=cumsum(x_gp3(:,:,2))*dx;

DIF1 = pos_gp1 - pos;
DIF2 = pos_gp2 - pos;
DIF3 = pos_gp3 - pos;

mae_1 = zeros(1, 4001);
mae_2 = zeros(1, 4001);
mae_3 = zeros(1, 4001);
for ll=1:length(t)
    mae_1(ll) = mean(abs(DIF1(:, ll)));
    mae_2(ll) = mean(abs(DIF2(:, ll)));
    mae_3(ll) = mean(abs(DIF3(:, ll)));
end