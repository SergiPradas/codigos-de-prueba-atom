clear all; close all; 
%% CREACION DE HISTOGRAMAS DE LAS POSICIONES DE LAS PARTICULAS
%{
%% EXTRACCION DE LOS VALORES Y LECTURA DE NPOP, NP, N1, N2
fileID = fopen('resultados_1.out');
C = textscan(fileID,'%f%f%f%f');
C1 = cell2mat(C(1)); C2 = cell2mat(C(2)); C3 = cell2mat(C(3)); C4 = cell2mat(C(4));
npop = C1(1); np = C2(1); N1 = C3(1); N2 = C4(1);
iter_mass = 31; iter_momento = 10;
num_fig = 1;

%% AHORA TENEMOS QUE HACER UN BUCLE PARA SEPARAR CADA UNA DE LAS ITERACIONES CON DIFERENTES VALORES DE MASA Y MOMENTO MAGNETICO
vec_beta =[]; vec_theta = [];
mat_x = []; %Colocamos en cada columna la coordenada x1 de cada particula
             %Columna diferente implica iteracion diferente
mat_y = []; 
mat_xident = [];

mat_r = []; %Lo mismo que las posiciones, pero ahora con:
            % r = sqrt(x1^2+x2^2)
mat_r_1 = []; %Matriz solo para las de tipo 1
mat_r_2 = []; %Matriz solo para las de tipo 2 
for kk = 1:5
    beta = C1(2 + 201*(kk-1)); theta = C2(2 + 201*(kk-1));
    vec_beta(kk) = beta; vec_theta(kk) = theta;
    
    mat_x(:,kk) = C1((3 + 201*(kk-1)) : (1+201*kk)); 
    mat_y(:,kk) = C2((3 + 201*(kk-1)) : (1+201*kk));
    for jj = 1:length(mat_x(:,kk))
        mat_r(jj,kk) = sqrt((mat_x(jj,kk))^2 +(mat_y(jj,kk))^2);
    end
    %mat_r(:,kk) = (mat_x(:,kk).*mat_x(:,kk) + mat_y(:,kk).*mat_y(:,kk)).^(1/2);
    mat_xident(:,kk) = C3((3 + 201*(kk-1)) : (1+201*kk));  
    % Ahora separamos las de tipo 1 de las de tipo 2:
    ii_r_1 = 1; ii_r_2 = 1;
    for ii = 1: length(mat_x(:,kk))
     if mat_xident(ii,kk) == 1
         mat_r_1(ii_r_1,kk) = mat_r(ii,kk);
         ii_r_1 = ii_r_1 + 1;
     else 
         mat_r_2(ii_r_2,kk) = mat_r(ii,kk);
         ii_r_2 = ii_r_2 + 1;
     end 
    end
    % Y ahora empezamos con el calculo de los histogramas:
    figure(num_fig)
    h1 = histogram(mat_r_1(:,kk),'FaceColor','blue','FaceAlpha', 0.5);
    hold on
    pd1 = fitdist(mat_r_1(:,kk),'Kernel','Support','positive');
    r_values=[0.001:0.001:11];
    %plot(r_values,pdf(pd1,r_values),'-b','LineWidth',2)
    hold on
    %plot(r_values, N1*pdf(pd1,r_values),'-k','LineWidth',2)
    %hold on
    h2 = histogram(mat_r_2(:,kk),'FaceColor','yellow','FaceAlpha',0.5);
    hold on
    pd2 = fitdist(mat_r_2(:,kk),'Kernel','Support','positive');
    %plot(r_values,N2*pdf(pd2,r_values),'-yellow','LineWidth',2)
    
    h1.BinWidth = 0.125;
    h2.BinWidth = 0.125;
    %h1.Normalization = 'countdensity';
    %h2.Normalization = 'countdensity'; 
    %mod = repmat(mat_r_1(:,kk),4,1); %Ho repetim 4 vegades.
    %hprova = histogram(mod,'FaceColor','black');
    %hprova.BinWidth = 0.25;
    %hold on
    %pd1 = fitdist(mod,'Kernel');
    %r_values=[0:0.01:11];
    %plot(r_values,pdf(pd1,r_values),'-b','LineWidth',2)
    
    %legend('Type 1 particles', 'Kernel Fit 1', 'Type 2 particles', 'Kernel Fit 2')
    legend('Type 1 particles', 'Type 2 particles')
    title(['\beta = ',num2str(beta),', \theta = ', num2str(theta)])
    xlabel('r')
 %{   
    num_fig = num_fig + 1;
    vec_center_bin_1 = linspace(h1.BinLimits(1) + h1.BinWidth/2,h1.BinLimits(2)-h1.BinWidth/2,length(h1.Values));
    h3values = (1/(2*pi*h1.BinWidth))*(h1.Values./vec_center_bin_1);
    figure(num_fig)
    plot(vec_center_bin_1,h3values,'b--o','LineWidth',2)
    hold on
    vec_center_bin_2 = linspace(h2.BinLimits(1) + h2.BinWidth/2,h2.BinLimits(2)-h2.BinWidth/2,length(h2.Values));
    h4values = (1/(2*pi*h2.BinWidth))*(h2.Values./vec_center_bin_2);
    plot(vec_center_bin_2,h4values,'y--o','LineWidth',2)
%}    
    
    num_fig = num_fig + 1;
    
    figure(num_fig)
    %plot(r_values,pdf(pd1,r_values),'-b','LineWidth',2)
    plot(r_values,(N1/(2*pi))*pdf(pd1,r_values)./r_values,'-b','LineWidth',2)
    hold on
    %plot(r_values,pdf(pd2,r_values),'-y','LineWidth',2)
    plot(r_values,(N2/(2*pi))*pdf(pd2,r_values)./r_values,'-y','LineWidth',2)
    
    num_fig = num_fig + 1;
   
end
%}
%% AHORA VOY A PROBAR CON EL FICHERO PRUEBA.OUT
N1 = 100; N2 = 100;
fileID = fopen('posiciones.out');
C = textscan(fileID,'%f%f%f');
C1 = cell2mat(C(1)); C2 = cell2mat(C(2)); C3 = cell2mat(C(3));

mat_x = []; %Colocamos en cada columna la coordenada x1 de cada particula
             %Columna diferente implica iteracion diferente
mat_y = []; 
mat_xident = [];

mat_r = []; %Lo mismo que las posiciones, pero ahora con:
            % r = sqrt(x1^2+x2^2)
mat_r_1 = []; %Matriz solo para las de tipo 1
mat_r_2 = []; %Matriz solo para las de tipo 2 
num_fig = 1;
for kk = 1:50
    mat_x(:,kk) = C1((1 + 200*(kk-1)) : (200*kk)); 
    mat_y(:,kk) = C2((1 + 200*(kk-1)) : (200*kk));
    for jj = 1:length(mat_x(:,kk))
        mat_r(jj,kk) = sqrt((mat_x(jj,kk))^2 +(mat_y(jj,kk))^2);
    end
    %mat_r(:,kk) = (mat_x(:,kk).*mat_x(:,kk) + mat_y(:,kk).*mat_y(:,kk)).^(1/2);
    mat_xident(:,kk) = C3((1 + 200*(kk-1)) : (200*kk));  
    % Ahora separamos las de tipo 1 de las de tipo 2:
    ii_r_1 = 1; ii_r_2 = 1;
    for ii = 1: length(mat_x(:,kk))
     if mat_xident(ii,kk) == 1
         mat_r_1(ii_r_1,kk) = mat_r(ii,kk);
         ii_r_1 = ii_r_1 + 1;
     else 
         mat_r_2(ii_r_2,kk) = mat_r(ii,kk);
         ii_r_2 = ii_r_2 + 1;
     end 
    end
    % Y ahora empezamos con el calculo de los histogramas:
%    figure(num_fig)
%    h1 = histogram(mat_r_1(:,kk),'FaceColor','blue','FaceAlpha', 0.5);
%    hold on
    pd1 = fitdist(mat_r_1(:,kk),'Kernel','Support','positive');
    r_values=[0.001:0.001:20];
    mat_pd1(:,kk) = (N1/(2*pi))*pdf(pd1,r_values)./r_values;
    %plot(r_values,pdf(pd1,r_values),'-b','LineWidth',2)
%    hold on
    %plot(r_values, N1*pdf(pd1,r_values),'-k','LineWidth',2)
    %hold on
%    h2 = histogram(mat_r_2(:,kk),'FaceColor','yellow','FaceAlpha',0.5);
%    hold on
    pd2 = fitdist(mat_r_2(:,kk),'Kernel','Support','positive');
    mat_pd2(:,kk) = (N2/(2*pi))*pdf(pd2,r_values)./r_values;
    %plot(r_values,N2*pdf(pd2,r_values),'-yellow','LineWidth',2)
    
%    h1.BinWidth = 0.125;
%    h2.BinWidth = 0.125;
    %h1.Normalization = 'countdensity';
    %h2.Normalization = 'countdensity'; 
    %mod = repmat(mat_r_1(:,kk),4,1); %Ho repetim 4 vegades.
    %hprova = histogram(mod,'FaceColor','black');
    %hprova.BinWidth = 0.25;
    %hold on
    %pd1 = fitdist(mod,'Kernel');
    %r_values=[0:0.01:11];
    %plot(r_values,pdf(pd1,r_values),'-b','LineWidth',2)
    
    %legend('Type 1 particles', 'Kernel Fit 1', 'Type 2 particles', 'Kernel Fit 2')
%    legend('Type 1 particles', 'Type 2 particles')
%    title('m_2/m_1 = 4, \mu_2/\mu_1 = 2')
%    xlabel('r')
 %{   
    num_fig = num_fig + 1;
    vec_center_bin_1 = linspace(h1.BinLimits(1) + h1.BinWidth/2,h1.BinLimits(2)-h1.BinWidth/2,length(h1.Values));
    h3values = (1/(2*pi*h1.BinWidth))*(h1.Values./vec_center_bin_1);
    figure(num_fig)
    plot(vec_center_bin_1,h3values,'b--o','LineWidth',2)
    hold on
    vec_center_bin_2 = linspace(h2.BinLimits(1) + h2.BinWidth/2,h2.BinLimits(2)-h2.BinWidth/2,length(h2.Values));
    h4values = (1/(2*pi*h2.BinWidth))*(h2.Values./vec_center_bin_2);
    plot(vec_center_bin_2,h4values,'y--o','LineWidth',2)
%}    
    
%    num_fig = num_fig + 1;
    
%    figure(num_fig)
%    %plot(r_values,pdf(pd1,r_values),'-b','LineWidth',2)
%    plot(r_values,(N1/(2*pi))*pdf(pd1,r_values)./r_values,'-b','LineWidth',2)
%    hold on
%    %plot(r_values,pdf(pd2,r_values),'-y','LineWidth',2)
%    plot(r_values,(N2/(2*pi))*pdf(pd2,r_values)./r_values,'-y','LineWidth',2)
    
%    num_fig = num_fig + 1;
   
end
%{
for kk = 1:20
figure(num_fig)
plot(r_values(500:end),mat_pd1(500:end,kk),'-b','LineWidth',2)
hold on
plot(r_values(500:end),mat_pd2(500:end,kk),'-y','LineWidth',2)
hold on
end
legend('Type 1 particles', 'Type 2 particles')
title('m_2/m_1 = 4, \mu_2/\mu_1 = 2')    
xlabel('r')
ylabel('\rho(r)')
%}

av_pd1= zeros(length(mat_pd1(:,10)),1); 
av_pd2 = zeros(length(mat_pd1(:,10)),1);
for ii = 1:length(mat_pd1(:,10))
    for jj = 1:50
      av_pd1(ii) = av_pd1(ii) + mat_pd1(ii,jj);
      av_pd2(ii) = av_pd2(ii) + mat_pd2(ii,jj);
    end
    av_pd1(ii) = av_pd1(ii)/length(mat_pd1(ii,:));
    av_pd2(ii) = av_pd2(ii)/length(mat_pd2(ii,:));
end
%plot(r_values(500:end),av_pd1(500:end),'-k','LineWidth',2)
%hold on
%plot(r_values(500:end),av_pd2(500:end),'--k','LineWidth',2)
%legend('Type 1 particles', 'Type 2 particles')
%title('m_2/m_1 = 4, \mu_2/\mu_1 = 2')    
%xlabel('r')
%ylabel('\rho(r)')


%num_fig = num_fig + 1;

figure(num_fig)
subplot(1,2,1)
plot(r_values(1:end),av_pd1(1:end),'-k','LineWidth',2)
hold on
plot(r_values(1:end),av_pd2(1:end),'--k','LineWidth',2)
legend('Type 1 particles', 'Type 2 particles')    
xlabel('r')
ylabel('\rho(r)')
subplot(1,2,2)
semilogx(r_values(1:end),av_pd1(1:end),'-k','LineWidth',2)
hold on
semilogx(r_values(1:end),av_pd2(1:end),'--k','LineWidth',2)
legend('Type 1 particles', 'Type 2 particles')   
xlabel('r')
ylabel('\rho(r)')

sgtitle('Average density profile. m_2/m_1 = 2, \mu_2/\mu_1 = 4')

num_fig = num_fig + 1;
%% Prueba distribucion posiciones
r1t1 = 0; r1t2 = 0;
iter = 5;
for ii = 1: length(mat_x(:,iter))
   if (mat_xident(ii,iter) == 1)
     figure(num_fig)
     plot(mat_x(ii,iter),mat_y(ii,iter),'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b') 
     hold on
     r1t1 = r1t1 + sqrt(mat_x(ii,iter)*mat_x(ii,iter)+mat_y(ii,iter)*mat_y(ii,iter));
  else 
     plot(mat_x(ii,iter),mat_y(ii,iter),'Marker','o','MarkerEdgeColor','m','MarkerFaceColor','m')
     hold on
     r1t2 = r1t2 + sqrt(mat_x(ii,iter)*mat_x(ii,iter)+mat_y(ii,iter)*mat_y(ii,iter));
  end
  xlabel('x1')
  ylabel('x2')
  title('Posicions de les partícules')
  %ylim([-15,15])
  %xlim([-15,15])
end
r1t1 = r1t1/N1
r1t2 = r1t2/N2

num_fig = num_fig + 1;

%{
%% 2D FITTING OF THE KERNEL FUNCTION 
x_values = [-15:0.01:15];
y_values = [-15:0.01:15];
%I els average els guardarem a:
av_pdxy1 = zeros(size(meshgrid(x_values,y_values)));
av_pdxy2 = zeros(size(meshgrid(x_values,y_values)));
%I ara podem fer els calculs:
for kk = 1: 50
    ii_1 = 1; ii_2 = 1;
    for ii = 1: length(mat_x(:,kk))
     if mat_xident(ii,kk) == 1
         mat_x1(ii_1,kk) = mat_x(ii,kk);
         mat_y1(ii_1,kk) = mat_y(ii,kk);
         ii_1 = ii_1 + 1;
     else 
         mat_x2(ii_2,kk) = mat_x(ii,kk);
         mat_y2(ii_2,kk) = mat_y(ii,kk);
         ii_2 = ii_2 + 1;
     end 
    end
    pdx1 = fitdist(mat_x1(:,kk),'Kernel');
    mat_pdx1(:,kk) = pdf(pdx1,x_values);
    pdy1 = fitdist(mat_y1(:,kk),'Kernel');
    mat_pdy1(:,kk) = pdf(pdy1,y_values);
    
    pdx2 = fitdist(mat_x2(:,kk),'Kernel');
    mat_pdx2(:,kk) = pdf(pdx2,x_values);
    pdy2 = fitdist(mat_y2(:,kk),'Kernel');
    mat_pdy2(:,kk) = pdf(pdy2,x_values);
    
    [xx1,yy1]     = meshgrid(x_values,y_values);
    [pdxx1,pdyy1] = meshgrid(mat_pdx1(:,kk),mat_pdy1(:,kk));
    
    [xx2,yy2]     = meshgrid(x_values,y_values);
    [pdxx2,pdyy2] = meshgrid(mat_pdx2(:,kk),mat_pdy2(:,kk));
 
    %Calculem la pdf total asumint que la x i la y son independents:
    pdxy1 = N1*pdxx1.*pdyy1;
    pdxy2 = N2*pdxx2.*pdyy2;

%{    
    %I amb això ja ho podem representar. Abans de fer-ho fixarem els limits
    %de la colorbar, per a que sigui compartida entre els dos tipus de
    %particules:
    bottom = min(min(min(pdxy1)),min(min(pdxy2)));
    top  = max(max(max(pdxy1)),max(max(pdxy2)));
    
    %I ara sí que podem fer la representació:
    colormap(winter)
    figure(num_fig)
    s1 = surf(xx1,yy1,pdxy1);
    set(s1,'LineStyle','none');
    xlabel('x')
    ylabel('y')
    title('\rho_1(x,y)')
    ylim([y_values(1),y_values(end)])
    xlim([x_values(1),x_values(end)])
    caxis manual
    caxis([bottom top]);
    colorbar
    view(2)
    
    
    num_fig = num_fig + 1;
    
    
    figure(num_fig)
    s2 = surf(xx2,yy2,pdxy2);
    set(s2,'LineStyle','none');
    xlabel('x')
    ylabel('y')
    title('\rho_2(x,y)')
    ylim([y_values(1),y_values(end)])
    xlim([x_values(1),x_values(end)])
    caxis manual
    caxis([bottom top]);
    colorbar
    view(2)

    
    num_fig = num_fig + 1;
    
%}
    %Aprofitem per anar aconseguint el average:
    av_pdxy1 = av_pdxy1 + pdxy1;
    av_pdxy2 = av_pdxy2 + pdxy2;
        
end

%I per acabar d'obtenir el average result:
av_pdxy1 = av_pdxy1/kk;
av_pdxy2 = av_pdxy2/kk;

%I llavors ja podem fer el plot:

    %bottom = min(min(min(av_pdxy1)),min(min(av_pdxy2)));
    %top  = max(max(max(av_pdxy1)),max(max(av_pdxy2)));
    
    %I ara sí que podem fer la representació:
    figure(num_fig)
    subplot(2,1,1)
    s1 = surf(xx1,yy1,av_pdxy1);
    set(s1,'LineStyle','none');
    %xlabel('x')
    ylabel('y')
    title('Average \rho_1(x,y)')
    ylim([y_values(1),y_values(end)])
    xlim([x_values(1),x_values(end)])
    %caxis manual
    %caxis([bottom top]);
    colormap(flipud(hot))
    colorbar
    view(2)
  

    subplot(2,1,2)
    s2 = surf(xx2,yy2,av_pdxy2);
    set(s2,'LineStyle','none');
    xlabel('x')
    ylabel('y')
    title('Average \rho_2(x,y)')
    ylim([y_values(1),y_values(end)])
    xlim([x_values(1),x_values(end)])
    %caxis manual
    %caxis([bottom top]);
    colormap(flipud(hot))
    colorbar
    view(2)

    
    num_fig = num_fig + 1;
    
%}

