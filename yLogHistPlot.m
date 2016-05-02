function [ph1,ph2] = yLogHistPlot(list1, title1, list2,title2, x)
if nargin<3
hist(list1,150)
ph1 = get(gca,'children');
% Determine number of histogram patches
N_patches = length(ph1);
for i = 1:N_patches
      % Get patch vertices
      vn = get(ph1(i),'Vertices');
      % Adjust y location
      vn(:,2) = vn(:,2) + 1;
      % Reset data
      set(ph1(i),'Vertices',vn)
end
set(gca,'yscale','log')
title(title1)    
grid on
grid minor
end

if nargin>3
    if nargin==4
        xlimit=[min(min(list1),min(list2)),max(max(list1),max(list2))];
    end
    if nargin == 5
        xlimit=x;
    end
subplot(2,1,1)
hist(list1,150)
ph1 = get(gca,'children');
% Determine number of histogram patches
N_patches = length(ph1);
for i = 1:N_patches
      % Get patch vertices
      vn = get(ph1(i),'Vertices');
      % Adjust y location
      vn(:,2) = vn(:,2) + 1;
      % Reset data
      set(ph1(i),'Vertices',vn)
end
xlim(xlimit);
set(gca,'yscale','log')
title(title1)
grid on
grid minor

subplot(2,1,2)
hist(list2,150)
ph2 = get(gca,'children');
% Determine number of histogram patches
N_patches = length(ph2);
for i = 1:N_patches
      % Get patch vertices
      vn = get(ph2(i),'Vertices');
      % Adjust y location
      vn(:,2) = vn(:,2) + 1;
      % Reset data
      set(ph2(i),'Vertices',vn)
end
xlim(xlimit);
set(gca,'yscale','log')
title(title2)
grid on
grid minor
end
end