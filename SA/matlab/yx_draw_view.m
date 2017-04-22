% %function [ g ] = read_graph( f_path )
clc;
clear;
f_path='../';
f_name='b.out';
% 读入文件
topo_file = fopen([f_path,f_name], 'r')

for i = 1 : 3
    tline = fgetl(topo_file) ; %fgetl 一次仅读入一行，并返回这一行的字符串 
    if ~ischar(tline), break, end  %检查是否读入成功，若失败则退出循环
   disp(tline)  %显示读入的行的内容，你也可以把它存储到某个字符串变量中
end  %循环结束，此时fid 停留在fgetl最后读的一行末端。
%读完前7行后，用fscanf继续读
% Gen = 2, Average Cost = 168959.670000, Min Cost = 157376, Current MinCost = 157376, Current worstCost = 209434
%a= 'Gen = 67, PopScale = 50, Average Cost = 1995.680000, Average Fit = 335.032384, Min Cost = 1967, Current MinCost = 1967, Current worstCost = 2207';
A = fscanf(topo_file, '%*s %*c %d%*c %*s %*s %*c %f%*c %*s %*s %*c %d%*c %*s %*s %*c %d%*c %*s %*s %*c %d',[5 300]);
x=A(1,:);
y=A(2,:);
minCost=min(A(3,:));

plot(x,y);
hold on
plot(x,A(4,:),'r-');
hold on
plot(x,A(5,:),'g-');
legend('Average','Current Best','Current Worst');
xlabel('Generation')
ylabel('Cost')
title(minCost);