% %function [ g ] = read_graph( f_path )
clc;
clear;
f_path='../';
f_name='b.out';
% �����ļ�
topo_file = fopen([f_path,f_name], 'r')

for i = 1 : 3
    tline = fgetl(topo_file) ; %fgetl һ�ν�����һ�У���������һ�е��ַ��� 
    if ~ischar(tline), break, end  %����Ƿ����ɹ�����ʧ�����˳�ѭ��
   disp(tline)  %��ʾ������е����ݣ���Ҳ���԰����洢��ĳ���ַ���������
end  %ѭ����������ʱfid ͣ����fgetl������һ��ĩ�ˡ�
%����ǰ7�к���fscanf������
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