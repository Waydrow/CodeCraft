% %function [ g ] = read_graph( f_path )

% f_path='/';
% f_name='a.txt';
% % �����ļ�
% topo_file = fopen([f_path,f_name], 'r')
% 
% for i = 1 : 4
%     tline = fgetl(topo_file) ; %fgetl һ�ν�����һ�У���������һ�е��ַ��� 
%     if ~ischar(tline), break, end  %����Ƿ����ɹ�����ʧ�����˳�ѭ��
%    disp(tline)  %��ʾ������е����ݣ���Ҳ���԰����洢��ĳ���ַ���������
% end  %ѭ����������ʱfid ͣ����fgetl������һ��ĩ�ˡ�
% %����ǰ7�к���fscanf������
% A = fscanf(topo_file, '%*s %*c %d%*c  %*s %*s %*c %f%*c %*s %*s %*c %*f%*c %*s %*s %*c %d',[3 400]);
% x=A(1,:);
% y=A(2,:);
% minCost=min(A(3,:));
% subplot(1,2, 1);
% plot(x,y);

clc;
clear;
f_path='/';
f_name='b.out';
file=fopen([f_path f_name],'r');
for i=1:4
    tline=fgetl(file);
    if ~ischar(tline),break;end
    disp(tline);
end
A=[];
B=[];
for i=1:24726-4
    if mod(i,123)==1
     j=1; 
    temp=fgetl(file);
    A=[A;sscanf(temp,'%*s %*c %d%*c  %*s %*s %*c %f%*c %*s %*s %*c %*f%*c %*s %*s %*c %d',[1,3])];
 
    else
      %temp=fgetl(file);
        if(mod(j,2)==0&&j<=120)
             temp=fgetl(file)
             
             B=[B;sscanf(temp,'%*s %d %*s %d  %*s %*s %*d',[1,2])];
        else
            temp=fgetl(file);
        end
         j=j+1;
   end
     
       if ~ischar(tline),break;end
    
    
end


%   
  num=[];
   for i=1:60:12060
       
      
        C=[];
        for j=i:i+59
        C=[C;B(j,:)];
       
        end
         table = tabulate(C(:,1));
         [maxCount,idx] = max(table(:,2));
         num=[num ;table(idx)];

       
       
   end
   plot(1:201,num');











