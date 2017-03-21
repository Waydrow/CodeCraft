% %function [ g ] = read_graph( f_path )

% f_path='/';
% f_name='a.txt';
% % 读入文件
% topo_file = fopen([f_path,f_name], 'r')
% 
% for i = 1 : 4
%     tline = fgetl(topo_file) ; %fgetl 一次仅读入一行，并返回这一行的字符串 
%     if ~ischar(tline), break, end  %检查是否读入成功，若失败则退出循环
%    disp(tline)  %显示读入的行的内容，你也可以把它存储到某个字符串变量中
% end  %循环结束，此时fid 停留在fgetl最后读的一行末端。
% %读完前7行后，用fscanf继续读
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











