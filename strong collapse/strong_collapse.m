%% deploy nodes
close all;
clc;
clear;

% lambda is the intensity 
lambda=10;

nmb=poissrnd(lambda);

x = rand(1,nmb);
y = rand(1,nmb);
x = 4*x + 1;
y = 4*y + 1;

%coordinates of fence nodes
%fence_x1 = [0,0,0,0,0, 0, 0,2, 2,4, 4,6, 6,8, 8,10,10,12,12,12,12,12,12,12];
%fence_y1 = [0,2,4,6,8,10,12,0,12,0,12,0,12,0,12, 0,12, 0, 2, 4, 6, 8,10,12];
fence_x1=[0,0,0,0,2,2,4,4,6,6,6,6];
fence_y1=[0,2,4,6,0,6,0,6,0,2,4,6];

%inner_x1 = [1,1,1,1,1, 1,3, 3,5, 5,7, 7,9, 9,11,11,11,11,11,11];
%inner_y1 = [1,3,5,7,9,11,1,11,1,11,1,11,1,11, 1, 3, 5, 7, 9,11];
inner_x1=[1,1,1,3,3,5,5,5];
inner_y1=[1,3,5,1,5,1,3,5];

%coordinates of all the nodes
node_x = [x, inner_x1, fence_x1];
node_y = [y, inner_y1, fence_y1]; 

node_coor = [node_x; node_y];

plot(node_x, node_y, '.'); hold on;

% figure
Theta=[0:0.005:1]*2*pi;
Xcircle=1.0*cos(Theta);   % coordinates of the 200 points of a typical coverage disc
Ycircle=1.0*sin(Theta);
for i=1:length(node_coor)        % ����fence node֮�����нڵ��Rs  /24
    Xc=Xcircle+node_x(i);
    Yc=Ycircle+node_y(i);
    plot(Xc,Yc,'k');
    fill(Xc,Yc,'g','facealpha',0.5);  %��͸����ʾ
    
    axis square; %��������������ϵ
    xlim([0 6]); %x,y������������
    ylim([0 6]);
    hold on;
end

hold on;

%��ӡ�ڵ����
for i = 1: length(node_coor) 
    text(node_x(i)+0.1, node_y(i), num2str(i));
end

for i=1: length(node_coor) 
    node(i).neighbor = [];
    node(i).weight=0; 
    node(i).simp=[];
    node(i).status=1;
    
    if i > length(node_coor) - 12
        node(i).fence_flag = 1;
    else
        node(i).fence_flag = 0;    
    end
    
    for j=1: length(node_coor)
        if (j==i) continue;
        else
            distance = dist2(node_coor(:,i),node_coor(:,j));
            if distance <= 2 
                node(i).neighbor = [node(i).neighbor, j];%�ҵ� node��i)�������ھӽڵ�  ��������1-simplex
            end
        end
    end
end

% ������������������
for i = 1 : length(node_coor)
    new_simp = struct('vert', [i], 'neighb', node(i).neighbor);
    node(i).simp{1} = new_simp;
end
%Ȼ�����εݹ� �ҳ�ÿ�������k-simplex
for i = 1 : length(node_coor)
    while(1)
        index_max = size(node(i).simp, 2);%���ĵ���,���ؾ��������

        no_max_simp = size(node(i).simp{index_max}, 2); %????
        %����ε������Ŀ
%����ÿ����������ά������no_max_simp��   ����ÿ�����ά����
% ����
        for j = 1 : no_max_simp
            vert_set = node(i).simp{index_max}(j).vert;%��j�����ά���εĶ���ļ��� ����j-1���εĶ��㼯��    (��0-simplex��ʼ���㣩
            neighb_set = node(i).simp{index_max}(j).neighb;%��j�����ά���ε��ھӼ���

            if ~isempty(neighb_set)
                
                no_neighb = length(neighb_set);
                for k = 1 : no_neighb     %�����j�����ά���� �ھӼ��ϵ�ÿ������
                    new_node = neighb_set(k);%���ڹ��ɸõ�j�����ε��ھӼ��ϵĵ�k������
                    if new_node < max(setdiff(vert_set, i))  %setdiff(vert_set, i)����vert_set�д��ڶ�i�в����ڵ�Ԫ��
                        continue;
                    else
                       
                        vert_set_temp = union(vert_set, new_node);
                        neighb_set_temp = intersect(neighb_set, node(new_node).neighbor);%��ǰi������ھӼ���  ��  ��ǰi���� ��k���ھӶ�����ھӼ���  �Ľ���   
                        simp_temp = struct('vert', vert_set_temp, 'neighb', neighb_set_temp);%�ɵ�i������͵�k�����㹹��  �Լ����ǹ�ͬ���ھӶ��㹹�ɵ�higher simplex

                        if index_max == size(node(i).simp, 2) %this is the first index_max+1 new higher simplex �����index_max����ά�ĸ��������ڸ�ά��һ��
                            node(i).simp{index_max+1}(1) = simp_temp; %��һ������ά�ĵ��Σ�ֱ�����
                        else                                    %this is the next index_max+1 higher simplex   ��ά�ĸ��Σ���
                            no_new_simp = size(node(i).simp{index_max+1}, 2);%�����ǵ�һ��������Ҫ�����±�Ӧ���Ƕ���
                            node(i).simp{index_max+1}(no_new_simp+1) = simp_temp;%index_max+1  �����index_max����ά�ĸ��������ڸ�ά��no_new_simp+1��
                        end
                    end
                end
            end
        end

        if index_max == size(node(i).simp, 2) % there is no higher simplex
            break;
        end
    end
   
end

% %weightȨ�ؼ���
% for i=1:nmb
%     node(i).weight=0;
%     index_max=size(node(i).simp,2);
%     if index_max>=2&&~isempty(node(i).simp{2}) 
%         edges_index=size(node(i).simp{2},2);
%         for j=1:edges_index
%             if isempty(node(i).simp{2}(j).neighb)
%                 node(i).weight=0;
%                 
%                 break;
%             else if index_max>=3&&~isempty(node(i).simp{3})
%                     tri_index=size(node(i).simp{3},2);
%                     for t=1:tri_index
%                         if isempty(node(i).simp{3}(t).neighb)
%                             node(i).weight=2;
%                             break;
%                         else
%                             node(i).weight=3;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% %����weight_max
% weight_max=0;
% for i=1:nmb
%     weight(i)=node(i).weight;
% end