function [m,num_seg] = tri3stri( DH0 )
%TRISTRI Summary of this function goes here
%   Detailed explanation goes here
m = DH0;
n = size(DH0,1);
num_seg=zeros(n);
% n=length(m);
for i=1:n
    for j=1:n
        if m(i,j)==0
            m(i,j)=1000000; % 将没有数值的地方用1000000填满
        end
    end
end
m=m-1000000*eye(n); % 去掉对角线，将对角线置零
      
%%
% it depends on the size of the real protein;
inf=1000000;
%%
for i=1:n
    for j=1:n
        if m(j,i)<inf
            %% 关于这里if语句中给num_seg(i,j)赋值1的说明：为什么是1而不是其他大于1
              %的值:这里用反证法：若m(i,j)<inf且 i~=j，num_seg(i,j)==0则num_seg(i,j)=1
              %若num_seg(i,j)>1则推出num_seg(i,j)已被赋值，等价于m(i,j)已经被计算过，则推出num_seg(i,j)~=0，这与等于零矛盾。
              %下边同理。
            if i~=j && num_seg(i,j)==0
                num_seg(i,j)=1; %  将有值得位置除了非对角线元素 置1
            end
            %%
            for k=1:n
                if m(i,k)<inf
                    %%
                    if i~=k && num_seg(i,k)==0
                        num_seg(i,k)=1;
                    end
                    %%
                     s=(m(j,i)+m(i,k));

                    if s<m(j,k)
                        m(j,k)=s;
                        %%
                        num_seg(j,k)=(num_seg(i,j)+num_seg(i,k));%这里得到的数始终对应最小路径的number of segments;1始终对应测定的距离。
  
%                         if num_seg(i,j)>=num_seg(i,k)
%                            num_seg(j,k)=num_seg(i,j)+1+(num_seg(i,k)-1)*rand;
%                         else
%                            num_seg(j,k)=num_seg(i,k)+1+(num_seg(i,j)-1)*rand; 
%                         end
                        %%
                    end
                end
            end
        end
    end
end
% figure;
% mesh(m);

end

