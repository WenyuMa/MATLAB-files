% This function is used for finding boundary cycle of each hole

function [boundary_cycle1]=find_cycle1()
 
global boundary_edge cycle_time node1
boundary_cycle1=[];

%% delete isolated edges
boundary_edge0=unique(boundary_edge);
for i=1:length(boundary_edge0)
   c=sum(boundary_edge==boundary_edge0(i));
   c1=0;
   if c>1
       continue;
   else
       i1=find(boundary_edge==boundary_edge0(i));
       j=ceil(i1/2);
       cycle_time(j)=0;
       if j==i1/2
          i1=i1-1; 
       elseif  2*j-1==i1
          i1=i1+1;
       end
       index=find(boundary_edge==boundary_edge(i1));
       if length(index)<3
              for i2=1:length(index)
                  j1=ceil(index(i2)/2);
                  if cycle_time(j1)>0
                     c1=c1+1;
                     break;
                  end  
              end  
              while(c1==1)
                   c1=0;
                   cycle_time(j1)=0;
                   if j1==index(i2)/2
                       i1=index(i2)-1;
                   elseif 2*j1-1==index(i2)
                       i1=index(i2)+1;
                   end   
                   index=find(boundary_edge==boundary_edge(i1));
                   if length(index)<3
                      for i2=1:length(index)
                          j1=ceil(index(i2)/2);
                          if cycle_time(j1)>0
                             c1=c1+1;
                             break; 
                          end     
                      end    
                   end  
              end  
       end 
   end 
end

%% boundary cycle dectect
next=[];
last=[];
ini_edge=[0,0];
over=find(cycle_time>0);
flag1=length(over); 
hole_num=1;

% find initial edge
while flag1
    ini_old=ini_edge;
    k1=0;
    k2=0;
    k3=0;
    a=0;
    x=find(cycle_time==2);
    if x
       for i=1:length(x)
          index1=find(boundary_edge==boundary_edge(2*x(i)-1));
          index2=find(boundary_edge==boundary_edge(2*x(i)));
          if length(index1)==3 
             for i1=1:3
                 if index1(i1)==2*x(i)-1
                    continue;
                 else
                     num1=ceil(index1(i1)/2);
                     if cycle_time(num1)==1
                         a=a+1;
                     end
                 end
             end
             if a==2
                ini_edge=[boundary_edge(2*x(i)-1),boundary_edge(2*x(i))] ;
                ini_weight=cycle_time(x(i));
                k1=1;
                break;
             end
          end
          if length(index2)==3 
              a=0;
              for i1=1:3
                 if index2(i1)==2*x(i)
                    continue;
                 else
                     num1=ceil(index2(i1)/2);
                     if cycle_time(num1)==1
                         a=a+1;
                     end
                 end
              end
              if a==2
                 ini_edge=[boundary_edge(2*x(i)),boundary_edge(2*x(i)-1)] ;
                 ini_weight=cycle_time(x(i));
                 k1=1;
                 break;
              else
                  continue;
              end
          else
              continue;
          end
       end
       if k1==0
              for i=1:length(x)
                  index1=find(boundary_edge==boundary_edge(2*x(i)-1));
                  index2=find(boundary_edge==boundary_edge(2*x(i)));
                  if length(index1)==3 
                      ini_edge=[boundary_edge(2*x(i)-1),boundary_edge(2*x(i))] ;
                      ini_weight=cycle_time(x(i));
                      k2=1;
                      break;
                  elseif length(index2)==3
                      ini_edge=[boundary_edge(2*x(i)),boundary_edge(2*x(i)-1)] ;
                      ini_weight=cycle_time(x(i));
                      k2=1;
                      break;
                  else 
                      continue;
                  end 
              end
              if k2==0
                 number=unidrnd(length(x));
                 y=x(number);
                 ini_edge=[boundary_edge(2*y-1),boundary_edge(2*y)];
                 ini_weight=cycle_time(y);          
              end 
       end
       elseif find(cycle_time==1)  
              x=find(cycle_time==1);
              for i=1:length(x)
                  index1=find(boundary_edge==boundary_edge(2*x(i)-1));
                  index2=find(boundary_edge==boundary_edge(2*x(i)));
                  if length(index1)==2 && length(index2)==2
                     ini_edge=[boundary_edge(2*x(i)-1),boundary_edge(2*x(i))] ;
                     ini_weight=cycle_time(x(i));
                     k3=1;
                     break;
                  else 
                      continue;
                  end 
              end 
              if k3==0
                 number=unidrnd(length(x));
                 y=x(number);
                 ini_edge=[boundary_edge(2*y-1),boundary_edge(2*y)];
                 ini_weight=cycle_time(y);          
              end     
    end 
    if ini_edge==ini_old
        y=find(cycle_time>0);
        number=unidrnd(length(y));
        x=y(number);
        ini_edge=[boundary_edge(2*x-1),boundary_edge(2*x)];
        ini_weight=cycle_time(x);
    end 
    
    %broadcast
    ini=ini_edge(1);
    node1(ini).msg=ini;
    next=ini_edge(2);
    start=ini_edge(2);
    node1(start).msg=[ini,start];
    flag2=length(next);
    while flag2
       last=next;
       next=[];  %
       for i=1:length(last)
           index=find(boundary_edge==last(i));
           flag_find=0;      %独立边标志位
           for j=1:length(index)
                   no_edge=ceil(index(j)/2);
               if cycle_time(no_edge)>0
                   if mod(index(j),2)==1          %奇数偶数判断
                      next_jump=boundary_edge(index(j)+1);
                   else  
                      next_jump=boundary_edge(index(j)-1); 
                   end
                   l=length(node1(last(i)).msg);
                   if l>1 && next_jump==node1(last(i)).msg(l-1)  %回到原来的边，continue
                       flag_re=1;
                       flag_ok=0;
                   else
                       flag_find=1;
                       if ini==next_jump    % return to the initial node?
                          flag_ok=1;
                       else  
                          flag_ok=0;
                       end
                       
                       if node1(next_jump).msg
                          if flag_ok==0 && node1(last(i)).msg(1)== node1(next_jump).msg(1) && node1(last(i)).msg(2)== node1(next_jump).msg(2)  % old message?
                             flag_re=1;
                          elseif flag_ok==1 
                              if l>2 && length(node1(next_jump).msg)>1
                                 if node1(last(i)).msg(1)== node1(next_jump).msg(1) && node1(last(i)).msg(2)== node1(next_jump).msg(2)&&node1(last(i)).msg(3)==node1(next_jump).msg(3)
                                    flag_re=1;
                                 else 
                                    flag_re=0;
                                 end
                              else
                                  flag_re=0;
                              end
                          else  
                             flag_re=0;
                          end 
                       else 
                           flag_re=0;
                       end 
                   end%
                       if flag_ok==0 && flag_re==0   
                          next=[next,next_jump];
                          node1(next_jump).msg=[node1(last(i)).msg,next_jump];
                       elseif flag_ok==0 && flag_re==1
                              continue; 
                       elseif flag_ok==1 && flag_re==1
                              continue;
                       elseif flag_ok==1 && ini_weight>0
                           node1(next_jump).msg=[node1(last(i)).msg,next_jump];   
                           hole=[node1(last(i)).msg,next_jump];
                           boundary_cycle1{hole_num}=hole;   %%%%%
                           hole_num=hole_num+1;
                             % disp(hole);%%%%%%%%
                              for n=1:length(hole)-1
                                  le=find(boundary_edge==hole(n));
                                  for m=1:length(le)
                                      if hole(n)>hole(n+1) 
                                         if hole(n+1)==boundary_edge(le(m)-1)
                                            num=le(m)/2;%
                                            break;
                                         else   
                                            continue;
                                         end  
                                      else 
                                          if hole(n+1)==boundary_edge(le(m)+1)
                                             num=(le(m)+1)/2;%
                                             break;
                                          else  
                                             continue;
                                          end  
                                      end 
                                  end
                                  cycle_time(num)=cycle_time(num)-1;
                              end
                              ini_weight=ini_weight-1;
                       end  
               end 
           end
           if ini_weight==0
               next=[];
               break;
           end
           if flag_find==0               %找不到下一条边
               error=node1(last(i)).msg;
               for n1=1:length(error)-1
                   le=find(boundary_edge==error(n1));
                   for m1=1:length(le)
                       if error(n1)>error(n1+1) 
                          if error(n1+1)==boundary_edge(le(m1)-1)
                             num2=le(m1)/2;%
                             break;
                          else    
                             continue;
                          end   
                       else  
                           if error(n1+1)==boundary_edge(le(m1)+1)
                              num2=(le(m1)+1)/2;%
                              break;
                           else   
                              continue;
                           end   
                       end   
                   end 
                      cycle_time(num2)=cycle_time(num2)-1; 
               end
           end
       end
       flag2=length(next);
    end
    
   % find next initial edge
   over=find(cycle_time>0);
   flag1=length(over);
end