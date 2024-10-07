
close all;

clear all;
input_image = imread('C:\Users\kaush\Downloads\MonuSeg-20200319T073151Z-001\MonuSeg\HeLa-train\19.tif'); 
%input_image=imresize(input_image,[824 1048]);
%input_image=background_detectionb(input_image);
figure
imshow(input_image);
input_image=(input_image);
input_image=imcrop(input_image);
I=input_image;
qz=I;
figure
imshow(input_image);
figure
imshow(input_image);
%qz=rgb2gray(qz);
HSV = rgb2hsv( input_image ) ;
imshow(HSV(:,:,3));
I=HSV(:,:,3);
level = graythresh(I);
    imgbw =  im2bw(I,0.07);
    figure
    imshow(imgbw);
   % imgbw=imcomplement(imgbw);
    BW=imgbw;
bw=BW;
bw7=BW;
bw3=BW;

lab=bwlabel(bw3,8);
hj=bwlabel(bw3,8);
bgt=bw3;

imwrite(bgt,'C:\Users\kaush\Downloads\debopria da\first_phase_potential clumpOverlay_5.tif');
stats = regionprops('table',hj,'Centroid',...
    'MajorAxisLength','MinorAxisLength','Area');
%se=strel('disk',6);
%bw3=imerode(bw3,se);
arer=cat(1,stats.Area);
M = mean(arer,'all');
avg=0;
sum=0;
j=0;
for io=1:numel(arer)
  if(arer(io)>100)
     sum=sum+arer(io);
     j=j+1;
  end
end
avg=sum/j;
 sd3 = cat(1,stats.MajorAxisLength);
 sd4=cat(1,stats.MinorAxisLength);
 arr=  cat(1,stats.Area);
 cent=cat(1,stats.Centroid);
 x_index=cent(:,1);
 y_index=cent(:,2); 
for ui=1:numel(stats)
   if (ui<=numel(arr))
    if((arr(ui)/sd3(ui))>25)
      gh(ui,1)=arr(ui)/sd3(ui);
      gh(ui,2)=x_index(ui);
      gh(ui,3)=y_index(ui);
   else
       gh(ui,1)=0;
       gh(ui,2)=x_index(ui);
      gh(ui,3)=y_index(ui);
   end
   end
end

for ui=1:numel(stats)
   if (ui<=numel(arr))
    if((sd3(ui)-sd4(ui))>5)
      gh(ui,1)=arr(ui)/sd3(ui);
      gh(ui,2)=x_index(ui);
      gh(ui,3)=y_index(ui);
   else
       gh(ui,1)=0;
       gh(ui,2)=x_index(ui);
      gh(ui,3)=y_index(ui);
   end
   end
end

for ui=1:numel(stats)
   if (ui<=numel(arr))
    if(arer(ui)>1000)
      gh(ui,1)=arr(ui)/sd3(ui);
      gh(ui,2)=x_index(ui);
      gh(ui,3)=y_index(ui);
   else
       gh(ui,1)=0;
       gh(ui,2)=x_index(ui);
      gh(ui,3)=y_index(ui);
   end
   end
end



bwe=bw3;


for k = 1:numel(stats)
    if(k<=numel(arr) && gh(k,1)==0)
     j9=ceil(x_index(k));
    j10=ceil(y_index(k));
    j11=ceil(sd3(k));
    
    ccv=j11/2;
    ccv=ceil(ccv);
    
   
    s3=0;
    s1=0;
    s2=0;
    ko=1;
    
    
    for xz=j10-ccv:j10+ccv
        
        for xz1=j9-ccv:j9+ccv
            %bv(xz,xz1)=1;
       if(xz<=size(bw3,1) && xz1<=size(bw3,2) && xz>0 && xz1>0 && hj(j10,j9)==hj(xz,xz1) )
       %disp(xz);
       %disp(xz1);
           bw3(xz,xz1)=0;
      
          
        
      % end
        
        
        end
        
        end
    end
    end
end

figure(9)
imshow(bw3);
bwe1=bwe-bw3;
BW=bw3;
BW=imfill(BW,'holes');
se=strel('disk',2);
 

%332BW=imerode(BW,se);


D = -bwdist(~BW,'quasi-euclidean');
%figure
%imshow(D,[])
%saveas(gcf,'/Users/kaushikiroy/Downloads/cody 2/3c/Overlay_3_seeds.tif');

Ld = watershed(D);

%figure
%imshow(label2rgb(Ld))
%saveas(gcf,'/Users/kaushikiroy/Downloads/cody 2/3c/Overlay_3_seeds.tif');


bw2 = bw;
bw2(Ld == 0) = 0;
%figure
%imshow(bw2)
%saveas(gcf,'/Users/kaushikiroy/Downloads/cody 2/figures_4_b&c/Overlay_1_seed.tif');


mask = imextendedmin(D,2,8);
mask9=mask;
%figure
%imshow(mask);
se=strel('disk',1);
%mask=imdilate(mask,se);
figure(10)
imshow(mask);
imshowpair(BW,mask,'blend')
saveas(gcf,'C:\Users\kaush\Downloads\debopria da\Overlay_3_a.tif');
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
BW(Ld2 == 0) = 0;
%bw3 = bw;
%BW(Ld2 == 0) = 0;
%BW=bw3;
%saveas(gcf,'/Users/kaushikiroy/Downloads/cody 2/figures_4_b&c/Overlay_1_seed.tif');
figure(11)
imshow(BW)
imwrite(BW,'C:\Users\kaush\Downloads\debopria da\Overlay_3_BDT.tif');
%figure (12), imshow(BW);
[b,l]=bwlabel(BW,8);
[B,L,N,A] = bwboundaries(BW);
s = regionprops(BW,'centroid');
centroids = cat(1,s.Centroid);
vr=centroids(:,1);
vr1=centroids(:,2);
figure(13)
imshow(BW);
%BW=imcrop(BW);
BW=bwareaopen(BW,50);
lab=bwlabel(BW,8);
lab1=bwlabel(mask,8);
lastmask=mask;
stats=regionprops('table',lab,'Centroid',...
    'MajorAxisLength','MinorAxisLength', 'Eccentricity','Circularity'); 

stats1=regionprops(lab1,'Centroid',...
    'MajorAxisLength','MinorAxisLength', 'Eccentricity','Circularity'); 
%figure
%imshow(BW);
centroids = cat(1,stats1.Centroid);
vry=centroids(:,1);

vry1=centroids(:,2);
st = regionprops(lab,'Centroid','Orientation','MajorAxisLength');
for k = 1:numel(st)
    c = st(k).Centroid;
    c1 = stats1(k).Centroid;
    hold on;
   plot(c(:,1),c(:,2), 'r.')
   plot(c1(:,1),c1(:,2), 'b.')
hlen = floor(st(k).MajorAxisLength/2);
xCentre = st(k).Centroid(1);
yCentre = st(k).Centroid(2);
cosOrient = cosd(st(k).Orientation);
sinOrient = sind(st(k).Orientation);
xcoords = xCentre + hlen * [cosOrient -cosOrient];
ycoords = yCentre + hlen * [-sinOrient sinOrient];

%line(xcoords, ycoords);
hold off;
end
saveas(gcf,'C:\Users\kaush\Downloads\debopria da\seed_centroid.tif');

for yu=1:numel(st)
   c=st(yu).Centroid;
   c1=stats1(yu).Centroid;
   x_cord=floor(c(:,2));
   y_cord=floor(c(:,1));
   x_cord1=floor(c1(:,2));
   y_cord1=floor(c1(:,1));
  
   if(x_cord<=size(bw,1)&&x_cord1<=size(bw,1) && y_cord<=size(bw,2) && y_cord1<=size(bw,2))
        bn(yu)=lab(x_cord,y_cord);
   bn1(yu)=lab(x_cord1,y_cord1);
   end
end




for yu=1:numel(st)
   c=st(yu).Centroid;
   c1=stats1(yu).Centroid;
   x_cord=floor(c(:,2));
   y_cord=floor(c(:,1));
   x_cord1=floor(c1(:,2));
   y_cord1=floor(c1(:,1));
  
   if(x_cord<=size(bw,1)&&x_cord1<=size(bw,1) && y_cord<=size(bw,2) && y_cord1<=size(bw,2))
        bn(yu)=lab(x_cord,y_cord);
   
   for yu7=1:numel(bn1)
   if(bn(yu)==bn1(yu7))
       c2=stats1(yu7).Centroid;
       x_cord3=floor(c2(:,2));
   y_cord3=floor(c2(:,1));
       dt6=x_cord-x_cord3;
       bio(yu)=dt6;
       dt7=y_cord-y_cord3;
       bio1(yu)=dt7;
       if(dt6<0)
       dt6=dt6*(-1);
       
       end
       if(dt7<0)
           dt7=dt7*(-1);
       end
      
      dt6=dt6*dt6;
      dt7=dt7*dt7;
       totd=dt6+dt7;
       totd=sqrt(totd);
   end
   end
   end
aror(yu)=totd;
center_x(yu)=x_cord;
center_y(yu)=y_cord;
seed_x(yu)=x_cord3;
seed_y(yu)=y_cord3;
end
centroids = cat(1,st.Centroid);
vr=centroids(:,1);

vr1=centroids(:,2);
statsr=regionprops(lab,'MajorAxisLength');
ma=cat(1,statsr.MajorAxisLength);
statsr1=regionprops(lab,'MinorAxisLength');
mi=cat(1,statsr1.MinorAxisLength);
BW3=zeros(size(BW,1),size(BW,2));
BW3=BW;
for k = 1:length(st)
    if(aror(k)>=5)
     j9=ceil(vr(k));
    j10=ceil(vr1(k));
    j11=ceil(ma(k));
    j12=ceil(mi(k));
    ccv=j11/2;
    ccv=ceil(ccv);
    ccv1=j12/2;
    ccv1=ceil(ccv1);
    s3=0;
    s1=0;
    s2=0;
    ko=1;
    
    
    for xz=j10-ccv:j10+ccv
        
        for xz1=j9-ccv:j9+ccv
            if(xz<=size(BW,1) && xz1<=size(BW,2) && xz>0 && xz1>0 && lab(j10,j9)==lab(xz,xz1))
           BW3(xz,xz1)=0;
            end  
        end
    end
    end
end
BW10=BW3;
figure
imshow(BW10);
BW3=BW-BW3;


BW3=bwareaopen(BW3,50);
%figure(14)
%imshow(BW3);
r1=qz(:,:,1);
g1=qz(:,:,2);
b1=qz(:,:,3);
for ui=1:size(qz,1)
    for ui1=1:size(qz,2)
        
        if(BW3(ui,ui1)==0)
            
        r1(ui,ui1)=0;
         g1(ui,ui1)=0;
          b1(ui,ui1)=0;
        end
        
    end
end
qz=cat(3,r1,g1,b1);
figure(14)
imshow(qz);

imwrite(BW3,'C:\Users\kaush\Downloads\debopria da\final phase clump Overlay_3.tif');
labe=bwlabel(BW3,8);
BW4=BW3;
%BW=gh10;
BW=BW3;
[lab67 count]=bwlabel(BW,8);
lr=graythresh(qz);
bwe=im2bw(qz,0.06);
bwe=imfill(bwe,8,'holes');
figure
imshow(bwe);
bw3=bwe;
bw7=imclearborder(bw7);
D = -bwdist(~BW,'quasi-euclidean');
%figure
%imshow(D,[])
%saveas(gcf,'/Users/kaushikiroy/Downloads/cody 2/3c/Overlay_3_seeds.tif');

Ld = watershed(D);

%figure
%imshow(label2rgb(Ld))
%saveas(gcf,'/Users/kaushikiroy/Downloads/cody 2/3c/Overlay_3_seeds.tif');


bw2 = bw;
bw2(Ld == 0) = 0;
%figure
%imshow(bw2)
%saveas(gcf,'/Users/kaushikiroy/Downloads/cody 2/figures_4_b&c/Overlay_1_seed.tif');
figure (18)
imshow(BW4);

mask=imextendedmin(D,1,8);
se=strel('disk',1);
%mask=imdilate(mask,se);

gh9=zeros(size(BW3,1),size(BW3,2));
%[x,y]=getpts;
%x=floor(x);
%y=floor(y);
%gh9(x,y)=1;
%[x,y]=getpts;
%x=floor(x);
%y=floor(y);
%gh9(x,y)=1;

%eta change
%mask = gh10;
se=strel('disk',17);
%mask=imdilate(mask,se);
figure(19)
imshow(mask);
se=strel('disk',4);
%mask=imdilate(mask,se);
figure(20)
imshow(mask);
mask9=mask;
imshowpair(BW,mask,'blend')

D2 = imimposemin(D,mask);
Ld2 = watershed(D2);

%bw3 = bw;
BW(Ld2 == 0) = 0;
%saveas(gcf,'/Users/kaushikiroy/Downloads/cody 2/figures_4_b&c/Overlay_1_seed.tif');
figure(21)
imshow(BW)
BW3=BW;
[lab3,count]=bwlabel(BW,8);
mask=im2bw(mask);

stats11 = regionprops(mask,'centroid');
centroids11 = cat(1,stats11.Centroid);
   
vr1=centroids11(:,1);
vr1=ceil(vr1);

vr2=centroids11(:,2);
vr2=ceil(vr2);
ro=size(mask,1);
ro1=size(mask,2);
mask=im2bw(mask);
for ui=1:size(vr1)
    for ui1=1:size(vr2)
        yu=vr1(ui);
        
        yu1=vr2(ui1);
        if(mask(yu1,yu)==1)
        for tp=yu1-5:yu1+5
            for tp1=yu-5:yu+5
               if(tp>0 && tp<size(mask,1) && tp1>0 && tp1<size(mask,2))
                  mask(tp,tp1)=1; 
               end
            end    
            end
        end
        
        
    end
end
for ui=1:size(mask9,1)
    for ui1=1:size(mask9,2)
        if(qz(ui,ui1)==0)
            mask9(ui,ui1)=0;
            
        end
        
        
    end
end
mask=mask9;

mask=mask9;
mask=im2bw(mask);
%mask=imresize(mask,[256 256]);
%qz=imresize(qz,[256 256]);
stats11 = regionprops(mask,'centroid');
stats121 = regionprops(mask,'PixelIdxList');
centroids11 = cat(1,stats11.Centroid);
   
vr1=centroids11(:,1);
vr1=ceil(vr1);

vr2=centroids11(:,2);
vr2=ceil(vr2);
ro=size(mask,1);
ro1=size(mask,2);
mask=im2bw(mask);
figure
imshow(mask),title("mask");
mask=im2bw(mask);
[la,count]=bwlabel(mask,8);

lk=1;
mask1=zeros(size(mask,1),size(mask,2));
for ty=1:count
    mask1=zeros(size(mask,1),size(mask,2));
    
 for t=1:size(mask,1)
    for t1=1:size(mask,2)
        
        if(la(t,t1)==lk)
           mask1(t,t1)=1;
           m9{lk,1}=mask1; 
           % lk=lk+1;
        tr(lk,1)=lab67(t,t1);
        end
        
    end
    
 end
lk=lk+1; 
end



figure
imshow(I);

rt=find(mask(:,1)==1);
rt1=find(mask(:,2)==1);
seg1=ones(size(mask,1),size(mask,2));
seg3=zeros(size(mask,1),size(mask,2));
s=0;
flag=0;
po=1;
%[seg] = region_seg_1(qz,mask,200);

figure
imshow(qz)

for j = 1:length(m9)
    
    flag9=0;
      seg=ones(size(seg3,1),size(seg3,2));        
   seg1=ones(size(seg3,1),size(seg3,2));
          phi1=m9{j,1};
          
          yu=vr1(j);
        
          yu1=vr2(j);
          %po=po+1;
          
          for j1 = 1:length(m9)
              flag9=0;
          if(j1~=j) 
            yu2=vr1(j1);
        
          yu3=vr2(j1);  
          d1=yu2-yu;
          d2=yu3-yu1;
          if (d1<0)
              d1=d1*(-1);
          end
          
          if(d2<0)
              
              
              d2=d2*(-1);
              
          end
          s=d1+d2;
          if (s<200 && tr(j)==tr(j1) )
          phi_2= m9{j1,1};
         for yu=1:size(ui,1)
             if ((j==ui(yu,1)&&j1==ui(yu,2))||(j1==ui(yu,1)&&j==ui(yu,2)))
             flag9=1;
             pos5=yu;
            
             end
         end
         if(flag9==0)
         ui(po,1)=j;
          ui(po,2)=j1;
          %[seg] = region_seg(qz,phi1,200);
          [seg,seg2] = region_seg(qz,phi1,phi_2,200);
          seg90{po,1}=seg2;
         po=po+1;
         end
         if (flag9==1)
            seg=seg90{pos5,1}; 
            disp(pos5);
         figure
         imshow(seg);
         end
          
          %figure
          seg1=seg1&seg ;
          %imshow(seg);
          
          flag=1;
         
          end
         
          
          end
          
           
          s=0;
          end
          if (flag==0)
          
          [seg] = region_seg_1(qz,phi1,200);
          seg1=seg1&seg;

           %[seg] = region_seg(qz,phi1,phi_2,200); %-- Run segmentation


          
          end
        flag=0;     
figure
imshow(seg1); title('Global Region-Based Segmentation');
seg=imfill(seg1,8,'holes');
%imwrite(seg,'/Users/Chanlab/Downloads/codes&images/data_set_205_images/images/103_ty.tif');
%seg = bwconvhull(seg,'objects');
seg3=seg3+seg1;
%figure
%imshow(seg1);          
s=0;  
       
end
figure
imshow(seg3),title('final');
[lab,count]=bwlabel(seg3,8);
seg9=seg3+bwe1+BW10;
figure
imshow(seg9);



%figure
%imshow(BW11);
%imwrite(BW11,'C:\Users\kaush\Downloads\all_data\T-47D_out\Overlay_63.tif');