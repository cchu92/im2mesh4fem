%% read png and find countours

% img=imread('ex_blob.png');

%% make pretend bitmap image
img_size=512;

img=zeros(img_size,img_size);

% make two squares of different values
sq1_c1=[100,100];
sq1_c2=[400,400];

sq2_c1=[150,200];
sq2_c2=[250,300];

val1=255;
val2=150;

img(sq1_c1(1):sq1_c2(1),sq1_c1(2):sq1_c2(2))= val1;
img(sq2_c1(1):sq2_c2(1),sq2_c1(2):sq2_c2(2))= val2;

%% find contours and get edges

%assume image has been read in
img1=double(img == val1 | img == val2); % background layer
img2=double(img == val2); % inside layer

figure
subplot(1,3,1);
imshow(img,[]);
subplot(1,3,2);
imshow(img1,[]);
subplot(1,3,3);
imshow(img2,[]);


C=bwboundaries(img1);
C=C{1};
C1=DecimatePoly(C,[1e-4 1]);

figure
subplot(1,2,1);
hold on
imshow(img)
plot(C(:,2),C(:,1),'.-','MarkerSize',10), axis equal
axis('tight')
axis([0 img_size 0 img_size])

subplot(1,2,2)
hold on
imshow(img)
plot(C1(:,2),C1(:,1),'.-','MarkerSize',10), axis equal
axis('tight')
axis([0 img_size 0 img_size])


C=bwboundaries(img2);
C=C{1};
C2=DecimatePoly(C,[1e-4 1]); %from mathworks exchange

figure
subplot(1,2,1);

hold on
imshow(img2)
plot(C(:,2),C(:,1),'.-','MarkerSize',10), axis equal
axis('tight')
axis([0 img_size 0 img_size])

subplot(1,2,2)
hold on
imshow(img2)
plot(C2(:,2),C2(:,1),'.-','MarkerSize',10), axis equal
axis('tight')
axis([0 img_size 0 img_size])



%should be unique nodes and closed - i.e. first and last the same, so we
%dont need last one
node1=C1(1:end-1,:);


nn=size(node1,1);
edge1= [(1:nn-1)',(2:nn-0)'; nn,1] ;

node2=C2(1:end-1,:);
nn=size(node2,1);
edge2= [(1:nn-1)',(2:nn-0)'; nn,1] ;


figure;

hold on; axis image on;
patch('faces',edge1(:,1:2),'vertices',node1, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;
patch('faces',edge2(:,1:2),'vertices',node2, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;
axis([0 img_size 0 img_size])

%% make mesh with simple sizing function

% edge2=edge2+size(node1,1);
node=[node1; node2];
edge=[edge1;edge2+size(node1,1)];
partid=[ones(size(edge1,1),1);2*ones(size(edge2,1),1)]';

part{1} = [find(partid == 1) find(partid == 2) ] ;
part{2} = [find(partid == 2)] ;

%sizing function
hmax = +5 ;
[vlfs,tlfs,hlfs] = lfshfn2(node,edge,part) ;
% apply minimum
hlfs = min(hmax,hlfs) ;

[slfs] = idxtri2(vlfs,tlfs) ;
hfun = @trihfn2;
[vert,etri,tria,tnum] = refine2(node,edge,part,[],hfun,vlfs,tlfs,slfs,hlfs) ;
% [vert,etri,tria,tnum] = refine2(node,edge,part);

[vert,etri,tria,tnum] = smooth2(vert,etri,tria,tnum) ;

figure;
patch('faces',tria(tnum==1,1:3),'vertices',vert, ...
    'facecolor',[1.,1.,1.], ...
    'edgecolor',[.2,.2,.2]) ;
hold on; axis image on;
patch('faces',tria(tnum==2,1:3),'vertices',vert, ...
    'facecolor',[.9,.9,.9], ...
    'edgecolor',[.2,.2,.2]) ;
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;
title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(tria,1))]) ;

tricost(vert,etri,tria,tnum);
save('mesh2dex','vert','etri','tria','tnum');
%% Make mesh with custom sizing field
hfun = @hfuncirc ;

%make simple mesh to evaluate the sizing field
[vert,etri,tria,tnum] = refine2(node,edge,part,[],10);
[vert,etri,tria,tnum] = smooth2(vert,etri,tria,tnum) ;
% plot the sizing field results
hvrt = feval(hfun,vert);

figure;
patch('faces',tria(:,1:3),'vertices',vert , ...
    'facevertexcdata' , hvrt, ...
    'facecolor','interp', ...
    'edgecolor','none') ;
% hold on; axis image off;
hold on; axis image on;
patch('faces',tria(tnum==1,1:3),'vertices',vert, ...
    'facecolor','none', ...
    'edgecolor',[.2,.2,.2]) ;
patch('faces',tria(tnum==2,1:3),'vertices',vert, ...
    'facecolor','none', ...
    'edgecolor',[.2,.2,.2]) ;
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor','none', ...
    'linewidth',1.5) ;


%make new mesh using the real sizing field 
[vert,etri,tria,tnum] = refine2(node,edge,part,[],hfun) ;
% [vert,etri,tria,tnum] = refine2(node,edge,part);

[vert,etri,tria,tnum] = smooth2(vert,etri,tria,tnum) ;



figure;
patch('faces',tria(tnum==1,1:3),'vertices',vert, ...
    'facecolor',[1.,1.,1.], ...
    'edgecolor',[.2,.2,.2]) ;
hold on; axis image on;
patch('faces',tria(tnum==2,1:3),'vertices',vert, ...
    'facecolor',[.9,.9,.9], ...
    'edgecolor',[.2,.2,.2]) ;
patch('faces',edge(:,1:2),'vertices',node, ...
    'facecolor','w', ...
    'edgecolor',[.1,.1,.1], ...
    'linewidth',1.5) ;
title(['MESH-OPT.: KIND=DELFRONT, |TRIA|=', ...
    num2str(size(tria,1))]) ;

tricost(vert,etri,tria,tnum);

% save('mesh2dex','vert','etri','tria','tnum');

function [hfun] = hfuncirc(test)
%HFUN8 user-defined mesh-size function for DEMO-8.

hmax = 25 ;
hmin = 5 ;

xmid = 400 ;
ymid = 250 ;
radius=50;

hfun=hmax*ones(size(test,1),1);

% hcir = exp( -.5*(test(:,1)-xmid).^2 -2.*(test(:,2)-ymid).^2 ) ;
% hcir=0;

dist=test-[xmid ymid];
dist=sum(dist.^2,2).^0.5;

hfun(dist < radius) = hmin;


end