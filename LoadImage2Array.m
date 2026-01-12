 
%read image and convert to 3D-array
function [ImageArray, ImageArrayFULL] = LoadImage2Array(imagename,myXYZ,Xshift,Yshift,PixelSize)
    myImage=double(imread(imagename));

    marg.top=max(myXYZ(:,1));
    marg.bottom=min(myXYZ(:,1));
    marg.right=max(myXYZ(:,2));
    marg.left=min(myXYZ(:,2));

    maxIM = max(myImage, [],'all');
    Threshold = maxIM-maxIM/5;

    [ImgY,ImgX]=meshgrid(0:(size(myImage,1)-1), 0:(size(myImage,2)-1));
    intensity = [];
    for c= 1:1:size(myImage,1)
         row=myImage(c,:)';
         % row(row(:) > Threshold) = [Threshold];  
         intensity=cat(1,row, intensity);
    end

    ImgY=-1*(ImgY-size(myImage,1));
    ImgX=(ImgX+0)*PixelSize;
    ImgY=(ImgY-0.5)*PixelSize;
    ImgX=ImgX+Xshift;
    ImgY=ImgY+Yshift;


    % 
    % ImgX=ImgX-0.5;
    % ImgY=ImgY+0.5;
    % ImgX=(ImgX*PixelSize)+Xshift;
    % ImgYOffset=(size(myImage,1)*PixelSize)+Yshift;
    % ImgY=-1*((ImgY*PixelSize)-ImgYOffset);
    % ImgX=ImgX-PixelSize/2;
    % ImgY=ImgY-PixelSize/2;



    ImageArray=[ImgX(:),ImgY(:),intensity];
    ImageArrayFULL = ImageArray;
    ImageArray(ImageArray(:,2)<marg.left,:)=[];
    ImageArray(ImageArray(:,2)>marg.right,:)=[];
    ImageArray(ImageArray(:,1)>marg.top,:)=[];
    ImageArray(ImageArray(:,1)<marg.bottom,:)=[];

end


