Photogrammetry Algorithm
======================
A collection of the implementation of photogrammetry algorithm,including Resection, Calculation of elements of exterior orientation and Forward Intersection.
======================
Before compiling the source code, please make sure that you have configured the OpenCV open source library properly.       
To compile the source code, first compile the *readpixel.cpp* and produce the  *readpixel.exe*. Then, compile the *resection.cpp* and *forwardintersection.cpp*. Both these 2 solutions will call the *readpixel.exe* we produced previously.     

### Resection
Given the reference image as below, the solution of *Resection* is used to compute the length, width and height of the corridor.In the image below, the O,A,B,C are the 4 control points we set for the image.
![](./photo/reference_image1.png)    
     
### Forward Intersction
Given 2 reference images as below, the solution of Forward Intersction is used to compute the elements of exterior orientation, projection parameters and the coordinates of the image points in the real world. In the images below, A,B,C,D are the 4 control points we selected from the IMG_LEFT and A',B',C',D' are the corresponding points in IMG_RIGHT.       
In the source code, we hard code the default size for an image as 5616*3744. Also,we hard code the pixel coordinates of control points shown in the images below:       
IMG_LEFT: A(2783£¬1758),B(3116£¬1760),C(2778£¬2092),D(3116£¬2096)          
IMG_RIGHT(Corresponding Points): A'(2477£¬1768),B'(2809£¬1763),C'(2474£¬2102),D'(2811£¬2100)     
##### IMG_LEFT
![](./photo/reference_image2.png) 
##### IMG_RIGHT
![](./photo/reference_image3.png) 
