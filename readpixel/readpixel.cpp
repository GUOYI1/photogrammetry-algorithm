#include<cv.h>
#include<highgui.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include<fstream>
using namespace std;
struct Readpixel
{
	int x;
	int y;
};

int n = 0;
IplImage *img = 0;
Readpixel rp;
void on_mouse(int event, int x, int y, int flags, void *usfc)
{
	CvFont font;
	cvInitFont(&font, CV_FONT_HERSHEY_TRIPLEX, 4.0, 4.0, 0, 8, CV_AA);
	if (event == CV_EVENT_LBUTTONDOWN)
	{
		CvPoint pt = cvPoint(x, y);
		rp.x = pt.x;
		rp.y = pt.y;
		cout << "(" << pt.x << "," << pt.y << ")" << endl;
		cout << "OK" << endl;
		char temp[16];
		sprintf(temp, "  (%d,%d)", pt.x, pt.y);
		cvPutText(img, temp, pt, &font, cvScalar(0, 0, 0, 0));
		cvCircle(img, pt, 10, cvScalar(255, 0, 0), CV_FILLED, CV_AA, 0);
		cvShowImage("Image", img);
		cout << "Press Enter to reture the value of coordinate..." << endl;
	}
}

int main()
{
	char image_path[100];
	cout << "Please input the path of the image or Drag the picture into the frame." << endl;
	cin >> image_path;
	cout << "Tips£ºAfter loading an image,please notice the information on the console" << endl;
	system("Pause");
	cout << "================================" << endl;
	cout << "Select the target in the image" << endl;
	img = cvLoadImage(image_path);
	cvNamedWindow("Image", 0);
	cvResizeWindow("Image", 900, 600);
	cvSetMouseCallback("Image", on_mouse, 0);
	cvShowImage("Image", img);
	cvWaitKey(0);
	cvDestroyWindow("Image");
	cvReleaseImage(&img);
	ofstream outfile;
	outfile.open("temp.dat", ios::out);
	outfile << rp.x << " " << rp.y << endl;
	outfile.close();
	return 0;
}