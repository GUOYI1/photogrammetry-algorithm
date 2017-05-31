#include"common.h"
#include"image.h"
#include"matrix.h"
using namespace std;

int main()
{
	Image::parameter();
	Image p1(1),p2(2);
	cout << "Select the target point in image1:" << endl;
	p1.GetXY();
	cout << "Select the corresponding point in image2:" << endl;
	p2.GetXY();
	p1.Resection();														
	p2.Resection();
	Image::ForwardIntersection(p1, p2);
	system("Pause");
	return(0);
}
