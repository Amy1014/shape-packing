//#include "StdAfx.h"
#include "CmAffine.h"
#include "highgui.h"
#include <map>
#include <math.h>
using std::pair;
using namespace std;
#define PI 3.14159265358979323846

double sqrDist(cv::Point2d& p, cv::Point2d& q){
	cv::Vec2d v= p-q;
	return v.dot(v);
}
double round(double number)
{
	return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

bool CmAffine::ExtractBoundary(const CvMat *img8U, PointSet &pntSet)
{
	CvMat *cannyImg = cvCreateMat(img8U->rows, img8U->cols, CV_8UC1);
	cvCanny(img8U, cannyImg, 200, 400, 3); 
	pntSet.clear();

	for (int r = 0; r < cannyImg->rows; r++)
	{
		byte* data = cannyImg->data.ptr + cannyImg->step * r;   
		for (int c = 0; c < cannyImg->cols; c++)
		{
			if (data[c])
				pntSet.push_back(Point(c, r));
		}
	}
	cvReleaseMat(&cannyImg);

	return pntSet.size() > 0;
}

bool CmAffine::GetSampling(PointSet &pntSet, int samplenum)
{
	if (pntSet.size() < (size_t)samplenum)
		return false;
 
	int	startSize = min(int(pntSet.size()), samplenum * 3);

	vector<int>	idlist(pntSet.size());
	vector<bool> alive(pntSet.size(), true);

	for (unsigned int i = 0; i < idlist.size(); i ++)
		idlist[i] = i;

	random_shuffle(idlist.begin(), idlist.end());

	vector<pair<double, pair<int,int>>>	pwDist(startSize * (startSize - 1) / 2);
	int	cc = 0;

	for (int i = 0; i < startSize; i ++)
		for (int j = i + 1; j < startSize; j ++)
			pwDist[cc ++] = make_pair(sqrDist(pntSet[idlist[i]], pntSet[idlist[j]]), make_pair(idlist[i], idlist[j]));
	sort(pwDist.begin(), pwDist.end());

	for (int nodeLeft = startSize, i = 0; nodeLeft > samplenum && i < cc; i ++)
	{
		if (alive[pwDist[i].second.first] && alive[pwDist[i].second.second])
			alive[pwDist[i].second.first] = false, nodeLeft --;
	}

	PointSet	pntSetNew(samplenum);

	cc = 0;
	for (int i = 0; i < startSize; i ++)
	{
		if (alive[idlist[i]])
			pntSetNew[cc++] = pntSet[idlist[i]];
	}
	
	pntSet = pntSetNew;
	return true;
}

/************************************************************************/
/* srcP64F, dst64F: 2 by N matrix composed by (x, y) of N points		*/
/* S64F:    a 2 by 2 matrix, dst = S64F * srcP64F;						*/
/************************************************************************/
void CmAffine::Reduction(IN CvMat *srcP64F, OUT CvMat *S64F, OUT CvMat *dst64F)
{
	CvMat *covMat = cvCreateMat(2, 2, CV_64FC1);
	CvMat *EV = cvCreateMat(2,2,CV_64F);
	CvMat *ev = cvCreateMat(2,1,CV_64F);
	CvMat *temp = cvCreateMat(2,2,CV_64F);
	CvMat *mean = cvCreateMat(2, 1, CV_64F);

	cvCalcCovarMatrix((const CvArr **)&srcP64F, srcP64F->cols, covMat, mean, CV_COVAR_NORMAL | CV_COVAR_COLS);
	cvEigenVV(covMat, EV, ev);
	cvZero(S64F);
	CV_MAT_ELEM(*S64F, double, 0, 0) = 1 / sqrt(ev->data.db[0]);
	CV_MAT_ELEM(*S64F, double, 1, 1) = 1 / sqrt(ev->data.db[1]);

	cvGEMM(S64F, EV, 1, NULL, 0, temp, 0);
	cvGEMM(EV, temp, 1, NULL, 0, S64F, CV_GEMM_A_T);
	cvGEMM(S64F, srcP64F, 1, NULL, 0, dst64F, 0);

	cvReleaseMat(&mean);
	cvReleaseMat(&covMat);
	cvReleaseMat(&EV);
	cvReleaseMat(&ev);
	cvReleaseMat(&temp);
}


/************************************************************************/
/* Pts1 is 2 by N1 matrix												*/
/* Pts2 is 2 by N2 matrix												*/
/* AA is 2x2 matrix and shift is 2x1 matrix								*/
/* AA is the affine transform of the two point set that AA*Pts1+shift 	*/
/*       best approximate Pts2.											*/
/************************************************************************/
double CmAffine::GetAffine2D(IN CvMat *Pts1, IN CvMat *Pts2, OUT CvMat *AA, OUT CvMat *shift)
{
	int N1 = Pts1->width, N2 = Pts2->width;
	assert(Pts1->height == 2 && Pts2->height == 2);

	CvMat *PP1 = cvCreateMat(2, N1, CV_64FC1);
	CvMat *PP2 = cvCreateMat(2, N2, CV_64FC1);
	double xMean1, xMean2, yMean1, yMean2;
	// Centerize input array and produce output point sets: PP1 and PP2
	{
		CvMat subMatS, subMatD;
		CvScalar temp;
		cvGetRow(Pts1, &subMatS, 0);
		cvGetRow(PP1, &subMatD, 0);
		temp = cvAvg(&subMatS);
		cvSubS(&subMatS, temp, &subMatD);
		xMean1 = temp.val[0];

		cvGetRow(Pts1, &subMatS, 1);
		cvGetRow(PP1, &subMatD, 1);
		temp = cvAvg(&subMatS);
		cvSubS(&subMatS, temp, &subMatD);
		yMean1 = temp.val[0];

		cvGetRow(Pts2, &subMatS, 0);
		cvGetRow(PP2, &subMatD, 0);
		temp = cvAvg(&subMatS);
		cvSubS(&subMatS, temp, &subMatD);
		xMean2 = temp.val[0];

		cvGetRow(Pts2, &subMatS, 1);
		cvGetRow(PP2, &subMatD, 1);
		temp = cvAvg(&subMatS);
		cvSubS(&subMatS, temp, &subMatD);
		yMean2 = temp.val[0];
	}

	// Balance the shape in different directions
	CvMat *S1 = cvCreateMat(2, 2, CV_64FC1);
	CvMat *S2 = cvCreateMat(2, 2, CV_64FC1);
	CvMat *P1 = cvCreateMat(2, N1, CV_64FC1);
	CvMat *P2 = cvCreateMat(2, N2, CV_64FC1);
	Reduction(PP1, S1, P1);
	Reduction(PP2, S2, P2);

	//CvMat* showImg = cvCreateMat(800, 800, CV_8UC3);
	//cvZero(showImg);
	//DisplayPointSet(P1, cvScalar(255, 0, 0, 0), NULL, showImg);
	//DisplayPointSet(P2, cvScalar(0, 0, 255, 0), CmSprintf("%sP1P2R.jpg", outDir.c_str()), showImg);
	//cvReleaseMat(&showImg);

	// Third moment
	CvMat X1, Y1, X2, Y2;
	CvMat *M1_3 = cvCreateMat(4, 1, CV_64FC1);
	CvMat *M2_3 = cvCreateMat(4, 1, CV_64FC1);
	{
		cvScale(P1, P1, 1.0/N1);
		cvScale(P2, P2, 1.0/N2);

		cvGetRow(P1, &X1, 0);
		cvGetRow(P1, &Y1, 1);
		cvGetRow(P2, &X2, 0);
		cvGetRow(P2, &Y2, 1);

		// Find third moment for P1: M1_3
		CvMat *temp = cvCreateMat(1, N1, CV_64FC1);
		cvPow(&X1, temp, 3);
		M1_3->data.db[0] = cvSum(temp).val[0];
		cvPow(&X1, temp, 2);
		M1_3->data.db[1] = cvDotProduct(temp, &Y1);
		cvPow(&Y1, temp, 2);
		M1_3->data.db[2] = cvDotProduct(temp, &X1);
		cvPow(&Y1, temp, 3);
		M1_3->data.db[3] = cvSum(temp).val[0];
		cvReleaseMat(&temp);
		
		// Find third moment for P2: M2_3
		temp = cvCreateMat(1, N2, CV_64FC1);
		cvPow(&X2, temp, 3);
		M2_3->data.db[0] = cvSum(temp).val[0];
		cvPow(&X2, temp, 2);
		M2_3->data.db[1] = cvDotProduct(temp, &Y2);
		cvPow(&Y2, temp, 2);
		M2_3->data.db[2] = cvDotProduct(temp, &X2);
		cvPow(&Y2, temp, 3);
		M2_3->data.db[3] = cvSum(temp).val[0];
		cvReleaseMat(&temp);
	}

	// Get rotation
	double minCostTheta;
	{
		double minCost = numeric_limits<double>::max();
		double MM[4];
		double *d = M1_3->data.db;
		for (double r = -PI, rEnd = PI, rStep = PI / 180; r < rEnd; r += rStep)
		{
			double A = cos(r), B = -sin(r); 
			double C = -B, D = A;

			MM[0] = A*A*A*d[0] + 3*A*A*B*d[1] + 3*A*B*B*d[2] + B*B*B*d[3];
			MM[3] = C*C*C*d[0] + 3*C*C*D*d[1] + 3*C*D*D*d[2] + D*D*D*d[3];

			MM[1] = A*A*C*d[0] + (2*A*B*C+A*A*D)*d[1] + (2*A*B*D+B*B*C)*d[2] + B*B*D*d[3];
			MM[2] = C*C*A*d[0] + (2*C*D*A+C*C*B)*d[1] + (2*C*D*B+D*D*A)*d[2] + D*D*B*d[3];

			double GX = 0;
			for (int i = 0; i < 4; i++)
			{
				double KX = MM[i] - M2_3->data.db[i];
				GX += KX * KX;
			}
			if (GX < minCost)
			{
				minCost = GX;
				minCostTheta = r;
			}
		}
	}

	// Build Affine transform 
	{
		// Affine without translation
		CvMat *RR = cvCreateMat(2, 2, CV_64FC1);
		double *r = RR->data.db;
		r[0] = cos(minCostTheta);
		r[1] = -sin(minCostTheta);
		r[2] = -r[1];
		r[3] = r[0];

		CvMat *t = cvCreateMat(2, 2, CV_64FC1);
		cvMatMul(RR, S1, t);
		cvInv(S2, RR);
		cvMatMul(RR, t, AA);
		cvReleaseMat(&t);
		cvReleaseMat(&RR);

		// Build translation matrix: shift = Mean1 - A * Mean0
		double *A = AA->data.db;
		shift->data.db[0] = xMean2 - A[0] * xMean1 - A[1] * yMean1;
		shift->data.db[1] = yMean2 - A[2] * xMean1 - A[3] * yMean1;
	}

	cvReleaseMat(&M2_3);
	cvReleaseMat(&M1_3);
	cvReleaseMat(&P2);
	cvReleaseMat(&P1);
	cvReleaseMat(&S2);
	cvReleaseMat(&S1);
	cvReleaseMat(&PP1);
	cvReleaseMat(&PP2);

	return minCostTheta;
}

void CmAffine::GetAffine2D(IN PointSet& pts0, IN PointSet& pts1)
{
	// Initial input data
	CvMat *Pts0 = CreatePointSetMat(pts0);
	CvMat *Pts1 = CreatePointSetMat(pts1);
	CvMat *AA = cvCreateMat(2, 2, CV_64FC1);
	CvMat *shift = cvCreateMat(2, 1, CV_64FC1);
	GetAffine2D(Pts0, Pts1, AA, shift);

	CvMat* showImg = cvCreateMat(800, 800, CV_8UC3);
	cvZero(showImg);
	DisplayPointSet(Pts0, CV_RGB(1.0,0,0), NULL, showImg);
	DisplayPointSet(Pts1, CV_RGB(0.0,1.0,0), NULL, showImg);
	CvMat *PP = cvCreateMat(2, Pts0->cols, CV_64FC1);
	//cvMatMul(AA, Pts0, PP);

	FindTranslated(Pts0, AA, shift, PP);
	

	DisplayPointSet(PP, CV_RGB(1.0,1.0,0), NULL, showImg);
	cvReleaseMat(&showImg);
	cvReleaseMat(&PP);

	cvReleaseMat(&Pts0);
	cvReleaseMat(&Pts1);
	cvReleaseMat(&shift);
	cvReleaseMat(&PP);
	cvReleaseMat(&AA);
}

void CmAffine::FindTranslated(IN CvMat *PtsB, IN CvMat *AA, IN CvMat *shift, OUT CvMat *PtsA)
{
	cvMatMul(AA, PtsB, PtsA);
	double *x = PtsA->data.db;
	double *y = (double*)(PtsA->data.ptr + PtsA->step);
#pragma omp parallel for
	for (int c = 0; c < PtsA->cols; c++)
	{
		x[c] += shift->data.db[0];
		y[c] += shift->data.db[1];
	}
}

CvMat* CmAffine::CreatePointSetMat(PointSet &pntSet)
{
	CvMat *pMat = cvCreateMat(2, pntSet.size(), CV_64FC1);
	double* x = pMat->data.db;
	double* y = x + pMat->cols;
	for (int i = 0, iSize = pMat->cols; i < iSize; i++)
	{
		x[i] = pntSet[i].x;
		y[i] = pntSet[i].y;
	}
	return pMat;
}

/************************************************************************/
/* pntSet: an 2xN CV_64FC1 matrix. N point positions                    */
/* saveName: if not NULL, the showing image will be saved to this path	*/
/* showImg: if input image is 800x800 then no memory is allocated		*/
/************************************************************************/
void CmAffine::DisplayPointSet(IN CvMat *pntSet, IN const CvScalar &color, 
							   IN const char* saveName, OUT CvMat *showBuffer)
{
	assert(pntSet->rows == 2);

	double *x = pntSet->data.db;
	double *y = (double*)(pntSet->data.ptr + pntSet->step);
	double minX = x[0], minY = y[0];
	double maxX = minX, maxY = minY;
	for (size_t i = 1, iEnd = pntSet->cols; i < iEnd; i++)
	{
		minX = min(minX, x[i]);
		minY = min(minY, y[i]);
		maxX = max(maxX, x[i]);
		maxY = max(maxY, y[i]);
	}

	printf("Display Point Set: %s : (%g, %g)->(%g, %g)\n", saveName, minX, minY, maxX, maxY);

	int dSize = 800;
	CvMat *dspImg = NULL;
	
	if (showBuffer == NULL)
	{
		dspImg = cvCreateMat(dSize, dSize, CV_8UC3);
		cvZero(dspImg);
	}
	else
	{
		assert(showBuffer->width == dSize && showBuffer->height == dSize);
		assert(CV_MAT_TYPE(showBuffer->type) == CV_8UC3);
		dspImg = showBuffer;
	}
	
	
	double scale = 700.0 / max((maxX - minX), (maxY - minY));

	double shiftX = -(maxX + minX) * scale / 2 + 400;
	double shiftY = -(maxY + minY) * scale / 2 + 400;
	for (int i = 0; i < pntSet->cols; i++)
	{
		CvPoint pnt = cvPoint(round(x[i] * scale + shiftX), round(y[i] * scale + shiftY));
		cvCircle(dspImg, pnt, 3, color);
	}

	if (saveName != NULL)
		cvSaveImage(saveName, dspImg);
	if (showBuffer == NULL)
		cvReleaseMat(&dspImg);
}
