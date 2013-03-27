#pragma once

/************************************************************************/
/* Source code implementing:[09ICCV]An Algebraic Approach to Affine		*/
/*     Registration of Point Sets										*/
/************************************************************************/
#include <cv.h>
#include <Geex/basics/types.h>
#define IN 
#define OUT

using namespace cv;
struct CmAffine
{
	typedef std::vector<cv::Point2d> PointSet; 

public:

	/************************************************************************/
	/* srcP64F, dst64F: 2 by N matrix composed by (x, y) of N points		*/
	/* S64F:    a 2 by 2 matrix, dst = S64F * srcP64F;						*/
	/************************************************************************/
	static void Reduction(IN CvMat *srcP64F, OUT CvMat *S64F, OUT CvMat *dst64F);

	/************************************************************************/
	/* Pts1 is 2 by N1 matrix												*/
	/* Pts2 is 2 by N2 matrix												*/
	/* AA is 2x2 matrix and shift is 2x1 matrix								*/
	/* AA is the affine transform of the two point set that AA*Pts1+shift 	*/
	/*       best approximate Pts2.											*/
	/* Return the rotation angle  										    */
	/************************************************************************/
	static double GetAffine2D(IN CvMat *Pts1, IN CvMat *Pts2, OUT CvMat *AA, OUT CvMat *shift);

	/************************************************************************/
	/*                                                                      */
	/************************************************************************/
	static void GetAffine2D(IN PointSet& pts0, IN PointSet& pts1);

	/************************************************************************/
	/* Use translation AA and shift to deform positions in PtsB             */
	/* PtsB: 2xN matrix for points before translation						*/
	/* PtsA: 2xN matrix for points after translation						*/
	/************************************************************************/
	static void FindTranslated(IN CvMat *PtsB, IN CvMat *AA, IN CvMat *shift, OUT CvMat *PtsA);

	static bool ExtractBoundary(const CvMat *img8U, PointSet &pntSet);
	
	static bool GetSampling(PointSet& pntSet, int samplenum);
	static CvMat *CreatePointSetMat(PointSet &pntSet);


	/************************************************************************/
	/* pntSet: an 2xN CV_64FC1 matrix. N point positions                    */
	/* saveName: if not NULL, the showing image will be saved to this path	*/
	/* showImg: if input image is 800x800 then no memmory is alloced		*/
	/************************************************************************/
	static void DisplayPointSet(IN CvMat *pntSet, IN const CvScalar &color,
		IN const char* saveName = NULL, OUT CvMat *showBuffer = NULL);
};

double sqrDist(cv::Point2d& p, cv::Point2d& q);
double round(double number);