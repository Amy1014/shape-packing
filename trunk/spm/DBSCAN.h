#ifndef DBSCAN_H
#define DBSCAN_H
#include "tsANN/include/ANN/ANN.h"
enum{DBSCAN_CORE_POINT,DBSCAN_BORDER_POINT,DBSCAN_NOISE_POINT};
typedef double MY_FLOAT;
class DBSCAN{
public:
	DBSCAN(MY_FLOAT * coord, long point_n,int dim,MY_FLOAT neighbor_radius_threshold,MY_FLOAT cluster_radius_threshold,int core_threshold,int& neighbor_n_min,int& neighbor_n_max);
	~DBSCAN();
	void recognize_core_points(int& neighbor_n_min,int& neighbor_n_max);
	void clustering_core_points();
	void clustering_non_core_points();
	void change_parameters(MY_FLOAT neighbor_radius_threshold,MY_FLOAT cluster_radius_threshold,int core_threshold,int& neighbor_n_min,int& neighbor_n_max);
	int inline get_internal_label(int index)
	{
		return internal_label_[index];
	}
	bool inline is_core_point(int index)
	{
		return is_core_point_[index];
	}
	void debug(MY_FLOAT x,MY_FLOAT y, MY_FLOAT z);
public:
	long point_n_;
	ANNpointArray dataPts_;
	ANNkd_tree* kdTree_;
	bool* is_core_point_;
	int* internal_label_;
	int label_n_;
	int dim_;
	//can be adjusted
	int core_threshold_;
	MY_FLOAT neighbor_radius_threshold_;//
	MY_FLOAT cluster_radius_threshold_;//
	//////////////////////////////////////////////////////////////////////////debug
	int debug_index;
};
#endif