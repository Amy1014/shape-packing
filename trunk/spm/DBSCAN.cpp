#include "DBSCAN.h"
#include <set>
#include <list>
DBSCAN::DBSCAN(MY_FLOAT * coord, long point_n,int dim,MY_FLOAT neighbor_radius_threshold,MY_FLOAT cluster_radius_threshold,int core_threshold,int& neighbor_n_min,int& neighbor_n_max)
{
	point_n_=point_n;
	dim_=dim;
	neighbor_radius_threshold_=neighbor_radius_threshold;
	cluster_radius_threshold_=cluster_radius_threshold;
	core_threshold_=core_threshold;
	label_n_=0;
	dataPts_ = annAllocPts(point_n, dim); // allocate data points
	memcpy(dataPts_[0],coord,point_n*dim*sizeof(MY_FLOAT));
	is_core_point_=new bool[point_n];
	internal_label_=new int[point_n];
	kdTree_ = new ANNkd_tree(dataPts_, point_n, dim);
	recognize_core_points(neighbor_n_min,neighbor_n_max);
	clustering_core_points();
	clustering_non_core_points();
}
void DBSCAN::change_parameters(MY_FLOAT neighbor_radius_threshold,MY_FLOAT cluster_radius_threshold,int core_threshold,int& neighbor_n_min,int& neighbor_n_max)
{
	neighbor_radius_threshold_=neighbor_radius_threshold;
	cluster_radius_threshold_=cluster_radius_threshold;
	core_threshold_=core_threshold;
	recognize_core_points(neighbor_n_min,neighbor_n_max);
	clustering_core_points();
	clustering_non_core_points();
}
DBSCAN::~DBSCAN()
{
	annDeallocPts(dataPts_);
	if(is_core_point_!=NULL)
		delete [] is_core_point_;
	if(internal_label_!=NULL)
		delete [] internal_label_;
	delete kdTree_;
	annClose(); // done with ANN
}
void DBSCAN::recognize_core_points(int& neighbor_n_min,int& neighbor_n_max)
{
	neighbor_n_max=-1;
	neighbor_n_min=1e8;
	MY_FLOAT squared_radius_threshold=neighbor_radius_threshold_*neighbor_radius_threshold_;
	for (int i=0;i<point_n_;++i)
	{
		int pn_in_sphere=kdTree_->annkFRSearch(dataPts_[i],squared_radius_threshold,0)-1;
		if (neighbor_n_max<pn_in_sphere)neighbor_n_max=pn_in_sphere;
		if (neighbor_n_min>pn_in_sphere)neighbor_n_min=pn_in_sphere;		
		is_core_point_[i]=(pn_in_sphere>=core_threshold_);
	}
}
void DBSCAN::clustering_core_points()
{
	bool* checked=new bool[point_n_];
	int* cache_index=new int[100];
	MY_FLOAT* cache_dist=new MY_FLOAT[100];
	int cache_size=100;
	for (int i=0;i<point_n_;++i)
	{
		internal_label_[i]=-1;
		checked[i]=false;
	}
	int current_label=0;
	std::list<int> temp_stack;
	MY_FLOAT squared_radius_threshold=cluster_radius_threshold_*cluster_radius_threshold_;
	for(int current_index=0;current_index<point_n_;++current_index)//find a unchecked point
	{
		if (!checked[current_index]&&is_core_point(current_index))
		{
			internal_label_[current_index]=current_label;
			checked[current_index]=true;
			temp_stack.push_back(current_index);
			while (!temp_stack.empty())
			{
				int current_point=temp_stack.front();
				temp_stack.pop_front();
				int pn_in_sphere=kdTree_->annkFRSearch(dataPts_[current_point],squared_radius_threshold,0);
				if (pn_in_sphere>cache_size)
				{
					cache_size=pn_in_sphere*2;
					delete []cache_index;
					cache_index=new int[cache_size];
					delete []cache_dist;
					cache_dist=new MY_FLOAT[cache_size];
				}
				kdTree_->annkFRSearch(dataPts_[current_point],squared_radius_threshold,pn_in_sphere,cache_index,cache_dist);
				for (int i=1;i<pn_in_sphere;++i)
				{
					if (!checked[cache_index[i]]&&is_core_point(cache_index[i]))
					{
						internal_label_[cache_index[i]]=current_label;
						checked[cache_index[i]]=true;
						temp_stack.push_back(cache_index[i]);
					}					
				}
			}
			++current_label;
		}
	}
	label_n_=current_label;
	delete[] cache_index;
	delete[] cache_dist;
	delete[] checked;
}
void DBSCAN::clustering_non_core_points()//assign it to its nearest core point
{
	int* cache_index=new int[100];
	MY_FLOAT* cache_dist=new MY_FLOAT[100];
	int cache_size=100;
	int initial_query_n=20;
	for (int i=0;i<point_n_;++i)
	{
		if(internal_label_[i]<0)
		{
			bool need_continue=true;
			int current_query_n=initial_query_n;
			while(need_continue)
			{
				kdTree_->annkSearch(dataPts_[i],current_query_n,cache_index,cache_dist);
				int j=1;
				if (current_query_n>initial_query_n)j=current_query_n-initial_query_n;
				for (;j<current_query_n;++j)
				{
					if (is_core_point(cache_index[j])&&internal_label_[cache_index[j]]>=0)
					{
						internal_label_[i]=internal_label_[cache_index[j]];
						need_continue=false;
						break;
					}
				}
				if (need_continue)
				{
					current_query_n+=initial_query_n;
					if (current_query_n>cache_size)
					{
						cache_size=current_query_n*2;
						delete []cache_index;
						cache_index=new int[cache_size];
						delete []cache_dist;
						cache_dist=new MY_FLOAT[cache_size];
					}
				}
			}

		}
	}
	delete[] cache_index;
	delete[] cache_dist;
}
void DBSCAN::debug(MY_FLOAT x,MY_FLOAT y, MY_FLOAT z)
{
	ANNcoord p[]={x,y,z};
	ANNidx index[]={-1};
	ANNdist dist[]={0};
	kdTree_->annkSearch(p,1,index,dist);
	debug_index=index[0];
	MY_FLOAT squared_radius_threshold=neighbor_radius_threshold_*neighbor_radius_threshold_;
	int pn_in_sphere=kdTree_->annkFRSearch(dataPts_[debug_index],squared_radius_threshold,0)-1;
}