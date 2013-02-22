#pragma  once
#include "spm_cgal.h"
#include <fstream>
//#include <iostream>
#include <vector>
#include <string>

using namespace std;
namespace Geex{

	template<typename T> fstream& operator>>(fstream& infile, vector<T>& v)
	{
		unsigned int size;
		infile>>size;
		v.resize(size);
		for (unsigned int i = 0; i < size; i++)
			infile>>v[i];
		return infile;
	}
	template<typename T> fstream& operator<<(fstream& outfile, const vector<T>& v)
	{
		outfile<<v.size();
		for (unsigned int i = 0; i < v.size(); i++)
			outfile<<v[i];
		return outfile;
	}

class Archive
{
public:
	typedef enum {STORE, LOAD} type;
	Archive(const string& filename, type t = STORE);

	// binary IO
	template<typename T> fstream& write(const T& v)
	{
		fs.write((const char*)&v, sizeof(T));
		return fs;
	}
	template<typename T> fstream& read(T& v)
	{
		fs.read((char*)&v, sizeof(T));
		return fs;
	}
	fstream& write(const vec3& v);
	fstream& read(vec3& v);
	fstream& write(const Point_2& p);
	fstream& read(Point_2& p);
	fstream& write(const Point_3& p);
	fstream& read(Point_3& p);
	fstream& write(const Segment_3& s);
	fstream& read(Segment_3& s);
	fstream& write(const Polygon_2& pgn);
	fstream& read(Polygon_2& pgn);
	fstream& write(const vector<bool>& bitword);
	fstream& read(vector<bool>& bitword);
	
	template <typename T> fstream& write(const vector<T>& v)
	{
		unsigned int size(v.size());
		fs.write((const char*)&size, sizeof(unsigned int));
		for (unsigned int i = 0; i < size; i++)
			write(v[i]);
		return fs;
	}
	template <typename T> fstream& read(vector<T>& v)
	{
		unsigned int size;
		fs.read((char*)&size, sizeof(unsigned int));
		v.resize(size);
		for (unsigned int i = 0; i < size; i++)
			read(v[i]);
		return fs;
	}
	template<typename T> void append_var( const T& v )
	{
		fs.write((const char*)&v, sizeof(T));
	}
	template<typename T> void append_vector(const vector<T>& v)
	{
		//std::cout<<"Vector Len: "<<v.size()<<std::endl;
		//fs<<unsigned int(v.size())<<std::endl;
		//append_var(v.size());
		unsigned int size(v.size());
		fs.write((char*)&size, sizeof(unsigned int));
		for (unsigned int i = 0; i < v.size(); i++)
			write(v[i]);
	}
	template<typename T> void restore_var(T &v)
	{
		fs.read((char*)&v, sizeof(T));
	}
	template<typename T> void restore_vector(vector<T> &v)
	{
		unsigned int size;
		fs.read((char*)&size, sizeof(unsigned int));
		//std::cout<<"Vector length: "<<size<<std::endl;
		v.resize(size);
		for (unsigned int i = 0; i < size; i++)
			read(v[i]);
	}
	void close();
	~Archive();
private:
	fstream fs;
};

}