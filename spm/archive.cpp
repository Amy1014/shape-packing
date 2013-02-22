#include "archive.h"

namespace Geex{

	//specialized IO for each type
	//text 
	fstream& operator>>(fstream& infile, vec3& v)
	{
		infile>>v.x>>v.y>>v.z;
		return infile;
	}
	fstream& operator<<(fstream& outfile, const vec3& v)
	{
		outfile<<v.x<<v.y<<v.z;
		return outfile;
	}
	fstream& operator>>(fstream& infile, Segment_3& s)
	{
		double sx, sy, sz, tx, ty, tz;
		infile>>sx>>sy>>sz;
		infile>>tx>>ty>>tz;
		s = Segment_3(Point_3(sx, sy, sz), Point_3(tx, ty, tz));
		return infile;
	}
	fstream& operator<<(fstream& outfile, const Segment_3& s)
	{
		Point_3 src = s.source(), tgt = s.target();
		outfile<<src.x()<<src.y()<<src.z();
		outfile<<tgt.x()<<tgt.y()<<tgt.z();
		return outfile;
	}
	//binary
	fstream& Archive::write(const vec3& v)
	{
		fs.write((const char*)&v.x, sizeof(double));
		fs.write((const char*)&v.y, sizeof(double));
		fs.write((const char*)&v.z, sizeof(double));
		return fs;
	}
	fstream& Archive::read(vec3& v)
	{
		fs.read((char*)&v.x, sizeof(double));
		fs.read((char*)&v.y, sizeof(double));
		fs.read((char*)&v.z, sizeof(double));
		return fs;
	}
	fstream& Archive::write(const Point_2& p)
	{
		fs.write((const char*)&p.x(), sizeof(double));
		fs.write((const char*)&p.y(), sizeof(double));
		return fs;
	}
	fstream& Archive::read(Point_2 &p)
	{
		double x, y;
		fs.read((char*)&x, sizeof(double));
		fs.read((char*)&y, sizeof(double));
		p = Point_2(x, y);
		return fs;
	}
	fstream& Archive::write(const Point_3& p)
	{
		fs.write((const char*)&p.x(), sizeof(double));
		fs.write((const char*)&p.y(), sizeof(double));
		fs.write((const char*)&p.z(), sizeof(double));
		return fs;
	}
	fstream& Archive::read(Point_3 &p)
	{
		double x, y, z;
		fs.read((char*)&x, sizeof(double));
		fs.read((char*)&y, sizeof(double));
		fs.read((char*)&z, sizeof(double));
		p = Point_3(x, y, z);
		return fs;
	}
	fstream& Archive::write(const Segment_3& s)
	{
		write(s.source());
		write(s.target());
		return fs;
	}
	fstream& Archive::read(Segment_3& s)
	{
		Point_3 src, tgt;
		read(src);
		read(tgt);
		s = Segment_3(src, tgt);
		return fs;
	}
	fstream& Archive::write(const Polygon_2& pgn)
	{
		unsigned int size(pgn.size());
		fs.write((const char*)&size, sizeof(unsigned int));
		for (unsigned int i = 0; i < size; i++)
			write(pgn.vertex(i));
		return fs;
	}
	fstream& Archive::read(Polygon_2& pgn)
	{
		pgn.clear();
		unsigned int size;
		fs.read((char*)&size, sizeof(unsigned int));
		for (unsigned int i = 0; i < size; i++)
		{
			Point_2 p;
			read(p);
			pgn.push_back(p);
		}
		return fs;
	}
	fstream& Archive::write(const vector<bool>& bitword)
	{
		unsigned int size = bitword.size();
		fs.write((const char*)&size, sizeof(unsigned int));
		char t(1);
		char f(0);
		for (unsigned int i = 0; i < size; i++)
		{
			if (bitword[i])
				fs.write((const char*)&t, sizeof(char));
			else
				fs.write((const char*)&f, sizeof(char));
		}
		return fs;
	}
	fstream& Archive::read(vector<bool>& bitword)
	{
		bitword.clear();
		unsigned int size;
		fs.read((char*)&size, sizeof(unsigned int));
		for (unsigned int i = 0; i < size; i++)
		{
			char v;
			fs.read((char*)&v, sizeof(char));
			if (v)
				bitword.push_back(true);
			else
				bitword.push_back(false);
		}
		return fs;
	}
	
	Archive::Archive(const string& filename, type t) 
			//fs(filename.c_str(), fstream::out | /*fstream::in |*/ fstream::binary) 
			{
				if (t == STORE)
					fs.open(filename.c_str(), fstream::out | fstream::binary);
				else
				{
					fs.open(filename.c_str(), fstream::in | fstream::binary);
					if (fs.fail())
						std::cout<<"File opening for reading failed.\n";
				}
			}

	Archive::~Archive()
	{
		fs.close();
	}
	void Archive::close()
	{
		fs.close();
	}


}