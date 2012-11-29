#include "stdio.h"
#include <vector>
#include <map>
#include <algorithm>


struct Vertex{
	float x,y,z;
};

struct Triangle
{
	Vertex vertices[3];
};

struct VertexIndex
{
	float value; 
	int triangleIndex;
	char vertexIndex;
};

struct Segment
{
	int triangleIndex;
	Vertex vertices[2];
};

struct Layer
{
	float z;
	std::vector<int> triangleIndices;
	std::vector<Segment> segments;
};

bool VertexIndexLessThan(const VertexIndex& a, const VertexIndex& b)
{
	return a.value<b.value;	
};

void loadStl(const char* filename, std::vector<Triangle>& triangles)
{
	printf("Loading %s...\n",filename);
	FILE* file=fopen(filename,"r");	
	int i=0;
	Triangle t;
	char line[256];
	while(!feof(file)){
		if(fgets(line, sizeof(line), file)){
			Vertex p;	
			int found=sscanf(line," vertex %e %e %e",&p.x,&p.y,&p.z);
			if(found) {
				t.vertices[i%3]=p;
				i++;
				if(i%3==0) triangles.push_back(t);				
			}
		}
	}
	fclose(file);
	printf("Loading complete, got %d triangles.\n",(int)triangles.size());
}

int main(int argc, const char** argv)
{
	if(argc!=2) {
		printf("Usage: a.out <filename>\n");
		return 1;
	}

	std::vector<Triangle> triangles;

	loadStl(argv[1],triangles);

	std::vector<VertexIndex> by_z;
	for(int i=0; i<triangles.size(); i++)
		for(int j=0; j<3; j++){
			VertexIndex vi={
				triangles[i].vertices[j].z,
				i,
				j
			};
			by_z.push_back(vi);
		}

//	for(int i=0; i<by_z.size(); i++)
//		printf("Z Value: %f\n",by_z[i].value);	


	printf("Sorting %d vertices.\n",(int)by_z.size());

	std::sort(by_z.begin(), by_z.end(),VertexIndexLessThan);

//	for(int i=0; i<by_z.size(); i++)
//		printf("Z Value: %f\n",by_z[i].value);	

	std::vector<Layer> layers;
	std::map<int,int> triangleRuns;

	float layer_height=0.3f;
	float min_z=by_z.front().value;
	float max_z=by_z.back().value;
	printf("Slicing from %f to %f, layers: %d\n",min_z, max_z, (int)((max_z-min_z)/layer_height));
	float next_layer_z=min_z+layer_height;
	for(int i=0; i<by_z.size(); i++)
	{
		float z=by_z[i].value;
		while(z>next_layer_z) {
			Layer layer;
			layer.z=next_layer_z;
			for(std::map<int,int>::iterator j=triangleRuns.begin(); j!=triangleRuns.end(); ++j)
				layer.triangleIndices.push_back(j->first);
			layers.push_back(layer);
			next_layer_z+=layer_height;
		}

		int index=by_z[i].triangleIndex;
		int count;
		if(triangleRuns.count(index)==0) count=1;
		else                             count=triangleRuns[index]+1;
		triangleRuns[index]=count;
		if(count==3) triangleRuns.erase(index);		
	}				

	printf("Layers: %d, triangle runs left: %d\n",(int)layers.size(),(int)triangleRuns.size());

	for(int i=0; i<layers.size(); i++)
		printf("Layer %d: z: %f, triangles: %d\n",i,layers[i].z, (int)layers[i].triangleIndices.size());

	return 0;
}
