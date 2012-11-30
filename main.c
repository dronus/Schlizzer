#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <array>
#include <math.h>
#include <exception>
#include <fstream>

using namespace std;

void error(const char* message){
	printf("%s\n",message);
	throw 0;
}

struct Vertex{
	float x,y,z;
	
	bool operator==(const Vertex& b) const{
		return this->x==b.x && this->y==b.y && this->z==b.z;
	}
};

bool operator!=(Vertex& a, Vertex& b){
	return !(a==b);
}

bool operator<(const Vertex& a, const Vertex& b)
{
	return 
		a.z< b.z             ||
		a.z==b.z && a.y< b.y ||
		a.z==b.z && a.y==b.y && a.x<b.x;
}

struct Triangle
{
	Vertex* vertices[3];
	std::vector<Triangle*> edges;
};

struct VertexIndex
{
	float value; 
	Triangle* triangle;
	char vertexIndex;
};

struct Segment
{
	Triangle* triangle;
	std::array<Vertex,2> vertices;
	std::array<Segment*,2> neighbours;
	int orderIndex;
	
	bool operator<(const Segment& other)const{
		return this->orderIndex<other.orderIndex;
	}
};


/*
// segment left-to-right order assuming segments dont cross
struct SegmentLeftOf
{
	bool operator() (const Segment* a, const Segment* b) const
	{
	
		Vertex a_bottom=a->vertices[0], a_top=a->vertices[1];
		if(a_bottom.y>a_top.y) std::swap(a_bottom, a_top);
		
		Vertex b0=b->vertices[0], b1=b->vertices[1];
		
		// select one of b's verticies inside a's height interval
		Vertex bm;
		if      (b0.y>=a_bottom.y && b0.y<=a_top.y) bm=b0;
		else if (b1.y>=a_bottom.y && b1.y<=a_top.y) bm=b1;
		else {
			printf("%f %f %f %f\n", a_bottom.y, a_top.y, b0.y, b1.y);
			assert(!"aua");
		}
		
		// project chosen b vertex onto a
		float t=(bm.y-a_bottom.y)/(a_top.y-a_bottom.y);
		float x=a_bottom.x+t*(a_top.x-a_bottom.x);
		
		// if projected point on a is left of bm, a is left of bm.
		return x<bm.x;
	};
};
*/

struct Layer
{
	float z;
	std::vector<Triangle*> triangles;
	std::vector<Segment> segments;
};

bool VertexIndexLessThan(const VertexIndex& a, const VertexIndex& b)
{
	return a.value<b.value;	
};

class VertexLessThan
{
	public: 
	bool operator() (const Vertex& a, const Vertex& b) const
	{
		return 
			a.z< b.z             ||
			a.z==b.z && a.y< b.y ||
			a.z==b.z && a.y==b.y && a.x<b.x;
	}
};


float distance(Vertex& a, Vertex& b){
	float dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}



/*template <class T> void swap(T& a, T& b)
{
	T tmp=a;
	a=b;
	b=tmp;
}*/

void sortTriangleVertices(Triangle& t)
{
	Vertex** vs=t.vertices;
	if(vs[0]->z > vs[1]->z) std::swap(vs[0],vs[1]);
	if(vs[0]->z > vs[2]->z) std::swap(vs[0],vs[2]);
	if(vs[1]->z > vs[2]->z) std::swap(vs[1],vs[2]);
}

void loadStl(const char* filename, std::vector<Vertex>& vertices, std::vector<Triangle>& triangles)
{
	std::map<Vertex,int,VertexLessThan> uniqueVertices;
	std::vector<int> indices;

	printf("Loading %s...\n",filename);
	FILE* file=fopen(filename,"r");	
	char line[256];
	while(!feof(file)){
		if(fgets(line, sizeof(line), file)){
			Vertex p;
			int found=sscanf(line," vertex %e %e %e",&p.x,&p.y,&p.z);
			if(found) {
				int uniqueVertex;
				if(uniqueVertices.count(p)==1) {
					uniqueVertex=uniqueVertices[p];					
				}else{
					vertices.push_back(p);
					uniqueVertex=vertices.size()-1;
					uniqueVertices[p]=uniqueVertex;
				}
				indices.push_back(uniqueVertex);
			}
		}
	}
	fclose(file);

	if(indices.size()%3 != 0)
		assert(!"indices count not mod 3");

	for(int i=0; i<indices.size(); i+=3)
	{
		Triangle t;
		for(int j=0; j<3; j++)
			t.vertices[j]=&vertices[indices[i+j]];
		sortTriangleVertices(t);
		assert(t.vertices[0]->z<=t.vertices[1]->z && t.vertices[1]->z<=t.vertices[2]->z);
		triangles.push_back(t);	
	}

	// compute edges
/*	std::map<std::pair<Vertex*,Vertex*>,Triangle*> uniqueEdges;
	for(int i=0; i<triangles.size(); i++)
	{
		Triangle& t=triangles[i];
		for(int j=0; j<3; j++){
			Vertex* v0=t.vertices[j        ];
			Vertex* v1=t.vertices[(j+1) % 3];

			if(v0>v1) std::swap(v0,v1);
			std::pair<Vertex*,Vertex*> edge=std::pair<Vertex*,Vertex*>(v0,v1);
			if(uniqueEdges.count(edge)==1) {
				Triangle* t2=uniqueEdges[edge];
				t.edges.push_back(t2);
				t2->edges.push_back(&t);
			}else{
				uniqueEdges[edge]=&t;
			}
		}
	}

	for(int i=0; i<triangles.size(); i++)
		assert(triangles[i].edges.size()==3);
*/	
	printf("Loading complete: %d vertices read, %d unique, %d triangles.\n",(int)indices.size(),(int)vertices.size(),(int)triangles.size());
}

void buildLayers(std::vector<Triangle>& triangles, std::vector<Layer>& layers)
{
	
	std::vector<VertexIndex> by_z;
	for(int i=0; i<triangles.size(); i++)
		for(char j=0; j<3; j++){
			VertexIndex vi={
				triangles[i].vertices[j]->z,
				&triangles[i],
				j
			};
			by_z.push_back(vi);
		}
	std::sort(by_z.begin(), by_z.end(),VertexIndexLessThan);

	std::map<Triangle*,int> triangleRuns;

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
			for(std::map<Triangle*,int>::iterator j=triangleRuns.begin(); j!=triangleRuns.end(); ++j)
				layer.triangles.push_back(j->first);
			layers.push_back(layer);
			next_layer_z+=layer_height;
		}

		Triangle* triangle=by_z[i].triangle;
		int count;
		if(triangleRuns.count(triangle)==0) count=1;
		else                                count=triangleRuns[triangle]+1;
		triangleRuns[triangle]=count;
		if(count==3) triangleRuns.erase(triangle);		
	}				
	assert(triangleRuns.size()==0);

	printf("Layers: %d\n",(int)layers.size());
}

Vertex computeIntersection(Vertex& a, Vertex& b, float z)
{
	float t=(z-a.z)/(b.z-a.z);
	Vertex intersection={
		a.x+t*(b.x-a.x),
		a.y+t*(b.y-a.y),
		z
	};

	return intersection;
}

std::array<Vertex,2> computeSegment(Triangle& t, float z)
{
	Vertex** vs=t.vertices;
	std::array<Vertex,2> segment;

	// triangle vertices are always ordered by z

	// ensure the triangles are correctly assigned to the layers
	assert(z>=vs[0]->z);
	assert(z<=vs[2]->z);

	// so we just need to check the second vertex to decide which edges
	// get intersected.	
	if(z<vs[1]->z){
		segment[0]=computeIntersection(*vs[0],*vs[1],z);
		segment[1]=computeIntersection(*vs[0],*vs[2],z);
	}else{
		segment[0]=computeIntersection(*vs[1],*vs[2],z);
		segment[1]=computeIntersection(*vs[0],*vs[2],z);
	}
	
	return segment;
}

void fill(Layer& layer, std::map<Vertex,std::vector<Segment*>,VertexLessThan> segmentsByVertex)
{
	float grid_spacing=0.5f;
	
	std::set<Vertex,VertexLessThan> sweepVertices;
	for(std::map<Vertex,std::vector<Segment*>>::iterator i=segmentsByVertex.begin(); i!=segmentsByVertex.end(); ++i){
		assert(i->first.z==layer.z);
		sweepVertices.insert(i->first);
	}
	
	std::vector<Segment> infill;
//	infill.reserve(1000);
	
	std::set<Segment*> sweepHeap;	
	float y=sweepVertices.begin()->y+grid_spacing;
	for(std::set<Vertex>::iterator i=sweepVertices.begin(); i!=sweepVertices.end(); ++i)
	{
		const Vertex& v=*i;
		
		while(v.y>y){
			std::vector<Vertex> intersections;
			for(std::set<Segment*>::iterator j=sweepHeap.begin(); j!=sweepHeap.end(); ++j){
				Segment& s=**j;
				Vertex& a=s.vertices[0], &b=s.vertices[1];
				float t=(y-a.y)/(b.y-a.y);
				if(t>=0 && t<1){
					Vertex intersection={
						a.x+t*(b.x-a.x),
						y,  // exact intersection 
						a.z // z const in layer
					};
					intersections.push_back(intersection);
				}

			}
//			assert(intersections.size() % 2 == 0);

			std::sort(intersections.begin(), intersections.end());

//			printf("\t\tSweep at y: %f, x's: ",y);
			for(int j=0; j<intersections.size(); j++){
				if(j%2==1) {
					Segment s={
						NULL,
						{intersections[j-1],intersections[j]},
						{NULL,NULL},
						(j+1)*0xFFF
					};
					infill.push_back(s);
				}
			}
//				printf("%f ",intersections[j].x);
//			printf("\n");	
			y+=grid_spacing;
		}
		
		std::vector<Segment*>& ss=segmentsByVertex[v];
		char segments_in_heap=0;
		char segment_index=-1;

		// ignore non manifold unconnected segments
		if(ss.size()!=2) continue;
		
		for(int j=0; j<ss.size(); j++){
			assert(ss[j]->vertices[0]==v || ss[j]->vertices[1]==v);
			if(sweepHeap.count(ss[j])==1) {
				segments_in_heap++;
				segment_index=j;
			}
		}

		if(segments_in_heap==0){
			// a new perimeter is encountered. add it to the heap.
			sweepHeap.insert(ss[0]);
			sweepHeap.insert(ss[1]);
		}else if(segments_in_heap==1){
			// an existing perimeter continues, segment_index must be it's index.
			sweepHeap.erase (ss[  segment_index]);
			sweepHeap.insert(ss[1-segment_index]);
		}else if(segments_in_heap==2){
			// an existing perimeter ends. 
			sweepHeap.erase (ss[0]);
			sweepHeap.erase (ss[1]);
		}		
		// DIRTY: commented out to ignore non manifolds
		//else
		//	assert(segments_in_heap<=2); // will fail. 
			
		
		
	}

	// DIRTY: commented out to ignore non manifolds
	//	assert(sweepHeap.size()==0);
	
	layer.segments.insert(layer.segments.end(),infill.begin(),infill.end());
}


void buildSegments(Layer& layer)
{

	// generate segments by intersecting the triangles touching this layer
	for(int i=0; i<layer.triangles.size(); i++)
	{
		Triangle* t=layer.triangles[i];
		Segment s={
			t,
			computeSegment(*t,layer.z),
			{NULL,NULL},
			-1
		};
		if(s.vertices[0]!=s.vertices[1]) {
		
	//		s.orderIndex=i; // TODO do better ordering

			layer.segments.push_back(s);
		}
	}
	
	// unify segment end vertices
	std::map<Vertex,std::vector<Segment*>,VertexLessThan> segmentsByVertex;
	for(int i=0; i<layer.segments.size(); i++)
	{
		Segment& s=layer.segments[i];
		for(int j=0; j<2; j++){
			Vertex& v=s.vertices[j];
			segmentsByVertex[v].push_back(&s);
		}
	}
		
	// link segments by neighbour pointers using the unique vertex map 
	for(std::map<Vertex,std::vector<Segment*>>::iterator i=segmentsByVertex.begin(); i!=segmentsByVertex.end(); ++i)
	{
		std::vector<Segment*>& ss=i->second;
		
		if(ss.size()==1) error("Unconnected segment");
//		if(ss.size()>2 ) error("Non manifold segment");
		
		Vertex v=i->first;
		
		// as we don't know the direction of each segment in the final trajectory,
		// we just link them in the same order as they list their vertices.
		// use two indices for the corresponding neighbour pointers
		// TODO maybe we should make this simpler and just use the first free neighbour pointer,
		// however errors are harder to track than.
		int index0, index1;
		if       (ss[0]->vertices[0]==v) index0=1;
		else if  (ss[0]->vertices[1]==v) index0=0;
		else     assert(!"evil index0");

		if       (ss[1]->vertices[0]==v) index1=1;
		else if  (ss[1]->vertices[1]==v) index1=0;
		else     assert(!"evil index1");
		
		// now index0, index1 should point to a free end of the segment
		assert(ss[0]->neighbours[index0]==NULL);
		assert(ss[1]->neighbours[index1]==NULL);
		
		// finally link both segments
		ss[0]->neighbours[index0]=ss[1];
		ss[1]->neighbours[index1]=ss[0];

		/*		
		//alternative one step neighbour linking 
		if(segmentsByVertex.count(v)==1) {
			Segment* neighbour=segmentsByVertex[v];
			assert(s.neighbours[j]==NULL);
			s.neighbours[j]=neighbour;
			assert(neighbour->neighbours[(j+1)%2]==NULL);			
			neighbour->neighbours[(j+1)%2]=&s;
		}else{
			segmentsByVertex[v]=&s;
		}
		*/
	}

	// check for dangling segments (caused by disconnected triangles)
/*	for(int i=0; i<layer.segments.size(); i++)
		for(int j=0; j<2; j++)
			if(layer.segments[i].neighbours[j]==NULL) {
				printf("Unconnected segment: %d %d\n",i,j);
				throw 0;
			}
*/	
		
	// now order the segments into consecutive loops. 
	int loops=0;
	int orderIndex=0;
	for(int i=0; i<layer.segments.size(); i++){
		Segment& segment=layer.segments[i];

		// only handle new loops
		if(segment.orderIndex!=-1) continue;
				
		// collect a loop
		Segment* s2=&segment;
		while(true){
			s2->orderIndex=orderIndex++;
			// DIRTY: check for NULL neighbours resulting from non manifolds
			if     (s2->neighbours[0] != NULL && s2->neighbours[0]->orderIndex==-1)
				s2=s2->neighbours[0];
			else if(s2->neighbours[1] != NULL && s2->neighbours[1]->orderIndex==-1)
				s2=s2->neighbours[1];
			else break;
		};
		
		// the loop should be closed:
		// DIRTY: ignore check to accept non manifolds		
//		assert(s2->neighbours[0]==&segment || s2->neighbours[1]==&segment);
		
		loops++;
	}
	
	// printf("\tTriangles: %d, segments: %d, vertices: %d, loops: %d\n",(int)layer.triangles.size(),(int)layer.segments.size(),(int)segmentsByVertex.size(),loops);

	fill(layer,segmentsByVertex);
	// caution: the neighbour[..] and other segment pointers are invalid now!	

	std::sort(layer.segments.begin(), layer.segments.end());
}


void saveGcode(const char* filename, std::vector<Layer>& layers)
{


	printf("Saving Gcode...\n");
	FILE* file=fopen(filename,"w");	
/*
G92 E0
G1 Z8.000 F7800.000
G1 X90.263 Y90.263
G1 F1800.000 E1.00000
G1 X109.737 Y90.263 F1800.000 E1.57267
G1 X109.737 Y109.737 E2.14534
G1 X90.263 Y109.737 E2.71801
G1 X90.263 Y90.787 E3.27524
G1 F1800.000 E2.27524

*/
	float extrusionFactor=0.0294f;
	float retractLength=1.f;

	Vertex offset={75.f, 75.f, 0.f};
	Vertex position={0,0,0};

	int travels=0, extrusions=0;
	int travelsSkipped=0, extrusionsSkipped=0;
	float travelled=0, extruded=0;

	float skipDistance=.01;

	for(int i=0; i<layers.size(); i++){
		Layer& l=layers[i];
		fprintf(file, "G92 E0\n");
		fprintf(file, "G1 Z%f F7800.000\n",l.z);
		float extrusion=retractLength;
		fprintf(file, "G1 F1800.000 E%f\n",extrusion); // retraction restart
		for(int j=0; j<l.segments.size(); j++){
			Vertex& v0=l.segments[j].vertices[0];
			Vertex& v1=l.segments[j].vertices[1];
			assert(v0.z==l.z);
			assert(v1.z==l.z);
			
			if(distance(v1,position)<distance(v0,position))
				swap(v0,v1);
			
			if(distance(v0,position)>skipDistance){
				fprintf(file,"G1 X%f Y%f\n",v0.x+offset.x,v0.y+offset.y);
				travels++;
				travelled+=distance(v0,position);
			}else
				travelsSkipped++;
			position=v0;
			extrusion+=extrusionFactor*distance(v0,v1);
			if(distance(v1,position)>skipDistance){
				fprintf(file,"G1 X%f Y%f E%f\n",v1.x+offset.x,v1.y+offset.y,extrusion);
				extrusions++;
				extruded+=distance(v1,position);
			}else
				extrusionsSkipped++;
			position=v1;
		}
		extrusion-=retractLength;
		fprintf(file, "G1 F1800.000 E%f\n",extrusion); // retraction restart		
	}	
	printf("Saving complete. %d bytes written. %d travels %.0f mm, %d extrusions %.0f mm, %d travel skips, %d extrusion skips\n",ftell(file),travels, travelled, extrusions, extruded, travelsSkipped, extrusionsSkipped);
}

int main(int argc, const char** argv)
{
	if(argc!=3) {
		printf("Usage: %s <.stl file> <.gcode file>\n",argv[0]);
		return 1;
	}

	std::vector<Vertex>   vertices;
	std::vector<Triangle> triangles;

	loadStl(argv[1], vertices, triangles);

	std::vector<Layer> layers;
	buildLayers(triangles, layers);

	for(int i=0; i<layers.size(); i++)
		buildSegments(layers[i]);

	/*for(int i=0; i<layers.size(); i++){
		Layer& l=layers[i];
		printf("Layer %d: z: %f, triangles: %d, segments: %d\n",i,l.z, (int)l.triangles.size(),(int)l.segments.size());
		for(int j=0; j<l.segments.size(); j++){
			Vertex& v0=l.segments[j].vertices[0];
			Vertex& v1=l.segments[j].vertices[1];
			printf("\tSegment (%f,%f,%f) to (%f,%f,%f)\n",v0.x,v0.y,v0.z,v1.x,v1.y,v1.z);
		}
	}*/
	
	saveGcode(argv[2],layers);
	//std::map<int,std::vector<int>> trianglesByVertex;

	return 0;
}



