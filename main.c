/*
 * Schlizzer
 *
 * A early experimental Slic3r replacement in C++0x
 *
 * Usage: Sschlizzer <input.stl> <output.gcode>
 * 
 * This program loads a given .stl (stereolithography data, actually triangle data) file 
 * and generates a .gcode (RepRap machine instructions) file that can be printed on a RepRap
 * machine.
 *
 * This process is mostly referred as "slicing", as the printers make objects layer by layer,
 * so the 3d model has to be virtually sliced into thin layers. This however is only the first
 * step, after which perimeters and a crosshatch filling pattern are computed for every layer.
 *
 * The .stl file should describe a clean mesh that is a 2-manifold. That means:
 * * every triangle should touch exactly three triangles along its three edges
 * * on every edge, two triangles share exactly two vertices
 * * every triangle has some positive surface area
 *  a triangle may touch further triangles at its vertices.
 *
 *  currently, many checks are disabled to accept objects breaking any of this rules. 
 *  the results however are undefined, albeit sometimes usable. 
 * 
 * in contrast to other slicers, this one does not need:
 * * triangles don't have to carry valid normals
 * * triangle vertices don't need a defined clockwise or counterclockwise order
 * 
 * Paul Geisler 2012 
 */

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


// data structures and operations

// a 3d vertex
struct Vertex{
	float x,y,z;	
};
bool operator==(const Vertex& a, const Vertex& b){
	return a.x==b.x && a.y==b.y && a.z==b.z;
}
bool operator!=(Vertex& a, Vertex& b){
	return !(a==b);
}
// dot product
float dot(const Vertex& a, const Vertex& b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}
// distance of two vertices
float distance(Vertex& a, Vertex& b){
	float dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
// difference of two vertices
Vertex operator-(const Vertex& a, const Vertex& b)
{
	Vertex r={a.x-b.x,a.y-b.y,a.z-b.z};
	return r;
}
// simple half ordering in z,y,x order
// automatically used by orderes containers like std::set, std::map
bool operator<(const Vertex& a, const Vertex& b)
{
	return 
		a.z< b.z             ||
		a.z==b.z && a.y< b.y ||
		a.z==b.z && a.y==b.y && a.x<b.x;
}
// directed order for sorting vertices along a given 3d direction
struct VertexSweepOrder
{
	Vertex dir;
	VertexSweepOrder(Vertex _dir){
		dir=_dir;
	}
	bool operator() (const Vertex& a, const Vertex& b) const
	{
		return dot(a,dir)<dot(b,dir); 
	}
};

// a triangle referencing three vertices
struct Triangle
{
	Vertex* vertices[3];
};
void sortTriangleVertices(Triangle& t)
{
	Vertex** vs=t.vertices;
	if(vs[0]->z > vs[1]->z) std::swap(vs[0],vs[1]);
	if(vs[0]->z > vs[2]->z) std::swap(vs[0],vs[2]);
	if(vs[1]->z > vs[2]->z) std::swap(vs[1],vs[2]);
}

// helper struct for sorting vertices by some value
struct VertexIndex
{
	float value; 
	Triangle* triangle;
};
bool operator<(const VertexIndex& a, const VertexIndex& b)
{
	return a.value<b.value;	
};

// a line segment 
// usually resulting of the intersection from a triangle with a z plane
struct Segment
{
	std::array<Vertex,2> vertices;     // the two endpoints of this segment
	// these are used temporary and are invalidated later. do not use them:
	std::array<Segment*,2> neighbours; // pointers to this segment adjacent ones
	long orderIndex;  // an index generated to order segments for efficient printing
};
bool operator<(const Segment& a, const Segment& b){
	return a.orderIndex<b.orderIndex;
}

// a layer holding the segments build by intersecting the mesh with a z plane
struct Layer
{
	float z; // z plane 
	std::vector<Triangle*> triangles;  // triangles touching this layer
	std::vector<Segment> segments;     // segments generated for printing
};

// a map holding configuration values read from config.ini
std::map<std::string, float> config;

// read a config value
float get(const char* parameter){
	// ensure sure the value was set by the config file
	assert(config.count(parameter)==1);
	
	return config[parameter];
}


// functions

// load config file in Slic3r format
void loadConfig(const char* filename){
	printf("Loading config %s...\n",filename);
	FILE* file=fopen(filename,"r");	
	char line[256];
	while(!feof(file)){
		if(fgets(line, sizeof(line), file)){
			char name[256];
			float value;
			// scan for one 'name = value' entry
			// currenlty only handles single float values
			int found=sscanf(line,"%s = %e",name, &value);
			if(found) 
				config[name]=value; // store to config map			
		}
	}
	fclose(file);
}

// load an ASCII .stl file 
// fill the vertices and triangle list. the vertices are unified while loading.
void loadStl(const char* filename, std::vector<Vertex>& vertices, std::vector<Triangle>& triangles)
{
	// as .stl stores unconnected triangles, any vertex found is usually repeated in
	// several more triangles. to remesh that heap of triangles, we unify those vertices.
	
	// collection of unique vertices and their index number 
	std::map<Vertex,int> uniqueVertices;
	
	// collection of vertex indicies to reconstruct triangles after vertex merging
	std::vector<int> indices;

	printf("Loading %s...\n",filename);
	FILE* file=fopen(filename,"r");	
	char line[256];
	while(!feof(file)){
		// read file line by line
		if(fgets(line, sizeof(line), file)){
			Vertex p;
			// we only scan for vertex definitions, their triangle linking is given by 
			// groups of three consecutive definitions.
			int found=sscanf(line," vertex %e %e %e",&p.x,&p.y,&p.z);
			if(found) { 
				// we found a vertex declaration
				// check if this vertex is already known
				int index;
				if(uniqueVertices.count(p)==1) {
					// we know the vertex, so get its index
					index=uniqueVertices[p];					
				}else{
					// this is a new vertex, so store it 
					vertices.push_back(p);
					index=vertices.size()-1; // the new vertex is the last element
					uniqueVertices[p]=index; // store index
				}
				indices.push_back(index); // store index in triangle order
			}
		}
	}
	fclose(file);

	// if we read triangles only, there are triangles*3 vertex indices
	assert(indices.size()%3 == 0);

	// create triangles
	for(int i=0; i<indices.size(); i+=3)
	{
		Triangle t;
		// store vertex pointers 
		for(int j=0; j<3; j++)
			t.vertices[j]=&vertices[indices[i+j]];
		// sort vertices bottom-up for later operations
		sortTriangleVertices(t);
		// add triangle
		triangles.push_back(t);	
	}

	printf("Loading complete: %d vertices read, %d unique, %d triangles.\n",(int)indices.size(),(int)vertices.size(),(int)triangles.size());
}

// create initialized layers and assign triangles to them
void buildLayers(std::vector<Triangle>& triangles, std::vector<Layer>& layers)
{

	// we compute the infill by using a 'plane sweep'.
	// see http://en.wikipedia.org/wiki/Sweep_line_algorithm
	// for that we create an index of vertices sorted by z and iterate it while keeping a heap of 
	// triangles currently touching the sweep plane. As there is no topological change between the 
	// sweep locations, any layer inbetween can be initialized with triangles intersected by that layer. 

	// create a index into the vertices sorted by z
	std::vector<VertexIndex> by_z;
	for(int i=0; i<triangles.size(); i++)
		for(char j=0; j<3; j++){
			VertexIndex vi={
				triangles[i].vertices[j]->z,
				&triangles[i]
			};
			by_z.push_back(vi);
		}
	std::sort(by_z.begin(), by_z.end());

	// sweep heap: updated list of triangles touched by the current sweep plane
	// and the number of triangle vertices already passed by the sweep plane
	std::map<Triangle*,int> activeTriangles;

	// print geometric height
	float min_z=by_z.front().value;
	float max_z=by_z.back().value;
	float layer_height=get("layer_height");
	assert(layer_height>0);
	printf("Slicing from %f to %f\n",min_z, max_z);
	
	// now do the sweep over all vertices, interrupted at every next_layer_z to fill
	float next_layer_z=min_z+layer_height;
	for(int i=0; i<by_z.size(); i++)
	{
		float z=by_z[i].value;
		// add all layers passed by the sweep so far
		while(z>next_layer_z) {
			// create new layer
			Layer layer;
			layer.z=next_layer_z;
			// copy triangle pointers to the layer
			for(std::map<Triangle*,int>::iterator j=activeTriangles.begin(); j!=activeTriangles.end(); ++j)
				layer.triangles.push_back(j->first);
			// add layer to list
			layers.push_back(layer);
			// advance to next layer height
			next_layer_z+=layer_height;
		}
		// now the layers are on par with the sweep
		
		// update the heap by the current vertex
		// get triangle this vertex is of
		Triangle* triangle=by_z[i].triangle;
		int verticesPassed; // how many vertices of the triangle we passed
		if(activeTriangles.count(triangle)==0) 
			// ta new triangle is encountered, we just see its first vertex
			verticesPassed=1;  
		else
			// we already know this triangle, so we pass another of its vertices
			verticesPassed=activeTriangles[triangle]+1;
		// store new count
		activeTriangles[triangle]=verticesPassed;
		// if we reach a third vertex, it is the triangle's top vertex
		// so we passed it completely. remove it. 
		if(verticesPassed==3) 
			activeTriangles.erase(triangle);			
		// the heap now contains an up to date collection of active triangles
	}
	// the sweep is over, any triangle should have been passed and removed.
	assert(activeTriangles.size()==0);

	// print amount of layers found.
	printf("Layers: %d\n",(int)layers.size());
}

// compute intersection of a segment given by two vertices with a z plane
Vertex computeIntersection(Vertex& a, Vertex& b, float z)
{
	float t=(z-a.z)/(b.z-a.z); // projective length
	Vertex intersection={
		a.x+t*(b.x-a.x),
		a.y+t*(b.y-a.y),
		z
	};
	
	return intersection;
}

// compute intersection of a triangle with a z plane
// the triangle vertices must be ordered in z 
std::array<Vertex,2> computeSegment(Triangle& t, float z)
{
	Vertex** vs=t.vertices;
	
	// two vertices to return
	std::array<Vertex,2> segment;

	// triangle vertices are always ordered by z
	assert(vs[0]->z<=vs[1]->z);
	assert(vs[1]->z<=vs[2]->z);	

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

// compute 'infill', a hatching pattern to fill the inner area of a layer
// it is made by a line grid alternating between +/-45 degree on odd and even layers
void fill(int layerIndex, Layer& layer, std::map<Vertex,std::vector<Segment*>> segmentsByVertex)
{

	// we compute the infill by using a 'plane sweep'.
	// see http://en.wikipedia.org/wiki/Sweep_line_algorithm
	// for that the vertices are ordered in the fill pattern hatching direction 
	// the vertices are then iterated one by one and a heap of active segments is maintained
	// that can then be used to efficiently intersect the pattern lines at the given cut
	
	// place grid lines by nozzle diameter for 100% infill
	float grid_spacing=get("nozzle_diameter");
	
	// 45 degree hatching pattern directions
	Vertex dir          ={sqrt(2.f)/2,sqrt(2.f)/2,0};
	Vertex dirOrthogonal={dir.y,-dir.x,0};
	
	// every even layer we swap the directions to get a plywood like 3d pattern
	if(layerIndex%2==0) 
		std::swap(dir,dirOrthogonal);
	
	// initialize two orders that sort vertices along dir or dirOrthogonal
	VertexSweepOrder order          (dir);
	VertexSweepOrder orderOrthogonal(dirOrthogonal);
	
	// create ordered vertex list for the sweep
	std::vector<Vertex> sweepVertices;
	for(std::map<Vertex,std::vector<Segment*>>::iterator i=segmentsByVertex.begin(); i!=segmentsByVertex.end(); ++i){
		assert(i->first.z==layer.z); 
		sweepVertices.push_back(i->first);
	}
	std::sort(sweepVertices.begin(),sweepVertices.end(),VertexSweepOrder(dir));

	// the list of infill line segments		
	std::vector<Segment> infill;
	
	// the plane sweep heap.
	// in every sweep step, this is updated to contain the segments that interact with a hatching line
	std::set<Segment*> sweepHeap;
	
	// the sweep progress distance in direction dir. 
	// initialize by the first vertex in dir
	float sweepT=dot(*sweepVertices.begin(),dir)+grid_spacing;
	
	// ascending index written to the segments to sort them later
	long orderIndex=0;
	
	// the plane sweep, hopping from vertex to vertex along dir
	for(std::vector<Vertex>::iterator i=sweepVertices.begin(); i!=sweepVertices.end(); ++i)
	{
		const Vertex& v=*i;
		
		// check if the sweep has passed the next crosshatch line
		// and fill lines until the current sweep vertex is reached
		while(dot(v,dir)>sweepT){	
		
			// collect intersections	
			std::vector<Vertex> intersections;
			for(std::set<Segment*>::iterator j=sweepHeap.begin(); j!=sweepHeap.end(); ++j){
				Segment& s=**j;
				Vertex& a=s.vertices[0], &b=s.vertices[1];
				
				// compute intersection length on segment
				float aInDir=dot(a,dir);
				float bInDir=dot(b,dir);
				assert(std::min(aInDir,bInDir)<sweepT+.00001f);
				// assert(sweepT<max(aInDir,bInDir)); // disabled to accept non manifolds
				float d=bInDir-aInDir;
				float t=(sweepT-aInDir)/(bInDir-aInDir);
				
				// add intersection
				if(t>=0 && t<=1){ // for manifolds this should always be true
					Vertex intersection={
						a.x+t*(b.x-a.x),
						a.y+t*(b.y-a.y),
						a.z // z const in layer
					};
					intersections.push_back(intersection);
				}
			}
			// assert(intersections.size() % 2 == 0);  // disabled to accept non manifolds

			// sort intersections in dirOrthogonal, perpendicular to the sweep direction
			std::sort(intersections.begin(), intersections.end(),orderOrthogonal);
			
			// add fill line segments
			// the filling toggles on every intersection, starting with the leftmost outline 
			// a pathIndex is used to keep the generated segments ordered by the path taken first
			// otherwise the printer would need to fill the disconnected hatch line using useless travels
			// TODO as pathes are not identifiable, senseless travels occur where pathes appear and vanish
			long pathIndex=1;
			for(int j=0; j<intersections.size(); j++){
				if(j%2==1 && distance(intersections[j-1],intersections[j])>1) {
					Segment s={
						{intersections[j-1],intersections[j]},
						{NULL,NULL},
						pathIndex*0x10000L + orderIndex // order by path, then by cut index
					};

					// add segment
					infill.push_back(s);
					
					// advance mayor sort index for every path
					pathIndex++;
				}
			}
		
			sweepT+=grid_spacing;
			orderIndex++; // advance minor sort index
		}
		// the filling is on par, now update sweep heap 
		
		std::vector<Segment*>& ss=segmentsByVertex[v]; // the segments touching this sweep point

		// count segments linking the current vertex and already in the heap
		char segments_in_heap=0; // number of segments already in heap
		char segment_index=-1;   // store segment already there
		// ignore non manifold unconnected segments
		if(ss.size()!=2) continue;
		for(int j=0; j<ss.size(); j++){
			assert(ss[j]->vertices[0]==v || ss[j]->vertices[1]==v);
			if(sweepHeap.count(ss[j])==1) {
				segments_in_heap++;
				segment_index=j;
			}
		}

		// assert(segments_in_heap<=2); // disabled to accept non manifolds

		// now we can decide what happens at this sweep coordinate
		if(segments_in_heap==0){
			// a new perimeter is encountered. add it to the heap.
			// as perimeters are closed, there are two segments spawning from this new vertex
			sweepHeap.insert(ss[0]);
			sweepHeap.insert(ss[1]);
		}else if(segments_in_heap==1){
			// an existing perimeter continues, segment_index must be it's index.
			// so remove this old segment and add the new one
			sweepHeap.erase (ss[  segment_index]);
			sweepHeap.insert(ss[1-segment_index]);
		}else if(segments_in_heap==2){
			// an existing perimeter ends. 
			// as perimeters are closed, there are two segments ending here.
			sweepHeap.erase (ss[0]);
			sweepHeap.erase (ss[1]);
		}		
	}
	// the sweep is over, if the mesh was manifold the heap should be empty again
	// assert(sweepHeap.size()==0); // disabled to accept non manifolds
	
	layer.segments.insert(layer.segments.end(),infill.begin(),infill.end());
}

// build segments to be printed for a layer
// first, the contour gained by intersecting the triangles with it's z plane 
// second, the infill as generated by fill(..)
void buildSegments(int layerIndex, Layer& layer)
{
	// we try to build closed loops of sements for efficient printing

	// generate segments by intersecting the triangles touching this layer
	for(int i=0; i<layer.triangles.size(); i++)
	{
		Triangle* t=layer.triangles[i];
		Segment s={
			computeSegment(*t,layer.z),
			{NULL,NULL},
			-1
		};
		if(s.vertices[0]!=s.vertices[1]) 
			layer.segments.push_back(s);
	}
	
	// unify segment vertices
	std::map<Vertex,std::vector<Segment*>> segmentsByVertex;
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
		
		if(ss.size()==1) assert(!"Unconnected segment");
		// disabled to accept non manifolds
		// if(ss.size()>2 ) assert(!"Non manifold segment");
		
		Vertex v=i->first;
		
		// as we don't know the direction of each segment in the final trajectory,
		// we just link them in the same order as they list their vertices.
		// use two indices for the corresponding neighbour pointers
		// TODO maybe we should make this simpler and just use the first free neighbour pointer,
		// however errors are harder to track than.
		int index0, index1;
		if       (ss[0]->vertices[0]==v) index0=1;
		else if  (ss[0]->vertices[1]==v) index0=0;
		else     assert(!"bad index0");

		if       (ss[1]->vertices[0]==v) index1=1;
		else if  (ss[1]->vertices[1]==v) index1=0;
		else     assert(!"bad index1");
		
		// now index0, index1 should point to a free end of the segment
		assert(ss[0]->neighbours[index0]==NULL);
		assert(ss[1]->neighbours[index1]==NULL);
		
		// finally link both segments
		ss[0]->neighbours[index0]=ss[1];
		ss[1]->neighbours[index1]=ss[0];
	}

	// check for dangling segments (caused by disconnected triangles)
	// disabled to accept non manifold meshes
	/*
	for(int i=0; i<layer.segments.size(); i++)
		for(int j=0; j<2; j++)
			if(layer.segments[i].neighbours[j]==NULL) {
				printf("Unconnected segment: %d %d\n",i,j);
				throw 0;
			}
	*/
		
	// now order the segments into consecutive loops. 
	int loops=0;
	long orderIndex=0;
	for(int i=0; i<layer.segments.size(); i++){
		Segment& segment=layer.segments[i];

		// only handle new loops
		if(segment.orderIndex!=-1) continue;
				
		// collect a loop
		Segment* s2=&segment;
		while(true){
			s2->orderIndex=orderIndex++;
			// DIRTY: check for NULL neighbours to survive non manifolds
			if     (s2->neighbours[0] != NULL && s2->neighbours[0]->orderIndex==-1)
				s2=s2->neighbours[0];
			else if(s2->neighbours[1] != NULL && s2->neighbours[1]->orderIndex==-1)
				s2=s2->neighbours[1];
			else break;
		};
		
		// the loop should be closed:
		// DIRTY: ignore check to accept non manifolds		
		// assert(s2->neighbours[0]==&segment || s2->neighbours[1]==&segment);
		
		loops++;
	}

	// debug output	
	// printf("\tTriangles: %d, segments: %d, vertices: %d, loops: %d\n",(int)layer.triangles.size(),(int)layer.segments.size(),(int)segmentsByVertex.size(),loops);

	fill(layerIndex, layer,segmentsByVertex);
	std::sort(layer.segments.begin(), layer.segments.end());
	// caution: the neighbour[..] and other segment pointers are invalid now! 
}

// save Gcode
// iterates over the previously generated layers and emit gcode for every segment
// uses some configuration values to decide when to retract the filament, how much 
// to extrude and so on.
void saveGcode(const char* filename, std::vector<Layer>& layers)
{

	printf("Saving Gcode...\n");
	FILE* file=fopen(filename,"w");	

	// segments shorter than this are ignored
	float skipDistance=.01;

	// compute extrusion factor, that is the amount of filament feed over extrusion length
	float dia=get("filament_diameter");
	float filamentArea=3.14159f*dia*dia/4;
	float extrusionVolume=get("nozzle_diameter")*get("layer_height");
	float extrusionFactor=extrusionVolume/filamentArea*get("extrusion_multiplier");
	
	// retract filament if traveling
	float retract_length=get("retract_length");
	float retract_before_travel=get("retract_before_travel");

	// statistical values shown to the user
	int travels=0, longTravels=0, extrusions=0;
	int travelsSkipped=0, extrusionsSkipped=0;
	float travelled=0, extruded=0;

	// offset of the emitted Gcode coordinates to the .stl ones
	Vertex offset={75,75,get("z_offset")};

	Vertex position={0,0,0};
	for(int i=0; i<layers.size(); i++){
		Layer& l=layers[i];
		fprintf(file, "G92 E0\n");                        // reset extrusion axis
		fprintf(file, "G1 Z%f F7800.000\n",l.z+offset.z); // move to layer's z plane
		float extrusion=retract_length; // extrusion axis position

		for(int j=0; j<l.segments.size(); j++){
			Vertex& v0=l.segments[j].vertices[0];
			Vertex& v1=l.segments[j].vertices[1];
			assert(v0.z==l.z);
			assert(v1.z==l.z);
	
			// reorder segment for shorter or zero traveling
			if(distance(v1,position)<distance(v0,position))
				std::swap(v0,v1);
			
			// check distance to decide if we need to travel
			float d=distance(v0,position);
			if(d>skipDistance){ 
				// the sements are not connected, so travel without extrusion
				if(d>retract_before_travel){
					// we travel some time, do retraction
					extrusion-=retract_length;
					fprintf(file,"G1 F1800.0 E%f\n",extrusion);
					//G92 E0
				}
				// emit G1 travel command
				fprintf(file,"G1 X%f Y%f\n",v0.x+offset.x,v0.y+offset.y);
				if(d>retract_before_travel){
					// we travelled some time, undo retraction
					extrusion+=retract_length;
					fprintf(file,"G1 F1800.0 E%f\n",extrusion);
					longTravels++;
				}
				travels++;
				travelled+=distance(v0,position);
				position=v0;
			}else   // the segments where connected or not far away
				travelsSkipped++;
			
			extrusion+=extrusionFactor*distance(v0,v1); // compute extrusion by segment length
			if(distance(v1,position)>skipDistance){
				// emit G1 extrusion command
				fprintf(file,"G1 X%f Y%f E%f\n",v1.x+offset.x,v1.y+offset.y,extrusion);
				extrusions++;
				extruded+=distance(v1,position);
				position=v1;
			}else   // the segment is to short to do extrusion
				extrusionsSkipped++;			
		}
	}	
	// print some statisitcs
	printf("Saving complete. %ld bytes written. %d travels %.0f mm, %d long travels, %d extrusions %.0f mm, %d travel skips, %d extrusion skips\n",
		ftell(file),travels, travelled, longTravels, extrusions, extruded, travelsSkipped, extrusionsSkipped);
}


// main entry

int main(int argc, const char** argv)
{
	if(argc!=3) {
		printf("Usage: %s <.stl file> <.gcode file>\n",argv[0]);
		return 1;
	}

	// load config file to fill config map
	loadConfig("config.ini");

	// load the .stl file
	std::vector<Vertex>   vertices;
	std::vector<Triangle> triangles;
	loadStl(argv[1], vertices, triangles);

	// create layers and assign touched triangles to them
	std::vector<Layer> layers;
	buildLayers(triangles, layers);

	// create printable segments for every layer
	for(int i=0; i<layers.size(); i++)
		buildSegments(i,layers[i]);

	// save filled layers in Gcode format	
	saveGcode(argv[2],layers);
	
	return 0;
}



