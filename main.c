#include "stdio.h"

struct vertex{
	float x,y,z;
};

struct triangle
{
	vertex vertexs[3];
};

int main(int argc, const char** argv)
{
	if(argc!=2) {
		printf("Usage: a.out <filename>\n");
		return 1;
	}

	const char* filename=argv[1];
	FILE* file=fopen(filename,"r");

	char line[256];

	int count_vertices=0;
	while(!feof(file)){
		if(fgets(line, sizeof(line), file)){
			vertex p;
			int found=sscanf(line," vertex %e %e %e",&p.x,&p.y,&p.z);
			if(found) count_vertices++; 
		}
	}
	printf("Vertices: %d\n",count_vertices);
	
	triangle* triangles=new triangle[count_vertexs/3];
	rewind(file);
	int i=0;
	int pi=0;
	triangle t;
	while(!feof(file)){
		if(fgets(line, sizeof(line), file)){
			vertex p;	
			int found=sscanf(line," vertex %e %e %e",&p.x,&p.y,&p.z);
			if(found) {
				t.vertexs[i%3]=p;
				i++;
				if(i%3==0) triangles[i/3]=t;
			}
		}
	}
	printf("Loading complete.\n");
	fclose(file);

	vertex** vertex_index=new vertex*[count_vertices];
	for(i=0; i<count_vertices; i++) vertex_index[i]=triangles[i/3].points[i%3];
	qsort();

	delete triangles;
	return 0;
}
