// #define _GNU_SOURCE //defined for extension for qsort_r

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h> // for bool data type


////extern
// extern scalar css_test,css_test2,css_test3;
// extern scalar areasg,areasl,arealg;
// extern face vector fss_test,fss_test2,fss_test3;


// Forward declarations for cross-referencing
struct Pointt;
struct Edge;
struct Surface;
struct Line;
struct EdgeVector;

// Definitions
typedef struct Pointt {
    int number;
    double coordinates[2];
    struct Edge** edges;
    int edge_count;
    struct Surface** surfaces;
    int surface_count;
} Pointt;

typedef struct Edge {
    int number;
    Pointt* start_pointt;
    Pointt* end_pointt;
    double length;
    struct Surface** surfaces;
    int surface_count;
    int type;
} Edge;

typedef struct Surface {
    int number;
    Edge** edges;
    int edge_count;
    Pointt** pointts;
    int pointt_count;
    double area;
    int cut_count;
    int pointts_reordered;  // Add this flag
    int type;
} Surface;

// Line structure definition
typedef struct {
    int number;
    double a;
    double b;
    double c;
    double normal[2];
} Line;

typedef struct {
    Edge** data;
    int size;
    int capacity;
} EdgeVector;

typedef struct {
    Surface** data;
    int size;
    int capacity;
} SurfaceVector;


//
typedef struct {
    Pointt* intersection_pointt;
    Edge* associated_edge;
} IntersectionInfo;

// Global variables to track all instances
Pointt** all_pointts = NULL;
int pointt_count = 0;

Edge** all_edges = NULL;
int edge_count = 0;

Surface** all_surfaces = NULL;
int surface_count = 0;

Line** all_lines = NULL;
int line_count = 0;

double centroid[2];

// Function prototypes
int get_next_id(int* existing_ids, int count);
void add_edge_to_pointt(Pointt *pointt, Edge *edge);
void add_surface_to_list(Surface ***surfaces, int *count, Surface *surface);

Pointt* Pointt_create(double x, double y);
Edge* Edge_create(Pointt* start_pointt, Pointt* end_pointt);
Surface* Surface_create(Edge** edges, int edge_count, Pointt** pointts, int pointt_count, int ordered);
double Surface_calculate_area(Surface* surface);

// Function prototypes (new)
const char* SurfaceType_to_string(int type);

void Pointt_print(Pointt* pointt);
void Edge_print(Edge* edge);
void Surface_print(Surface* surface);
void Pointt_belongs_to_print(Pointt* pointt);
void Edge_belongs_to_print(Edge* edge);

double polar_angle(Pointt* pointt, double centroid_x, double centroid_y);
// int compare_pointts(const void* a, const void* b, void* centroid);
int compare_pointts(const void* a, const void* b);
void Surface_reorder_pointts(Surface* surface);
EdgeVector* init_edge_vector(int initial_capacity);
void add_edge_to_vector(EdgeVector* vec, Edge* edge);
SurfaceVector* init_surface_vector(int initial_capacity);
void add_surface_to_vector(SurfaceVector* vec, Surface* surface);

//line
Line* Line_create(int number, double a, double b, double c);
void Line_print(Line* line);
double Line_at(Line* line, double x, double y, int flag);

//function
bool are_parallel(double a, double b, double x1, double y1, double x2, double y2);
Pointt** intersect(Line* line, Edge* edge);
void create_edges_from_pointts(Pointt** pointts, int pointts_length, EdgeVector* own_edges, EdgeVector* overlap_edges);
void cut_surface(Surface* original_surface, Line* cutting_line, SurfaceVector* new_surfaces);
void cut_surfaces(SurfaceVector* surfaces, Line* cutting_line, SurfaceVector* new_surfaces);
void SurfaceVector_print(SurfaceVector* surfaceVector);

//type get
// Function to get the center of a surface based on the coordinates of its pointts
// Function to check whether a pointt is above a line
int Pointt_above_line(double pointt_coordinates[2], Line* line);
void Surface_get_center(Surface* surface, double center_coordinates[2]);
void Surface_get_type(Surface* surface, Line* line1, Line* line2);
void Surface_get_type_3_cut(Surface* surface, Line* line1, Line* line2, Line* line3);

// Function to check if a pointt lies on a line defined by two pointts
int Pointt_on_line(Pointt* pointt, Pointt* line_start, Pointt* line_end);
// Function to check if edge1 lies on edge2
int Edge_on(Edge* edge1, Edge* edge2);
void Edge_get_type(Edge* edge, EdgeVector* square_edges, SurfaceVector* surfaces_input);


//free
void free_all_pointts(Pointt** all_pointts, int pointt_count);
void free_all_edges(Edge** all_edges, int edge_count);
void free_all_surfaces(Surface** all_surfaces, int surface_count);
void free_all_lines(Line** all_lines, int line_count);
void free_all_resources();

void line_intersection_range(coord n1, double* c_min, double* c_max); 


// Main function
//int main(){
int cut_line_test_in_basilisk(coord n1,double c1,coord n2,double c2, double* data2) {
    double a1 = n1.x;
    double b1 = n1.y;
    double a2 = n2.x;
    double b2 = n2.y;
    // Create some Pointts for testing
    // printf("pointt_count=%d\n",pointt_count);
    // printf("edge_count=%d\n",edge_count);
    // printf("surface_count=%d\n",surface_count);

all_pointts = NULL;
pointt_count = 0;

all_edges = NULL;
edge_count = 0;

all_surfaces = NULL;
surface_count = 0;

all_lines = NULL;
line_count = 0;


double css_test2_value=0.0;
double css_test_value=0.0;
double arealg_value=0.0;
double areasl_value=0.0;
double areasg_value=0.0;
double fss_test_small_x=0.0;
double fss_test_small_y=0.0;
double fss_test_big_x=0.0;
double fss_test_big_y=0.0;
double fss_test2_small_x=0.0;
double fss_test2_small_y=0.0;
double fss_test2_big_x=0.0;
double fss_test2_big_y=0.0;

    Pointt* p1 = Pointt_create(0.0, 0.0);
    // printf("finish create pointt p1\n");
    Pointt* p2 = Pointt_create(1.0, 0.0);
    // printf("finish create pointt p2\n");
    Pointt* p3 = Pointt_create(0.0, 1.0);
    // printf("finish create pointt p3\n");
    Pointt* p4 = Pointt_create(1.0, 1.0);
    // printf("finish create pointt p4\n");
    // printf("finish create pointt\n");

    // Print Pointts, Edges, and Surfaces
    // printf("List of Pointts:\n");
    // for (int i = 0; i < pointt_count; ++i) {
    //     Pointt_print(all_pointts[i]);
    // }


    // Create some Edges for testing
    Edge* e1 = Edge_create(p1, p2); //bottom
    Edge* e2 = Edge_create(p2, p4); //right
    Edge* e3 = Edge_create(p4, p3); //up
    Edge* e4 = Edge_create(p3, p1); //left
    // printf("finish create edge\n");

    // printf("\nList of Edges:\n");
    // for (int i = 0; i < edge_count; ++i) {
    //     Edge_print(all_edges[i]);
    // }


    // Create a Surface for testing
    Edge* edges[] = {e1, e2, e3, e4};
    Pointt* pointts[] = {p1, p2, p4, p3};
    // printf("surface_count=%d\n",surface_count);
    Surface* s1 = Surface_create(edges, 4, pointts, 4, 0);
    // printf("finish create surface\n");

    
    
    // printf("\nList of Surfaces:\n");
    // for (int i = 0; i < surface_count; ++i) {
    //     Surface_print(all_surfaces[i]); 
    // }

    // Line* line1 = Line_create(1, 1.0, 1.0, -1.2);
    Line* line1 = Line_create(1, a1, b1, c1);
    // Line_print(line1);



    //test for create_edges_from_pointts
    if (false){
        EdgeVector* own_edges = init_edge_vector(10);
        EdgeVector* overlap_edges = init_edge_vector(10);;
        // Create edges from pointts
        Pointt* pointts_test[3];
        pointts_test[0] = p2;
        pointts_test[1] = p4;
        pointts_test[2] = p3;
        create_edges_from_pointts(pointts_test, 3, own_edges, overlap_edges);
        // printf("\nList of edges after create_edges_from_pointts:\n");
        // for (int i = 0; i < edge_count; ++i) {
        //     Edge_print(all_edges[i]);
        // }
        free(own_edges->data);
        free(own_edges);
        free(overlap_edges->data);
        free(overlap_edges);
    }

   //test for cut_surface
   if(false){
        SurfaceVector* new_surfaces = init_surface_vector(2); // Initial capacity of 2
        cut_surface(s1, line1, new_surfaces);

        //  printf("List of new Pointts:\n");
        //     for (int i = 0; i < pointt_count; ++i) {
        //         Pointt_print(all_pointts[i]);
        //     }
        // printf("new_surfaces_size=%d\n",new_surfaces->size);

        // printf("\nList of new_surfaces:\n");
        // SurfaceVector_print(new_surfaces);
        for (int i = 0; i < new_surfaces->size; ++i) {
            // printf("Surface %d:\n", i + 1);
            Surface_print(new_surfaces->data[i]);
            // printf("surface number: %d\n",new_surfaces->data[i]->number);
            // printf("surface->edge_account: %d\n",new_surfaces->data[i]->edge_count);
            // // for(int ii=0; ii<new_surfaces->data[i]->edge_count;ii++){
            // //     printf("edge: %d",new_surfaces->data[i]->edges[ii]->number);
            // // }
            // printf("surface->pointt_account: %d\n",new_surfaces->data[i]->pointt_count);
            // for(int ii=0; ii<new_surfaces->data[i]->pointt_count;++ii){
            //     printf("pointt: %d\n",new_surfaces->data[i]->pointts[ii]->number);
            // }

         }
        
        // Pointt_belongs_to_print(all_pointts[5]);
        // Edge_belongs_to_print(all_edges[1]);

        free(new_surfaces->data);
        free(new_surfaces);
   }

   //test for cut_surfaces
   if(false){
        SurfaceVector* new_surfaces1 = init_surface_vector(2); // Initial capacity of 2
        cut_surface(s1, line1, new_surfaces1);

        //  printf("List of new Pointts:\n");
        //     for (int i = 0; i < pointt_count; ++i) {
        //         Pointt_print(all_pointts[i]);
        //     }
        // printf("new_surfaces1_size=%d\n",new_surfaces1->size);

        // printf("\nList of new_surfaces:\n");
        // // SurfaceVector_print(new_surfaces);
        // for (int i = 0; i < new_surfaces1->size; ++i) {
        //     // printf("Surface %d:\n", i + 1);
        //     Surface_print(new_surfaces1->data[i]);
        //     // printf("surface number: %d\n",new_surfaces->data[i]->number);
        //     // printf("surface->edge_account: %d\n",new_surfaces->data[i]->edge_count);
        //     // // for(int ii=0; ii<new_surfaces->data[i]->edge_count;ii++){
        //     // //     printf("edge: %d",new_surfaces->data[i]->edges[ii]->number);
        //     // // }
        //     // printf("surface->pointt_account: %d\n",new_surfaces->data[i]->pointt_count);
        //     // for(int ii=0; ii<new_surfaces->data[i]->pointt_count;++ii){
        //     //     printf("pointt: %d\n",new_surfaces->data[i]->pointts[ii]->number);
        //     // }

        //  }    
        // Pointt_belongs_to_print(all_pointts[5]);
        // Edge_belongs_to_print(all_edges[1]);


        Line* line2 = Line_create(2, 1.0, 0, -0.4);
        // Line_print(line2);
        SurfaceVector* new_surfaces2 = init_surface_vector(2); 
        cut_surfaces(new_surfaces1, line2, new_surfaces2);
        



        Line* line3 = Line_create(3, 0, 1.0, -0.4);
        // Line_print(line3);
        SurfaceVector* new_surfaces3 = init_surface_vector(2); 
        cut_surfaces(new_surfaces2, line3, new_surfaces3);


        // printf("result of the first cut\n");
        // printf("size=%d\n",new_surfaces1->size);
        // SurfaceVector_print(new_surfaces1);
        // printf("result of the second cut\n");
        // printf("size=%d\n",new_surfaces2->size);
        // SurfaceVector_print(new_surfaces2);
        // printf("result of the third cut\n");
        // printf("size=%d\n",new_surfaces3->size);
        // SurfaceVector_print(new_surfaces3);


        free(new_surfaces1->data);
        free(new_surfaces1);
        free(new_surfaces2->data);
        free(new_surfaces2);
        free(new_surfaces3->data);
        free(new_surfaces3);
   }
    
//test for type
    if(true){
        SurfaceVector* new_surfaces1 = init_surface_vector(2); // Initial capacity of 2
        cut_surface(s1, line1, new_surfaces1);

        // Line* line2 = Line_create(2, 1.0, 0, -0.4);
        Line* line2 = Line_create(2, a2, b2, c2);
        // Line_print(line2);
        SurfaceVector* new_surfaces2 = init_surface_vector(2); 
        cut_surfaces(new_surfaces1, line2, new_surfaces2);

        for(int i=0;i<new_surfaces2->size;++i){
            Surface_get_type(new_surfaces2->data[i], line1, line2);
        }

        EdgeVector* square_edges = init_edge_vector(4);
        add_edge_to_vector(square_edges,e4); //left
        add_edge_to_vector(square_edges,e2); //right
        add_edge_to_vector(square_edges,e1); //bottom
        add_edge_to_vector(square_edges,e3); //top

        for (int i=0;i<new_surfaces2->size;++i) {
            for (int ii=0;ii<new_surfaces2->data[i]->edge_count;++ii){
                Edge_get_type(new_surfaces2->data[i]->edges[ii], square_edges, new_surfaces2);
                // printf("edge %d type:%d\n",new_surfaces2->data[i]->edges[ii]->number,new_surfaces2->data[i]->edges[ii]->type);
            }
        }
       
        // printf("result of the second cut\n");
        // printf("size=%d\n",new_surfaces2->size);
        // SurfaceVector_print(new_surfaces2);


        // type = 1; // not solid and gas;type = 2; //not solid and liquid
        // type = 3; //solid and gas,type= 4; //solid and liquid
        for (int i=0;i<new_surfaces2->size;++i) {
            int type_temp = new_surfaces2->data[i]->type;
            if(type_temp == 1){ //gas fraction
                css_test2_value = new_surfaces2->data[i]->area;
            }else if(type_temp == 2){ //liquid fraction
                css_test_value = new_surfaces2->data[i]->area;
            }
        }


        // if ((s1 == 1 && s2 == 2) || (s2 == 1 && s1 == 2)) type_temp=1;
        // if ((s1 == 1 && s2 == 3) || (s2 == 1 && s1 == 3)) type_temp=2;
        // if ((s1 == 1 && s2 == 4) || (s2 == 1 && s1 == 4)) type_temp=3;
        // if ((s1 == 2 && s2 == 3) || (s2 == 2 && s1 == 3)) type_temp=4;
        // if ((s1 == 2 && s2 == 4) || (s2 == 2 && s1 == 4)) type_temp=5;
        // if ((s1 == 3 && s2 == 4) || (s2 == 3 && s1 == 4)) type_temp=6;

        //   if (s == 1) baseType = 7;
        // else if (s == 2) baseType = 11;
        // else if (s == 3) baseType = 15;
        // else if (s == 4) baseType = 19;

        // // printf("baseType=%d\n",baseType);
        // // Check which edge the surface intersects with (edge1, edge2, etc)
        // // Assuming `Edge_intersects` is a function that takes an Edge and Surface and
        // // returns whether they intersect
        // for (int i = 1; i <= 4; ++i) {
        //     if (Edge_on(edge, square_edges->data[i-1])) {
        //         type_temp = baseType + (i-1);  // Adjust the type based on the intersecting edge
        //     }
        // }

         for (int i=0;i<new_surfaces2->size;++i) {
            for (int ii=0;ii<new_surfaces2->data[i]->edge_count;++ii){
                //arealg areasg areasl
                int type_temp = new_surfaces2->data[i]->edges[ii]->type;
                double length_temp = new_surfaces2->data[i]->edges[ii]->length;
                if(type_temp==1){
                    arealg_value = length_temp;
                }
                if(type_temp==2){
                    areasg_value = length_temp;
                }
                if(type_temp==5){
                    areasl_value = length_temp;
                }

                //fss_test.x[], fss_test2.x[1], area_temp = 1; //not solid and liquid
                // gas 7(left),8(right),9(bottom),10(top)
                if(type_temp==7){
                    fss_test2_small_x = length_temp;
                }else if(type_temp==8){
                    fss_test2_big_x = length_temp;
                }else if(type_temp==9){
                    fss_test2_small_y = length_temp;
                }else if(type_temp==10){
                    fss_test2_big_y = length_temp;
                }

                // liquid 11(left),12(right),13(bottom),14(top)
                if(type_temp==11){
                    fss_test_small_x = length_temp;
                }else if(type_temp==12){
                    fss_test_big_x = length_temp;
                }else if(type_temp==13){
                    fss_test_small_y = length_temp;
                }else if(type_temp==14){
                    fss_test_big_y = length_temp;
                }


            }
        }

    data2[0]=css_test2_value;
    data2[1]=css_test_value;
    data2[2]=arealg_value;
    data2[3]=areasl_value;
    data2[4]=areasg_value;
    data2[5]=fss_test_small_x;
    data2[6]=fss_test_small_y;
    data2[7]=fss_test_big_x;
    data2[8]=fss_test_big_y;
    data2[9]=fss_test2_small_x;
    data2[10]=fss_test2_small_y;
    data2[11]=fss_test2_big_x;
    data2[12]=fss_test2_big_y;

        free(new_surfaces1->data);
        free(new_surfaces1);
        free(new_surfaces2->data);
        free(new_surfaces2);
        free(square_edges->data);
        free(square_edges);

    }
    // Do something with own_edges and overlap_edges



    // Free allocated memory
    // ... Add code for freeing all dynamically allocated memory
    // (Exercise for the reader)
    
    // Free all resources
    free_all_resources();
    all_pointts = NULL;
    pointt_count = 0;

    all_edges = NULL;
    edge_count = 0;

    all_surfaces = NULL;
    surface_count = 0;

    all_lines = NULL;
    line_count = 0;
    return 0;
}


//int main(){
// int cut_line_test_in_basilisk_3_times(double a1,double b1,double c1,double a2,double b2,double c2 ,
// double a3,double b3,double c3, double* moving_liquid,double* moving_gas ) {
int cut_line_test_in_basilisk_3_times(coord n1 ,double c1,coord n2,double c2 ,
coord n3,double c3, double* moving_liquid,double* moving_gas ) {
//line1 : solid // line2: fliud //line3: face
    double a1 = n1.x;
    double b1 = n1.y;
    double a2 = n2.x;
    double b2 = n2.y;
    double a3 = n3.x;
    double b3 = n3.y;
*moving_liquid = 0.0;
*moving_gas = 0.0;
    // Create some Pointts for testing
    // printf("pointt_count=%d\n",pointt_count);
    // printf("edge_count=%d\n",edge_count);
    // printf("surface_count=%d\n",surface_count);

all_pointts = NULL;
pointt_count = 0;

all_edges = NULL;
edge_count = 0;

all_surfaces = NULL;
surface_count = 0;

all_lines = NULL;
line_count = 0;


double css_test2_value=0.0;
double css_test_value=0.0;
double arealg_value=0.0;
double areasl_value=0.0;
double areasg_value=0.0;
double fss_test_small_x=0.0;
double fss_test_small_y=0.0;
double fss_test_big_x=0.0;
double fss_test_big_y=0.0;
double fss_test2_small_x=0.0;
double fss_test2_small_y=0.0;
double fss_test2_big_x=0.0;
double fss_test2_big_y=0.0;

    Pointt* p1 = Pointt_create(0.0, 0.0);
    // printf("finish create pointt p1\n");
    Pointt* p2 = Pointt_create(1.0, 0.0);
    // printf("finish create pointt p2\n");
    Pointt* p3 = Pointt_create(0.0, 1.0);
    // printf("finish create pointt p3\n");
    Pointt* p4 = Pointt_create(1.0, 1.0);
    // printf("finish create pointt p4\n");
    // printf("finish create pointt\n");

    // Print Pointts, Edges, and Surfaces
    // printf("List of Pointts:\n");
    // for (int i = 0; i < pointt_count; ++i) {
    //     Pointt_print(all_pointts[i]);
    // }


    // Create some Edges for testing
    Edge* e1 = Edge_create(p1, p2); //bottom
    Edge* e2 = Edge_create(p2, p4); //right
    Edge* e3 = Edge_create(p4, p3); //up
    Edge* e4 = Edge_create(p3, p1); //left
    // printf("finish create edge\n");

    // printf("\nList of Edges:\n");
    // for (int i = 0; i < edge_count; ++i) {
    //     Edge_print(all_edges[i]);
    // }


    // Create a Surface for testing
    Edge* edges[] = {e1, e2, e3, e4};
    Pointt* pointts[] = {p1, p2, p4, p3};
    // printf("surface_count=%d\n",surface_count);
    Surface* s1 = Surface_create(edges, 4, pointts, 4, 0);
    // printf("finish create surface\n");

    
    
    // printf("\nList of Surfaces:\n");
    // for (int i = 0; i < surface_count; ++i) {
    //     Surface_print(all_surfaces[i]); 
    // }

    // Line* line1 = Line_create(1, 1.0, 1.0, -1.2);
    Line* line1 = Line_create(1, a1, b1, c1);
    // Line_print(line1);



    if(true){
        SurfaceVector* new_surfaces1 = init_surface_vector(2); // Initial capacity of 2
        cut_surface(s1, line1, new_surfaces1);

        // Line* line2 = Line_create(2, 1.0, 0, -0.4);
        Line* line2 = Line_create(2, a2, b2, c2);
        // Line_print(line2);
        Line* line3 = Line_create(3, a3, b3, c3);
        // Line_print(line3);
        SurfaceVector* new_surfaces2 = init_surface_vector(2); 
        cut_surfaces(new_surfaces1, line2, new_surfaces2);
        SurfaceVector* new_surfaces3 = init_surface_vector(2);
        cut_surfaces(new_surfaces2, line3, new_surfaces3); 

        for(int i=0;i<new_surfaces3->size;++i){
            Surface_get_type_3_cut(new_surfaces3->data[i], line1, line2, line3);
        }


        // type = 1; // not solid and gas;type = 2; //not solid and liquid
        // type = 3; //solid and gas,type= 4; //solid and liquid
        
        double temp1=0.0;
        double temp2=0.0;
        for (int i=0;i<new_surfaces3->size;++i) {
            int type_temp = new_surfaces3->data[i]->type;
            if(type_temp == 1){ //gas fraction
                temp1 = new_surfaces3->data[i]->area;
            }else if(type_temp == 2){ //liquid fraction
                temp2 = new_surfaces3->data[i]->area;
            }
        }
        *moving_gas = temp1;
        *moving_liquid = temp2;

        free(new_surfaces1->data);
        free(new_surfaces1);
        free(new_surfaces2->data);
        free(new_surfaces2);
        free(new_surfaces3->data);
        free(new_surfaces3);
        // free(square_edges->data);
        // free(square_edges);

    }


    // Do something with own_edges and overlap_edges



    // Free allocated memory
    // ... Add code for freeing all dynamically allocated memory
    // (Exercise for the reader)
    
    // Free all resources
    free_all_resources();
    all_pointts = NULL;
    pointt_count = 0;

    all_edges = NULL;
    edge_count = 0;

    all_surfaces = NULL;
    surface_count = 0;

    all_lines = NULL;
    line_count = 0;
    return 0;
}

// Function to get next available ID
int get_next_id(int* existing_ids, int count) {
    int id = 1;
    while (1) {
        int found = 0;
        for (int i = 0; i < count; ++i) {
            if (existing_ids[i] == id) {
                found = 1;
                break;
            }
        }
        if (!found) return id;
        id++;
    }
}

void add_edge_to_pointt(Pointt *pointt, Edge *edge) {
    // Check if the edge is already in the list
    for (int i = 0; i < pointt->edge_count; ++i) {
        if (pointt->edges[i] == edge) {
            return;  // Edge already in list, so return without adding
        }
    }

    // Reallocate memory to add new edge
    pointt->edges = (Edge**)realloc(pointt->edges, (pointt->edge_count + 1) * sizeof(Edge*));
    if (pointt->edges == NULL) {
        fprintf(stderr, "Memory allocation failed for edges\n");
        exit(1);
    }

    // Add the new edge
    pointt->edges[pointt->edge_count] = edge;
    pointt->edge_count++;
}

void add_surface_to_list(Surface ***surfaces, int *count, Surface *surface) {
    // Check if the surface is already in the list
    for (int i = 0; i < *count; ++i) {
        if ((*surfaces)[i] == surface) {
            return;  // Surface already in list, so return without adding
        }
    }

    // Reallocate memory to add new surface
    *surfaces = (Surface**)realloc(*surfaces, (*count + 1) * sizeof(Surface*));
    if (*surfaces == NULL) {
        fprintf(stderr, "Memory allocation failed for surfaces\n");
        exit(1);
    }

    // Add the new surface
    (*surfaces)[*count] = surface;
    (*count)++;
}


// Pointt function implementations
Pointt* Pointt_create(double x, double y) {
    // Search for existing Pointt with the same coordinates
    // for (int i = 0; i < pointt_count; ++i) {
    //     printf("all_pointts[%d]->coordinates=[%g,%g]\n",i,all_pointts[i]->coordinates[0],all_pointts[i]->coordinates[1]);
    // }
    for (int i = 0; i < pointt_count; ++i) {
        // printf("all_pointts[i]->coordinates[0]=%g\n",all_pointts[i]->coordinates[0]);
        if (all_pointts[i]->coordinates[0] == x && all_pointts[i]->coordinates[1] == y) {
            return all_pointts[i];
        }
    }

    Pointt* new_pointt = (Pointt*)malloc(sizeof(Pointt));
    if (new_pointt == NULL) {
        // Handle the error, typically by terminating the program or
        // returning an error code.
        printf("new_pointt == NULL\n");
        exit(1);  // or return NULL or a specific error code
    }
    new_pointt->coordinates[0] = x;
    new_pointt->coordinates[1] = y;
    new_pointt->edges = NULL;
    new_pointt->edge_count = 0;
    new_pointt->surfaces = NULL;
    new_pointt->surface_count = 0;

    // Assign a new number to the pointt
    int existing_ids[pointt_count];
    for (int i = 0; i < pointt_count; ++i) {
        existing_ids[i] = all_pointts[i]->number;
    }
    new_pointt->number = get_next_id(existing_ids, pointt_count);

    // Add to global list of all pointts
    all_pointts = (Pointt**)realloc(all_pointts, (pointt_count + 1) * sizeof(Pointt*));
    all_pointts[pointt_count++] = new_pointt;

    return new_pointt;
}

// Edge function implementations
Edge* Edge_create(Pointt* start_pointt, Pointt* end_pointt) {
    // Search for existing Edge with the same start and end pointts
    for (int i = 0; i < edge_count; ++i) {
        if ((all_edges[i]->start_pointt == start_pointt && all_edges[i]->end_pointt == end_pointt) ||
            (all_edges[i]->start_pointt == end_pointt && all_edges[i]->end_pointt == start_pointt)) {
            return all_edges[i];
        }
    }

    Edge* new_edge = (Edge*)malloc(sizeof(Edge));
    if (new_edge == NULL) {
        // Handle the error, typically by terminating the program or
        // returning an error code.
        exit(1);  // or return NULL or a specific error code
    }
    new_edge->start_pointt = start_pointt;
    new_edge->end_pointt = end_pointt;
    new_edge->surfaces = NULL;
    new_edge->surface_count = 0;

    double dx = start_pointt->coordinates[0] - end_pointt->coordinates[0];
    double dy = start_pointt->coordinates[1] - end_pointt->coordinates[1];
    new_edge->length = sqrt(dx * dx + dy * dy);

    // Assign a new number to the edge
    int existing_ids[edge_count];
    for (int i = 0; i < edge_count; ++i) {
        existing_ids[i] = all_edges[i]->number;
    }
    new_edge->number = get_next_id(existing_ids, edge_count);

    // Add to global list of all edges
    all_edges = (Edge**)realloc(all_edges, (edge_count + 1) * sizeof(Edge*));
    all_edges[edge_count++] = new_edge;

    new_edge->type = -1;

    // Add the new edge to the pointts it connects
    add_edge_to_pointt(start_pointt, new_edge);
    add_edge_to_pointt(end_pointt, new_edge);

    return new_edge;
}

// Surface function implementations
Surface* Surface_create(Edge** edges, int edge_count, Pointt** pointts, int pointt_count, int ordered) {
    // printf("1\n");
    // Check if surface with the same set of edges already exists
    for (int i = 0; i < surface_count; ++i) {
        Surface* existing_surface = all_surfaces[i];
        if (existing_surface->edge_count != edge_count) {
            continue;  // If edge count is different, it's not the same surface
        }

        // Check if all edges match
        int matched = 1;
        for (int j = 0; j < edge_count; ++j) {
            int edge_found = 0;
            for (int k = 0; k < edge_count; ++k) {
                if (edges[j]->number == existing_surface->edges[k]->number) {
                    edge_found = 1;
                    break;
                }
            }

            if (!edge_found) {
                matched = 0;
                break;
            }
        }

        if (matched) {
            return existing_surface;  // Return the existing surface if found
        }
    }

    // If you reach here, it means a surface with the same set of edges doesn't exist.
    // Hence, create a new surface
// printf("2\n");
    Surface* new_surface = (Surface*)malloc(sizeof(Surface));
    if (new_surface == NULL) {
        fprintf(stderr, "Memory allocation failed for Surface\n");
        exit(1);
    }
// printf("3\n");
    new_surface->pointt_count = pointt_count;
    new_surface->pointts = (Pointt**)malloc(pointt_count * sizeof(Pointt*));
    if (new_surface->pointts == NULL) {
        fprintf(stderr, "Memory allocation failed for pointts in Surface\n");
        exit(1);
    }
    for (int i = 0; i < pointt_count; ++i) {
        new_surface->pointts[i] = pointts[i];  // Assuming deepCopyPointt exists
    }

// printf("4\n");
    new_surface->edge_count = edge_count;
    new_surface->edges = (Edge**)malloc(edge_count * sizeof(Edge*));
    if (new_surface->edges == NULL) {
        fprintf(stderr, "Memory allocation failed for edges in Surface\n");
        exit(1);
    }
    for (int i = 0; i < edge_count; ++i) {
        new_surface->edges[i] = edges[i];  // Assuming deepCopyEdge exists
    }
//   printf("5\n");  
    

    new_surface->pointts_reordered = ordered;  // Initialize to false (0)
    // Reorder pointts only if they haven't been reordered before
    if (!new_surface->pointts_reordered) {
        Surface_reorder_pointts(new_surface);
        new_surface->pointts_reordered = 1;  // Set flag to true (1)
        // fprintf(stderr,"In principle, only the first surface create needs to be reordered,so this note could only be appear once!!");
    }



// printf("6\n");
    new_surface->area = Surface_calculate_area(new_surface);
    new_surface->cut_count = 0;
    new_surface->type = -1;
    // Get the next available surface ID
    int existing_ids[surface_count];
    for (int i = 0; i < surface_count; ++i) {
        existing_ids[i] = all_surfaces[i]->number;
    }
    new_surface->number = get_next_id(existing_ids, surface_count);
// printf("7\n");
    // Add to the global list of all surfaces
    all_surfaces = (Surface**)realloc(all_surfaces, (surface_count + 1) * sizeof(Surface*));
    if (all_surfaces == NULL) {
        fprintf(stderr, "Memory reallocation failed for all_surfaces\n");
        exit(1);
    }
    all_surfaces[surface_count++] = new_surface;
// printf("8\n");
    // Add the new surface to each edge's list of surfaces
    for (int i = 0; i < edge_count; ++i) {
        add_surface_to_list(&(edges[i]->surfaces), &(edges[i]->surface_count), new_surface);
    }

    // Add the new surface to each pointt's list of surfaces
    for (int i = 0; i < pointt_count; ++i) {
        add_surface_to_list(&(pointts[i]->surfaces), &(pointts[i]->surface_count), new_surface);
    }
// printf("9\n");

    return new_surface;
}

double Surface_calculate_area(Surface* surface) {
    // Step 1: Compute centroid of the polygon
    double centroid_x = 0.0;
    double centroid_y = 0.0;

    int pointt_count = surface->pointt_count;
    for (int i = 0; i < pointt_count; ++i) {
        centroid_x += surface->pointts[i]->coordinates[0];
        centroid_y += surface->pointts[i]->coordinates[1];
    }
    centroid_x /= pointt_count;
    centroid_y /= pointt_count;

    // Step 2: For each edge, compute the triangle formed by the edge and the centroid
    double total_area = 0.0;
    for (int i = 0; i < pointt_count; ++i) {
        double x1 = surface->pointts[i]->coordinates[0];
        double y1 = surface->pointts[i]->coordinates[1];

        double x2 = surface->pointts[(i + 1) % pointt_count]->coordinates[0];
        double y2 = surface->pointts[(i + 1) % pointt_count]->coordinates[1];

        // Step 3: Compute area of the triangle using Shoelace formula
        total_area += 0.5 * fabs(centroid_x * (y1 - y2) + x1 * (y2 - centroid_y) + x2 * (centroid_y - y1));
    }

    // Step 4: Update the surface's area
    // surface->area = total_area;
    return total_area;
}

// ... Previous definitions and function prototypes



// ... Previous code

// Function implementations (new)

// Function to check if a pointt lies on a line defined by two pointts
int Pointt_on_line(Pointt* pointt, Pointt* line_start, Pointt* line_end) {
    double dx = line_end->coordinates[0] - line_start->coordinates[0];
    double dy = line_end->coordinates[1] - line_start->coordinates[1];
    // printf("pointt: %g %g\n",pointt->coordinates[0],pointt->coordinates[1]);
    // printf("start: %g %g\n",line_start->coordinates[0],line_start->coordinates[1]);
    // printf("end: %g %g\n",line_end->coordinates[0],line_end->coordinates[1]);
    // If the line is vertical
    if (fabs(dx)<1e-20) {
        int result= (pointt->coordinates[0] == line_start->coordinates[0] &&
                ((pointt->coordinates[1]-line_start->coordinates[1])*
                (pointt->coordinates[1]-line_end->coordinates[1])<=0) );
        // printf("dx==0, result=%d\n",result);
        return result;
    }
    
    // If the line is horizontal
    if (fabs(dy)<1e-20) {
        int result= (pointt->coordinates[1] == line_start->coordinates[1] &&
                ((pointt->coordinates[0]-line_start->coordinates[0])*
                (pointt->coordinates[0]-line_end->coordinates[0])<=0) );
        // printf("dy==0, result=%d\n",result);
        return result;
    }
    
    double t = (pointt->coordinates[0] - line_start->coordinates[0]) / dx;
    
    // Check if the pointt is within the line segment
    if (t >= 0 && t <= 1) {
        // Check if the pointt lies on the line
        int result = (pointt->coordinates[1] == line_start->coordinates[1] + t * dy);
        // printf("0<=t<=1, result=%d\n",result);
        return result;
    }
    
    // printf("000\n");
    return 0;
}

// Function to check if edge1 lies on edge2
int Edge_on(Edge* edge1, Edge* edge2) {
    // Check if both the start and end pointts of edge1 lie on edge2
    return Pointt_on_line(edge1->start_pointt, edge2->start_pointt, edge2->end_pointt) &&
           Pointt_on_line(edge1->end_pointt, edge2->start_pointt, edge2->end_pointt);
}

void Edge_get_type(Edge* edge, EdgeVector* square_edges, SurfaceVector* surfaces_input) {
    // Check the number of surfaces connected to the edge
    int type_temp = -1;
    int ii_lim=edge->surface_count;
    for(int ii=0;ii<=ii_lim;ii++){
       if(ii==0){
            // printf("Edge %d belongs to surface: %d ", edge->number, edge->surfaces[ii]->number);
       }else{
            if(ii==ii_lim){
                // printf("\n");
            }else if(ii==ii_lim-1){
                // printf(" %d", edge->surfaces[ii]->number);
            }
       }
    }
    int number1=0,number2=1;
    int valid_count=edge->surface_count;
    if(1==1){
        valid_count=0;
        for(int jj=0;jj<edge->surface_count;jj++){
            int temp_jj=edge->surfaces[jj]->number;
            for(int ii=0;ii<surfaces_input->size;ii++){
                int temp_ii=surfaces_input->data[ii]->number;
                if(temp_jj==temp_ii){
                    valid_count++;
                    if(valid_count==1){
                        number1 = jj;
                    }else if(valid_count==2){
                        number2 = jj;
                    }else{
                        printf("Error edge belongs to more than 3 surfaces!!!!!\n");
                    }
                    
                }
            }
        }
    }
    // printf("valid_count=%d\n",valid_count);
    // if (edge->surface_count == 2) {
    if (valid_count == 2) {
        // int s1 = edge->surfaces[0]->type;
        // int s2 = edge->surfaces[1]->type;
        int s1 = edge->surfaces[number1]->type;
        int s2 = edge->surfaces[number2]->type;

        // Check pairs of surfaces and assign the type
        if ((s1 == 1 && s2 == 2) || (s2 == 1 && s1 == 2)) type_temp=1;
        if ((s1 == 1 && s2 == 3) || (s2 == 1 && s1 == 3)) type_temp=2;
        if ((s1 == 1 && s2 == 4) || (s2 == 1 && s1 == 4)) type_temp=3;
        if ((s1 == 2 && s2 == 3) || (s2 == 2 && s1 == 3)) type_temp=4;
        if ((s1 == 2 && s2 == 4) || (s2 == 2 && s1 == 4)) type_temp=5;
        if ((s1 == 3 && s2 == 4) || (s2 == 3 && s1 == 4)) type_temp=6;
    }else if (valid_count == 1) {   
    // else if (edge->surface_count == 1) {
         
        // int s = edge->surfaces[0]->type;
        int s = edge->surfaces[number1]->type;
        // printf("s=%d\n",s);
        int baseType;
        
        if (s == 1) baseType = 7;
        else if (s == 2) baseType = 11;
        else if (s == 3) baseType = 15;
        else if (s == 4) baseType = 19;

        // printf("baseType=%d\n",baseType);
        // Check which edge the surface intersects with (edge1, edge2, etc)
        // Assuming `Edge_intersects` is a function that takes an Edge and Surface and
        // returns whether they intersect
        for (int i = 1; i <= 4; ++i) {
            if (Edge_on(edge, square_edges->data[i-1])) {
                type_temp = baseType + (i-1);  // Adjust the type based on the intersecting edge
            }
        }
    }

    edge->type = type_temp;
}



void Pointt_print(Pointt* pointt) {
    printf("Pointt %d: (%.2f, %.2f)\n", pointt->number, pointt->coordinates[0], pointt->coordinates[1]);
}

void Edge_print(Edge* edge) {
    printf("Edge %d: Pointts %d to %d\n", edge->number, edge->start_pointt->number, edge->end_pointt->number);
}

const char* SurfaceType_to_string(int type) {
    switch (type) {
        case 1: return "Not solid and gas";
        case 2: return "Not solid and liquid";
        case 3: return "Solid and gas";
        case 4: return "Solid and liquid";
        default: return "Unknown type";
    }
}

void Surface_print(Surface* surface) {
    char type_s[20];
    printf("Surface %d: Type: %d(%s), Area: %.2f, Edges [",  surface->number, surface->type, SurfaceType_to_string(surface->type), surface->area);
    for (int i = 0; i < surface->edge_count; ++i) {
        printf("%d", surface->edges[i]->number);
        if (i < surface->edge_count - 1) {
            printf(", ");
        }
    }

    printf("], Pointts [");
    for (int i = 0; i < surface->pointt_count; ++i) {
        printf("%d", surface->pointts[i]->number);
        if (i < surface->pointt_count - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}

void Pointt_belongs_to_print(Pointt* pointt) {
    printf("Pointt %d: belongs to Edge: [", pointt->number);
    for (int i=0;i<pointt->edge_count;i++) {
        if(i==0){
            printf("%d",pointt->edges[i]->number);
        }else{
            printf(", %d ",pointt->edges[i]->number);
        }
        
    }
    printf("]\n");
    printf("Pointt %d: belongs to surface: [", pointt->number);
    for (int i=0;i<pointt->surface_count;i++) {
        if(i==0){
            printf("%d",pointt->surfaces[i]->number);
        }else{
            printf(",%d",pointt->surfaces[i]->number);
        }
    }
    printf("]\n");
}

void Edge_belongs_to_print(Edge* edge) {
    printf("Edge %d: belongs to surface: [", edge->number);
    for (int i=0;i<edge->surface_count;i++) {
        if(i==0){
            printf("%d",edge->surfaces[i]->number);
        }else{
            printf(",%d",edge->surfaces[i]->number);
        }
    }
    printf("]\n");
}

// Function to compute polar angle with respect to a centroid
//counter-clockwise
double polar_angle(Pointt* pointt, double centroid_x, double centroid_y) {
    return atan2(pointt->coordinates[1] - centroid_y, pointt->coordinates[0] - centroid_x);
}

// Comparison function to be used in qsort
// int compare_pointts(const void* a, const void* b, void* centroid) {
//to avoid using qsort_r, set centroid to gloabal variable
int compare_pointts(const void* a, const void* b) {
    Pointt* pointt_a = *(Pointt**)a;
    Pointt* pointt_b = *(Pointt**)b;
    double centroid_x = ((double*)centroid)[0];
    double centroid_y = ((double*)centroid)[1];

    double angle_a = polar_angle(pointt_a, centroid_x, centroid_y);
    double angle_b = polar_angle(pointt_b, centroid_x, centroid_y);

    if (angle_a < angle_b) return -1;
    if (angle_a > angle_b) return 1;
    return 0;
}

// Function to reorder the pointts in a Surface object in a clockwise sequence
void Surface_reorder_pointts(Surface* surface) {
    // Step 1: Compute centroid of the pointts
    double centroid_x = 0.0, centroid_y = 0.0;
    for (int i = 0; i < surface->pointt_count; ++i) {
        centroid_x += surface->pointts[i]->coordinates[0];
        centroid_y += surface->pointts[i]->coordinates[1];
    }
    centroid_x /= surface->pointt_count;
    centroid_y /= surface->pointt_count;

    // Step 2: Sort pointts based on polar angle with respect to centroid
    // double centroid[2] = { centroid_x, centroid_y };
    centroid[0] = centroid_x;
    centroid[1] = centroid_y;
    // qsort_r(surface->pointts, surface->pointt_count, sizeof(Pointt*), compare_pointts, centroid);
    qsort(surface->pointts, surface->pointt_count, sizeof(Pointt*), compare_pointts);
}

//EdgeVector
EdgeVector* init_edge_vector(int initial_capacity) {
    EdgeVector* vec = malloc(sizeof(EdgeVector));
    if (vec == NULL) {
        // Handle memory allocation failure
        exit(1);
    }
    vec->data = malloc(initial_capacity * sizeof(Edge*));
    if (vec->data == NULL) {
        // Handle memory allocation failure
        exit(1);
    }
    vec->size = 0;
    vec->capacity = initial_capacity;
    return vec;
}

void add_edge_to_vector(EdgeVector* vec, Edge* edge) {
    if (vec->size == vec->capacity) {
        vec->capacity *= 2;  // Double the capacity
        vec->data = realloc(vec->data, vec->capacity * sizeof(Edge*));
        if (vec->data == NULL) {
            // Handle memory allocation failure
            exit(1);
        }
    }
    vec->data[vec->size] = edge;
    vec->size++;
}

//SurfaceVector
SurfaceVector* init_surface_vector(int initial_capacity) {
    SurfaceVector* vec = malloc(sizeof(SurfaceVector));
    if (vec == NULL) {
        // Handle memory allocation failure
        exit(1);
    }
    vec->data = malloc(initial_capacity * sizeof(Surface*));
    if (vec->data == NULL) {
        // Handle memory allocation failure
        exit(1);
    }
    vec->size = 0;
    vec->capacity = initial_capacity;
    return vec;
}

void add_surface_to_vector(SurfaceVector* vec, Surface* surface) {
    if (vec->size == vec->capacity) {
        vec->capacity *= 2;  // Double the capacity
        vec->data = realloc(vec->data, vec->capacity * sizeof(Surface*));
        if (vec->data == NULL) {
            // Handle memory allocation failure
            exit(1);
        }
    }
    vec->data[vec->size] = surface;
    vec->size++;
}


// Line constructor function
Line* Line_create(int number, double a, double b, double c) {
    Line* new_line = (Line*) malloc(sizeof(Line));
    if (new_line == NULL) {
        // Handle the error, typically by terminating the program
        exit(1);
    }
    new_line->number = number;
    new_line->a = a;
    new_line->b = b;
    new_line->c = c;
    new_line->normal[0] = a;
    new_line->normal[1] = b;

    all_lines = (Line**)realloc(all_lines, (line_count + 1) * sizeof(Line*));
    if (all_lines == NULL) {
        // Handle the error, typically by terminating the program or
        // returning an error code.
        exit(1);
    }
    all_lines[line_count++] = new_line;


    return new_line;
}

// Function to print a Line
void Line_print(Line* line) {
    printf("Line %d: %fx + %fy + %f = 0\n", line->number, line->a, line->b, line->c);
}

// Function to find either x or y on the line given the other coordinate
double Line_at(Line* line, double x, double y, int flag) {
    if (flag == 0 && line->b != 0) {
        // Compute y based on x
        return (-line->a * x - line->c) / line->b;
    } else if (flag == 1 && line->a != 0) {
        // Compute x based on y
        return (-line->b * y - line->c) / line->a;
    } else {
        // This is an error condition, you can decide how you want to handle it
        printf("Invalid arguments or line parameters\n");
        exit(1);
    }
}


bool are_parallel(double a, double b, double x1, double y1, double x2, double y2) {
    // Direction vector of the line through the pointts
    double dx = x2 - x1;
    double dy = y2 - y1;

    // If both vectors have x-coordinate as 0, they are parallel
    if (a == 0.0 && dx == 0.0) {
        return true;
    }
    
    // If only one of the vectors has x-coordinate as 0, they are not parallel
    if ((a == 0.0 && dx != 0.0) || (a != 0.0 && dx == 0.0)) {
        return false;
    }
    
    // Calculate the ratio of the x and y components of the two vectors
    double ratio_a = a / dx;
    double ratio_b = b / dy;

    // If the ratios are equal, the vectors are parallel
    return fabs(ratio_a - ratio_b) < 1e-9;
}


Pointt** intersect(Line* line, Edge* edge) {
    // printf("yeye: %p\n", edge);
    // printf("yeye2: %p %p\n", edge->start_pointt, edge->end_pointt);

    double a = line->a;
    double b = line->b;
    double c = line->c;

    double x1 = edge->start_pointt->coordinates[0];
    double y1 = edge->start_pointt->coordinates[1];
    double x2 = edge->end_pointt->coordinates[0];
    double y2 = edge->end_pointt->coordinates[1];

    double det = a * (x2 - x1) + b * (y2 - y1);
    // printf("yeye3: %f %f %f %f %f\n", x1, y1, x2, y2, det);

    if (fabs(det) < 1e-9) {
            Pointt** pointts = (Pointt**)malloc(2 * sizeof(Pointt*));
            if (fabs(a * x1 + b * y1 + c) < 1e-9) {
                pointts[0] = edge->start_pointt;
                pointts[1] = edge->end_pointt;
            } else {
                pointts[0] = NULL;
                pointts[1] = NULL;
            }
            return pointts;
        }

        double t = -(a * x1 + b * y1 + c) / det;
        // printf("t= %f\n", t);

        if (0 <= t && t <= 1) {
            double x = x1 + t * (x2 - x1);
            double y = y1 + t * (y2 - y1);
            // printf("need to create\n%f %f\n", x, y);
            Pointt** pointts = (Pointt**)malloc(2 * sizeof(Pointt*));
            pointts[0] = Pointt_create(x, y);  // Automatically assigns a new number
            pointts[1] = NULL;
            return pointts;
        }
        return NULL;
}


// void create_edges_from_pointts(Pointt** pointts, int pointts_length, Edge*** Existing_edges, Edge*** own_edges) {
void create_edges_from_pointts(Pointt** pointts, int pointts_length, EdgeVector* own_edges, EdgeVector* overlap_edges){
    int initial_total_edges = edge_count;
    
    for (int i = 0; i < pointts_length; ++i) {
        int next_i = (i + 1) % pointts_length;  // To wrap around to the first pointt

        Pointt* start_pointt = pointts[i];
        Pointt* end_pointt = pointts[next_i];

        // Assuming Edge_create either adds a new edge or returns an existing one
        Edge* edge = Edge_create(start_pointt, end_pointt);


        if (edge_count > initial_total_edges) {
            // The total number of edges has increased, which means this is a new edge
            add_edge_to_vector(own_edges, edge);
            initial_total_edges=edge_count;
        } else {
            // This is an existing edge
            add_edge_to_vector(overlap_edges, edge);
        }
    }
}


//cut surface

void cut_surface(Surface* original_surface, Line* cutting_line, SurfaceVector* new_surfaces) {
    // Initialize an array for intersections
    IntersectionInfo* intersections = NULL;
    int intersection_count = 0;

    // Go through each edge of the surface
    for (int i = 0; i < original_surface->edge_count; ++i) {
        Edge* edge = original_surface->edges[i];
        
        // Use your intersect function to find intersections
        Pointt** result = intersect(cutting_line, edge);

        if (result) {
            for (int j = 0; result[j] != NULL; ++j) {
                bool already_exists = false;

                // Check if the pointt is already in the intersections array
                for (int k = 0; k < intersection_count; ++k) {
                    if (intersections[k].intersection_pointt == result[j]) {
                        already_exists = true;
                        break;
                    }
                }

                // If the pointt does not already exist in the array, add it
                if (!already_exists) {
                    intersections = realloc(intersections, (intersection_count + 1) * sizeof(IntersectionInfo));
                    intersections[intersection_count].intersection_pointt = result[j];
                    intersections[intersection_count].associated_edge = edge;
                    intersection_count++;
                }
            }

            free(result);  // Free the memory allocated for 'result' if necessary
        }
    }

    // Create a SurfaceVector to hold the results
    // SurfaceVector* resulting_surfaces = init_surface_vector(2);  // Initialize with capacity 2, can be changed
    // printf("intersection_count=%d\n",intersection_count);

    if (intersection_count != 2) {
        // If there are not exactly 2 intersection pointts, add the original surface to the result
        add_surface_to_vector(new_surfaces, original_surface);
    }
    else {
            // Make a copy of the original pointts
        Pointt** original_pointts = malloc(original_surface->pointt_count * sizeof(Pointt*));
        for (int i = 0; i < original_surface->pointt_count; ++i) {
            original_pointts[i] = original_surface->pointts[i];
        }
        int original_pointt_count = original_surface->pointt_count;

        // Insert intersection pointts
        for (int i = 0; i < intersection_count; ++i) {
            Pointt* intersection_pointt = intersections[i].intersection_pointt;
            Edge* intersected_edge = intersections[i].associated_edge;
            
            int insert_index1 = -1;
            int insert_index2 = -1;
            
            for (int j = 0; j < original_pointt_count; ++j) {
                if (original_pointts[j] == intersected_edge->start_pointt) {
                    insert_index1 = j;
                }
                if (original_pointts[j] == intersected_edge->end_pointt) {
                    insert_index2 = j;
                }
            }

            int position = -1;
            if (abs(insert_index1 - insert_index2) == 1) {
                position = (insert_index1 > insert_index2) ? insert_index1 : insert_index2;
            } else if ((insert_index1 == 0 && insert_index2 == original_pointt_count - 1) ||
                    (insert_index2 == 0 && insert_index1 == original_pointt_count - 1)) {
                position = original_pointt_count;
            }

            if (position != -1) {
                original_pointt_count++;
                original_pointts = realloc(original_pointts, original_pointt_count * sizeof(Pointt*));
                for (int j = original_pointt_count - 1; j > position; --j) {
                    original_pointts[j] = original_pointts[j - 1];
                }
                original_pointts[position] = intersection_pointt;
            }
        }

        // Find indices of intersection pointts
        int index1 = -1;
        int index2 = -1;
        for (int i = 0; i < original_pointt_count; ++i) {
            if (original_pointts[i] == intersections[0].intersection_pointt) {
                index1 = i;
            }
            if (original_pointts[i] == intersections[1].intersection_pointt) {
                index2 = i;
            }
        }

        // Create new surface pointts assuming clockwise ordering
        Pointt** new_surface1_pointts = NULL;
        Pointt** new_surface2_pointts = NULL;
        int count1 = 0, count2 = 0;

        if (index1 < index2) {
            // Calculate the lengths
            count1 = index2 - index1 + 1;
            count2 = original_pointt_count - index2 + index1 + 1;

            // Allocate memory
            new_surface1_pointts = malloc(count1 * sizeof(Pointt*));
            new_surface2_pointts = malloc(count2 * sizeof(Pointt*));

            // Populate new_surface1_pointts
            for (int i = 0; i < count1; ++i) {
                new_surface1_pointts[i] = original_pointts[index1 + i];
            }

            // Populate new_surface2_pointts
            for (int i = 0; i < original_pointt_count - index2; ++i) {
                new_surface2_pointts[i] = original_pointts[index2 + i];
            }
            for (int i = 0; i < index1 + 1; ++i) {
                new_surface2_pointts[original_pointt_count - index2 + i] = original_pointts[i];
            }
        } else {
            // Calculate the lengths
            count1 = original_pointt_count - index1 + index2 + 1;
            count2 = index1 - index2 + 1;

            // Allocate memory
            new_surface1_pointts = malloc(count1 * sizeof(Pointt*));
            new_surface2_pointts = malloc(count2 * sizeof(Pointt*));

            // Populate new_surface1_pointts
            for (int i = 0; i < original_pointt_count - index1; ++i) {
                new_surface1_pointts[i] = original_pointts[index1 + i];
            }
            for (int i = 0; i < index2 + 1; ++i) {
                new_surface1_pointts[original_pointt_count - index1 + i] = original_pointts[i];
            }

            // Populate new_surface2_pointts
            for (int i = 0; i < count2; ++i) {
                new_surface2_pointts[i] = original_pointts[index2 + i];
            }
        }

        // for (int i=0; i<count2; ++i) {
        //     printf("surface 2 contains: %d\n",new_surface2_pointts[i]->number);
        // }
//
            EdgeVector* new_surface1_edges = init_edge_vector(10);  // Initialize with capacity 10
            EdgeVector* overlaps1 = init_edge_vector(10);  // Initialize with capacity 10

            // Assuming new_surface1_pointts is an array of Pointt* with count new_surface1_pointt_count
            create_edges_from_pointts(new_surface1_pointts, count1, new_surface1_edges, overlaps1);


            EdgeVector* new_surface2_edges = init_edge_vector(10);  // Initialize with capacity 10
            EdgeVector* overlaps2 = init_edge_vector(10);  // Initialize with capacity 10

            // Assuming new_surface2_pointts is an array of Pointt* with count new_surface2_pointt_count
            create_edges_from_pointts(new_surface2_pointts, count2, new_surface2_edges, overlaps2);

            // At this pointt, you have new_surface1_edges, overlaps1, new_surface2_edges, and overlaps2 filled
            // You would now go ahead to create your new Surface instances

            // Assuming Surface_create is a function you've defined to create a new Surface
            // Edge** is replaced by EdgeVector for dynamic sizing, adjust accordingly

            int totalEdgeCount1 = new_surface1_edges->size + overlaps1->size;
            Edge** allEdges1 = malloc(totalEdgeCount1 * sizeof(Edge*));

            // Copy over the edges from new_surface1_edges
            for (int i = 0; i < new_surface1_edges->size; ++i) {
                allEdges1[i] = new_surface1_edges->data[i];
            }

            // Copy over the edges from overlaps1
            for (int i = 0; i < overlaps1->size; ++i) {
                allEdges1[new_surface1_edges->size + i] = overlaps1->data[i];
            }

            // printf("totalEdgeCount1: %d \n", totalEdgeCount1);
            // for (int i=0;i<totalEdgeCount1;++i) {
            //     printf("allEdges1: %d \n", allEdges1[i]->number);
            // }
            // printf("count1: %d \n", count1);
            // for (int i=0; i<count1; ++i) {
            //     printf("new_surface1_pointts: %d\n",new_surface1_pointts[i]->number);
            // }
            // Now create the new surface
            Surface* new_surface1 = Surface_create(allEdges1, totalEdgeCount1, new_surface1_pointts, count1, 1);

            // for (int i=0;i<new_surface1->edge_count;++i) {
            //     printf("new_surface1_edges: %d \n", new_surface1->edges[i]->number);
            // }
            // printf("new_surface1->pointt_count: %d \n", new_surface1->pointt_count);
            // for (int i=0; i<new_surface1->pointt_count; ++i) {
            //     printf("new_surface1: %d\n",new_surface1->pointts[i]->number);
            // }
            // Don't forget to free allEdges1 when you're done with it

            free(allEdges1);


            int totalEdgeCount2 = new_surface2_edges->size + overlaps2->size;
            Edge** allEdges2 = malloc(totalEdgeCount2 * sizeof(Edge*));

            // Copy over the edges from new_surface1_edges
            for (int i = 0; i < new_surface2_edges->size; ++i) {
                allEdges2[i] = new_surface2_edges->data[i];
            }

            // Copy over the edges from overlaps1
            for (int i = 0; i < overlaps2->size; ++i) {
                allEdges2[new_surface2_edges->size + i] = overlaps2->data[i];
            }

            // Now create the new surface
            Surface* new_surface2 = Surface_create(allEdges2, totalEdgeCount2, new_surface2_pointts, count2, 1);


            //  for (int i=0;i<new_surface2->edge_count;++i) {
            //     printf("new_surface2_edges: %d \n", new_surface2->edges[i]->number);
            // }
            // printf("new_surface2->pointt_count: %d \n", new_surface2->pointt_count);
            // for (int i=0; i<new_surface2->pointt_count; ++i) {
            //     printf("new_surface2: %d\n",new_surface2->pointts[i]->number);
            // }

            // Don't forget to free allEdges1 when you're done with it

            //don't free
            free(allEdges2);



            // Surface* new_surface1 = Surface_create(new_surface1_edges->data, new_surface1_edges->size, new_surface1_pointts, count1);
            // Surface* new_surface2 = Surface_create(new_surface2_edges->data, new_surface2_edges->size, new_surface2_pointts, count2);

            // for (int i=0;i<new_surface1_edges->size;++i) {
            //     printf("new_surface1_edges: %d \n", new_surface1_edges->data[i]->number);
            // }
            // for (int i=0;i<new_surface2_edges->size;++i) {
            //     printf("new_surface2_edges: %d \n", new_surface2_edges->data[i]->number);
            // }

            // for (int i=0;i<new_surface1->pointt_count;++i) {
            //     printf("new_surfaces_in_1: %d \n", new_surface1->pointts[i]->number);
            // }



            //  for (int i=0;i<new_surface2->pointt_count;++i) {
            //     printf("new_surfaces_in_2: %d \n", new_surface2->pointts[i]->number);
            // }
            add_surface_to_vector(new_surfaces, new_surface1);
            add_surface_to_vector(new_surfaces, new_surface2);

            if (new_surfaces == NULL || new_surfaces->data == NULL) {
                printf("Null pointter encountered.\n");
                return;
            }
            // printf("111\n");
            
            // printf("111 new_surface2->pointt_count: %d \n", new_surfaces->data[1]->pointt_count);
            // for (int i=0; i<new_surfaces->data[1]->pointt_count; ++i) {
            //     printf("111 new_surface2: %d\n",new_surfaces->data[1]->pointts[i]->number);
            // }
            // printf("111 new_surface2->edge_count: %d \n", new_surfaces->data[1]->edge_count);
            // for (int i=0;i<new_surfaces->data[1]->edge_count;++i) {
            //     printf("111 new_surface2_edges: %d \n", new_surfaces->data[1]->edges[i]->number);
            // }

            free(new_surface1_pointts);
            free(new_surface2_pointts); 
            free(new_surface1_edges->data);
            free(new_surface1_edges);
            free(overlaps1->data);
            free(overlaps1);
            free(new_surface2_edges->data);
            free(new_surface2_edges);
            free(overlaps2->data);
            free(overlaps2);
        

        free(original_pointts);
          
    }

    // Clean up
    if (intersections) {
        free(intersections);
    }

    // return resulting_surfaces;
}


void cut_surfaces(SurfaceVector* surfaces, Line* cutting_line, SurfaceVector* new_surfaces) {
    for (int i = 0; i < surfaces->size; ++i) {
        Surface* current_surface = surfaces->data[i];
        SurfaceVector* cut_result = init_surface_vector(2);
        cut_surface(current_surface, cutting_line, cut_result);
        for (int j = 0; j < cut_result->size; ++j) {
            add_surface_to_vector(new_surfaces, cut_result->data[j]);
        }
        // printf("surface_print i=%d\n",i);
        // SurfaceVector_print(cut_result);
        free(cut_result->data);  // assuming cut_surface dynamically allocates memory for cut_result
        free(cut_result); 
    }
}



void SurfaceVector_print(SurfaceVector* surfaceVector) {
    if (surfaceVector == NULL) {
        printf("SurfaceVector is NULL.\n");
        return;
    }
    printf("SurfaceVector size: %d, capacity: %d\n", surfaceVector->size, surfaceVector->capacity);
    for (int i = 0; i < surfaceVector->size; ++i) {
        // printf("Surface %d:\n", i + 1);
        Surface_print(surfaceVector->data[i]);
    }
}

// Function to check whether a pointt is above a line
int Pointt_above_line(double pointt_coordinates[2], Line* line) {
    double result = line->a * pointt_coordinates[0] + line->b * pointt_coordinates[1] + line->c;
    if (result > 0) return 1;
    if (result == 0) return 0;
    return -1;
}

// Function to get the center of a surface based on the coordinates of its pointts
void Surface_get_center(Surface* surface, double center_coordinates[2]) {
    double x_sum = 0;
    double y_sum = 0;
    for (int i = 0; i < surface->pointt_count; ++i) {
        x_sum += surface->pointts[i]->coordinates[0];
        y_sum += surface->pointts[i]->coordinates[1];
    }
    center_coordinates[0] = x_sum / surface->pointt_count;
    center_coordinates[1] = y_sum / surface->pointt_count;
}



//type get
void Surface_get_type(Surface* surface, Line* line1, Line* line2) {
    double center[2];
    Surface_get_center(surface, center);

    //line1 is line between solid and fluid, ab of line1 is the normal out solid
    //line2 is line between gas and liquid, ab of line2 is the normal out liquid

    // Evaluate the position of the center relative to line1 and line2
    int above_line1 = Pointt_above_line(center, line1);  // Assuming you have a function that returns 1 if above, 0 if on, -1 if below
    int above_line2 = Pointt_above_line(center, line2);

    int type_temp=-1;
    // Determine the type based on the position
    if (above_line1 == 1 && above_line2 == 1) {
        type_temp = 1; // not solid and gas
    }
    else if (above_line1 == 1 && above_line2 == -1) {
        type_temp = 2; //not solid and liquid
    }
    else if (above_line1 == -1 && above_line2 == 1) {
        type_temp = 3; //solid and gas
    }
    else if (above_line1 == -1 && above_line2 == -1) {
        type_temp = 4; //solid and liquid
    }

    surface->type = type_temp;
}

//type get
void Surface_get_type_3_cut(Surface* surface, Line* line1, Line* line2, Line* line3) {
    double center[2];
    Surface_get_center(surface, center);

    //line1 is line between solid and fluid, ab of line1 is the normal out solid
    //line2 is line between gas and liquid, ab of line2 is the normal out liquid
    //line3 is face line parallel moving, ab of line2 is the normal out face

    // Evaluate the position of the center relative to line1 and line2
    int above_line1 = Pointt_above_line(center, line1);  // Assuming you have a function that returns 1 if above, 0 if on, -1 if below
    int above_line2 = Pointt_above_line(center, line2);
    int above_line3 = Pointt_above_line(center, line3);

    int type_temp=-1;
    // Determine the type based on the position
    if (above_line1 == 1 && above_line2 == 1 && above_line3==1) {
        type_temp = 1; // not solid and gas, moving part
    }
    else if (above_line1 == 1 && above_line2 == -1 && above_line3==1) {
        type_temp = 2; //not solid and liquid, moving part
    }
    else if (above_line1 == -1 && above_line2 == 1 && above_line3==1) {
        type_temp = 3; //solid and gas, moving part
    }
    else if (above_line1 == -1 && above_line2 == -1 && above_line3==1) {
        type_temp = 4; //solid and liquid, moving part
    }
    else if (above_line1 == 1 && above_line2 == 1 && above_line3==-1) {
        type_temp = 5; // not solid and gas, remaining part
    }
    else if (above_line1 == 1 && above_line2 == -1 && above_line3==-1) {
        type_temp = 6; //not solid and liquid, remaining part
    }
    else if (above_line1 == -1 && above_line2 == 1 && above_line3==-1) {
        type_temp = 7; //solid and gas, remaining part
    }
    else if (above_line1 == -1 && above_line2 == -1 && above_line3==-1) {
        type_temp = 8; //solid and liquid, remaining part
    }

    surface->type = type_temp;
}



//free

void free_pointt(Pointt* pointt) {
    // Free the edges array if needed
    if (pointt->edges != NULL) {
        free(pointt->edges);
    }

    // Free the surfaces array if needed
    if (pointt->surfaces != NULL) {
        free(pointt->surfaces);
    }

    // Free the pointt itself
    free(pointt);
}

void free_surface(Surface* surface) {
    // Free the edges array
    if (surface->edges != NULL) {
        free(surface->edges);
    }

    // Free the pointts array
    if (surface->pointts != NULL) {
        free(surface->pointts);
    }

    // Free the surface itself
    free(surface);
}

void free_edge(Edge* edge) {
    // Free the surfaces array if needed
    if (edge->surfaces != NULL) {
        free(edge->surfaces);
    }

    // Free the edge itself
    free(edge);
}



void free_all_pointts(Pointt** all_pointts, int pointt_count) {
    for (int i = 0; i < pointt_count; ++i) {
        free_pointt(all_pointts[i]);
    }
    // free(all_pointts);
}

void free_all_edges(Edge** all_edges, int edge_count) {
    for (int i = 0; i < edge_count; ++i) {
        free_edge(all_edges[i]);
    }
    // free(all_edges);
}

void free_all_surfaces(Surface** all_surfaces, int surface_count) {
    for (int i = 0; i < surface_count; ++i) {
        free_surface(all_surfaces[i]);
    }
    // free(all_surfaces);
}

void free_all_lines(Line** all_lines, int line_count) {
    for (int i = 0; i < line_count; ++i) {
        free(all_lines[i]);
    }
    // free(all_lines);
}

void free_all_resources() {
    free_all_pointts(all_pointts, pointt_count);
    free_all_edges(all_edges, edge_count);
    free_all_surfaces(all_surfaces, surface_count);
    free_all_lines(all_lines, line_count);
}


// void line_intersection_range(double a, double b, double* c_min, double* c_max) {
void line_intersection_range(coord n1, double* c_min, double* c_max) {
    // Vertices of the unit square
    double a = n1.x;
    double b = n1.y;
    double vertices[4][2] = {{0,0}, {1,0}, {0,1}, {1,1}};
    
    // Initialize c_min and c_max
    *c_min = 1.0e10;  // a very large positive number
    *c_max = -1.0e10; // a very large negative number
    
    // Check the special cases of vertical and horizontal lines
    if (a == 0 && b == 0) {
        // The line is a point and doesn't intersect the square
        *c_min = *c_max = 0;
        return;
    }
    
    // Compute c for each vertex and update c_min and c_max
    for(int i=0; i<4; i++) {
        double x = vertices[i][0];
        double y = vertices[i][1];
        double c = -a*x - b*y;
        if (c < *c_min) *c_min = c;
        if (c > *c_max) *c_max = c;
    }
}

// double binary_search_for_c(double a, double b, double c_min, double c_max, double target_area, double tol, int max_iter) {
//     int iter = 0;
//     double c_mid, area_mid;
    
//     while(iter < max_iter) {
//         c_mid = (c_min + c_max) / 2.0;
//         area_mid = area(a, b, c_mid);
        
//         if (fabs(area_mid - target_area) < tol) {
//             return c_mid; // Found c satisfying the condition
//         }
        
//         if ((area_mid < target_area) == (area(a, b, c_min) < target_area)) {
//             c_min = c_mid;
//         } else {
//             c_max = c_mid;
//         }
        
//         iter++;
//     }
    
//     return -1; // Indicate that c was not found within max_iter iterations
// }