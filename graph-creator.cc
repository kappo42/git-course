#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <list>

using namespace std;

#define pi 3.14159265358979323846

#define R 6371
#define TO_RAD (3.1415926536 / 180)

double dist(double th1, double ph1, double th2, double ph2)
{
	double dx, dy, dz;
	ph1 -= ph2;
	ph1 *= TO_RAD, th1 *= TO_RAD, th2 *= TO_RAD;
    
	dz = sin(th1) - sin(th2);
	dx = cos(ph1) * cos(th1) - cos(th2);
	dy = sin(ph1) * cos(th1);
	return asin(sqrt(dx * dx + dy * dy + dz * dz) / 2) * 2 * R;
}

// MAP for the graph
// <int,list<int>> (point of interest, list of points attacched to him)
typedef std::map<int, list<int> > graphcityMap;
graphcityMap grcM;
graphcityMap::iterator it_grcM;

// The structure location identifies geographicla position of users and consequently samples
struct location{
    int id;
    double lat;                  // Latitute  - Geo
    double lon;                  // Longitude - Geo
    double alt;                  // Altitude  - Geo
    int pos;
    
    location(int A, double B, double C, double D): id(A), lat(B), lon(C), alt(D){}
};

// List to include the file
list<location> locationsL;
list<location>::iterator it_locationsL;


/*********************************************************************************************/
// delimiter to read the file and global definition for vector of strings
char const row_delim = '\n';
char const field_delim = '\t';
typedef vector<vector<string> > Rows;
Rows rows;



   int main(){
     
   ifstream infile( "./map-points-center-lux.txt");
   ofstream graphfile; // this file prints the graph, only ids with the list of adjacencies
   graphfile.open("map-graph-points.txt");
   ofstream associationfile; // this file only prints
   associationfile.open("map-association.txt");
   
       int id_point=0;

   // in this file there's no dummy line as first line
   
   for (string row; getline(infile, row, row_delim); ) {

      // 1: longitude, 2: latitude, 3: altitude
       
      rows.push_back(Rows::value_type());
      istringstream ss(row);
      for (string field; getline(ss, field, field_delim); ) {
          rows.back().push_back(field);
      }
    
      id_point++;

      location temp_loc(id_point,atof(rows.back().at(1).c_str()),atof(rows.back().at(0).c_str()),atof(rows.back().at(2).c_str()));
      locationsL.push_back(temp_loc);
      
   }
   if (!infile.eof()){
   cerr << "Error in reading the file 'map-points-center-lux.txt'!"<<endl;
   }
   infile.close();
   
   // printing association file
   for (it_locationsL= locationsL.begin(); it_locationsL != locationsL.end(); ++it_locationsL) {
       associationfile<<(*it_locationsL).id<<" "<<(*it_locationsL).lat<<" "<<(*it_locationsL).lon<<" "<<(*it_locationsL).alt<<endl;
   }
       
   associationfile.close();
       
   // clustering
   for (it_locationsL= locationsL.begin(); it_locationsL != locationsL.end(); ++it_locationsL) {
      //it_locationsL=locationsL.begin();
      list<location>::iterator it_tmplocL;
       
      list<int> list_adjacencies;// temporary list in which we insert all points closest to the one considered
      
      for (it_tmplocL= locationsL.begin(); it_tmplocL != locationsL.end(); ++it_tmplocL) {
          
          // we compute the distance from the reference point it_locationsL and the current under
          // consideration it_tmplocL
          double meters   = dist((*it_locationsL).lat,(*it_locationsL).lon,(*it_tmplocL).lat,(*it_tmplocL).lon)*1000;
          //cout<<fixed<<meters<<"\t"<<(*it_tmplocL).lat<<", "<<(*it_tmplocL).lon<<endl;
          
          // check if this value is below our threshold 50 meters
          if(meters<50.00 && ((*it_locationsL).id)!=((*it_tmplocL).id)){
              // insertion in  map of adjacencies
              list_adjacencies.push_back((*it_tmplocL).id); 
          }
         else return; 
      }
       // actual insertion in map
       grcM.insert(pair <int, list <int> >((*it_locationsL).id,list_adjacencies));
   }
    
   // printing the map
    for (it_grcM= grcM.begin(); it_grcM != grcM.end(); ++it_grcM) {
        list<int>::iterator it_tmplst;
        list<int> temp_list=(*it_grcM).second;
        
        graphfile<<(*it_grcM).first<<" ";
        for (it_tmplst= temp_list.begin(); it_tmplst != temp_list.end(); ++it_tmplst) {
            graphfile<<(*it_tmplst)<<" ";
        }
        graphfile<<endl;
    }
       graphfile.close();
   
   return 0;
   }

