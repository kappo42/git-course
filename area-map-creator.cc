#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <map>
#include <list>
#include <iomanip>      // std::setprecision

using namespace std;

struct location{
    float xval;
    float yval;
    
    location(float A, float B): xval(A), yval(B){}
};

bool operator < (const location &l, const location &r) {
    
return l.xval < r.xval;
}// to build the map

/*********************************************************************************************/
// delimiter to read the file and global definition for vector of strings
char const row_delim = '\n';
char const field_delim = ' ';
typedef vector<vector<string> > Rows;
Rows rows;

int a=6378137;
float e2=0.000669;
float radius=0.015;

/*********************************************************************************************/
// DEFINITION OF POLYGONS COMPOSING AREAS

int area_one_sides=6;
float  area_oneLAT[]={49.6127,49.6151,49.6168,49.6138,49.6133,49.6122};
float  area_oneLON[]={6.12168,6.12361,6.12618,6.12778,6.12597,6.12583};

int area_two_sides=6;
float  area_twoLAT[]={49.6127,49.6122,49.6095,49.6074,49.6098,49.6114};
float  area_twoLON[]={6.12168,6.12583,6.12625,6.12179,6.12065,6.1201};

int area_three_sides=4;
float  area_threeLAT[]={49.6122,49.6095,49.6096,49.612};
float  area_threeLON[]={6.12583,6.12625,6.1296,6.12965};

int area_four_sides=10;
float  area_fourLAT[]={49.6145,49.6153,49.6155,49.6149,49.6127,49.6117,49.612,49.6122,49.6133,49.6138};
float  area_fourLON[]={6.12735,6.12696,6.12838,6.13027,6.13348,6.13414,6.12965,6.12583,6.12597,6.12778};

int area_five_sides=5;
float  area_fiveLAT[]={49.6096,49.612,49.6117,49.6106,49.6052};
float  area_fiveLON[]={6.1296,6.12965,6.13414,6.13421,6.13365};

// PRINT POLYGON DEFINING THE AREA
void printArea(std::ofstream& dumFile, int sides, float *areaLAT, float *areaLON){
    for (int i=0; i<sides; i++) {
        float lat=areaLAT[i];
        float lon=areaLON[i];
        
        float gamma=a/sqrt(1-e2*pow(sin(lat),2));
        //float x= ((gamma+atof(rows.back().at(2).c_str()))*cos(lat)*cos(lon))/50000;
        //float y= ((gamma+atof(rows.back().at(2).c_str()))*cos(lat)*sin(lon))/50000;
        
        float x= ((gamma)*cos(lat)*cos(lon))/50000;
        float y= ((gamma)*cos(lat)*sin(lon))/50000;
        
        if(i<(sides-1)){
            dumFile<<"("<<x<<","<<y<<")--";
        }else{
            dumFile<<"("<<x<<","<<y<<");";
        }
    }
}

// COMPARING FLOAT
#define EPSILON 0.0001
bool CompareFloat (float A, float B)
{
    float diff = A - B;
    return (diff < EPSILON) && (-diff < EPSILON);
}

/*********************************************************************************************/
// DETECTING IF A POINT BELONGS TO A GIVEN AREA
// solution from http://stackoverflow.com/questions/4287780/detecting-whether-a-gps-coordinate-falls-within-a-polygon-on-a-map
double PI = 3.14159265;
double TWOPI = 2*PI;

double Angle2D(double y1, double x1, double y2, double x2)
{
    double dtheta,theta1,theta2;
    
    theta1 = atan2(y1,x1);
    theta2 = atan2(y2,x2);
    dtheta = theta2 - theta1;
    while (dtheta > PI)
        dtheta -= TWOPI;
    while (dtheta < -PI)
        dtheta += TWOPI;
    
    return(dtheta);
}

// checks whether a point(latitude,longitude) is inside an area
bool coordinate_is_inside_area(double latitude, double longitude, float *lat_array, float *long_array,int sides)
{
    int i;
    double angle=0;
    double point1_lat;
    double point1_long;
    double point2_lat;
    double point2_long;
    int n = sides;
    
    for (i=0;i<n;i++) {
        point1_lat = lat_array[i] - latitude;
        point1_long = long_array[i] - longitude;
        point2_lat = lat_array[(i+1)%n] - latitude;
        point2_long = long_array[(i+1)%n] - longitude;
        angle += Angle2D(point1_lat,point1_long,point2_lat,point2_long);
    }
    
    if (abs(angle) < PI)
        return false;
    else
        return true;
}
 
// MAPS with samples per area
std::map<int,int> area_one_samplesM;
std::map<int,int>::iterator it_areaonesamples;

std::map<int,int> area_two_samplesM;
std::map<int,int>::iterator it_areatwosamples;

std::map<int,int> area_three_samplesM;
std::map<int,int>::iterator it_areathreesamples;

std::map<int,int> area_four_samplesM;
std::map<int,int>::iterator it_areafoursamples;

std::map<int,int> area_five_samplesM;
std::map<int,int>::iterator it_areafivesamples;

// =====>>>>> Just one iterator is sufficient!!!!!


int  time_intervals[12]={830,900,930,1000,1030,1100,1130,1200,1230,1300,1330,1400};

// Counting the number of samples per time interval   =====>>>>> Very likely to be redundant with next
void collectStatsTimeMap(std::map< int, int >  &tmpM, int time_c){
    if(time_c>800 && time_c <=830){
        if(tmpM.count(830)==0){
            tmpM[830]=1;
        }else{
            tmpM[830]+=1;
        }
    }
    if(time_c>830 && time_c <=900){
        if(tmpM.count(900)==0){
            tmpM[900]=1;
        }else{
            tmpM[900]+=1;
        }
    }
    if(time_c>900 && time_c <=930){
        if(tmpM.count(930)==0){
            tmpM[930]=1;
        }else{
            tmpM[930]+=1;
        }
    }
    if(time_c>930 && time_c <=1000){
        if(tmpM.count(1000)==0){
            tmpM[1000]=1;
        }else{
            tmpM[1000]+=1;
        }
    }
    if(time_c>1000 && time_c <=1030){
        if(tmpM.count(1030)==0){
            tmpM[1030]=1;
        }else{
            tmpM[1030]+=1;
        }
    }
    if(time_c>1030 && time_c <=1100){
        if(tmpM.count(1100)==0){
            tmpM[1100]=1;
        }else{
            tmpM[1100]+=1;
        }
    }
    if(time_c>1100 && time_c <=1130){
        if(tmpM.count(1130)==0){
            tmpM[1130]=1;
        }else{
            tmpM[1130]+=1;
        }
    }
    if(time_c>1130 && time_c <=1200){
        if(tmpM.count(1200)==0){
            tmpM[1200]=1;
        }else{
            tmpM[1200]+=1;
        }
    }
    if(time_c>1200 && time_c <=1230){
        if(tmpM.count(1230)==0){
            tmpM[1230]=1;
        }else{
            tmpM[1230]+=1;
        }
    }
    if(time_c>1230 && time_c <=1300){
        if(tmpM.count(1100)==0){
            tmpM[1300]=1;
        }else{
            tmpM[1300]+=1;
        }
    }
    if(time_c>1300 && time_c <=1330){
        if(tmpM.count(1130)==0){
            tmpM[1330]=1;
        }else{
            tmpM[1330]+=1;
        }
    }
    if(time_c>1330 && time_c <=1400){
        if(tmpM.count(1200)==0){
            tmpM[1400]=1;
        }else{
            tmpM[1400]+=1;
        }
    }
}


// Haversine distance
#define pi 3.14159265358979323846

#define R 6371
#define TO_RAD (3.1415926536 / 180)

double havdist(double th1, double ph1, double th2, double ph2)
{
	double dx, dy, dz;
	ph1 -= ph2;
	ph1 *= TO_RAD, th1 *= TO_RAD, th2 *= TO_RAD;
    
	dz = sin(th1) - sin(th2);
	dx = cos(ph1) * cos(th1) - cos(th2);
	dy = sin(ph1) * cos(th1);
	return asin(sqrt(dx * dx + dy * dy + dz * dz) / 2) * 2 * R;
}

std::vector<location> area_one_sV;
std::vector<location> area_two_sV;
std::vector<location> area_three_sV;
std::vector<location> area_four_sV;
std::vector<location> area_five_sV;

std::map<int, std::vector<location> > area_one_sM;
std::map<int, std::vector<location> > area_two_sM;
std::map<int, std::vector<location> > area_three_sM;
std::map<int, std::vector<location> > area_four_sM;
std::map<int, std::vector<location> > area_five_sM;
std::map<int, std::vector<location> >::iterator it_areasm;

// Storing samples per time
void collectSamplesTimeVector(std::map<int, std::vector<location> >  &tmpM, int time, location loc){
    if(time>800 && time <=830){
        if(tmpM.count(830)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(830,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(830);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(830);
            tmpM.insert(pair <int, std::vector<location> >(830,tmpV));
        }
    }
    if(time>830 && time <=900){
        if(tmpM.count(900)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(900,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(900);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(900);
            tmpM.insert(pair <int, std::vector<location> >(900,tmpV));
        }
    }
    if(time>900 && time <=930){
        if(tmpM.count(930)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(930,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(930);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(930);
            tmpM.insert(pair <int, std::vector<location> >(930,tmpV));
        }
    }
    if(time>930 && time <=1000){
        if(tmpM.count(1000)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1000,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1000);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1000);
            tmpM.insert(pair <int, std::vector<location> >(1000,tmpV));
        }
    }
    if(time>1000 && time <=1030){
        if(tmpM.count(1030)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1030,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1030);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1030);
            tmpM.insert(pair <int, std::vector<location> >(1030,tmpV));
        }
    }
    if(time>1030 && time <=1100){
        if(tmpM.count(1100)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1100,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1100);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1100);
            tmpM.insert(pair <int, std::vector<location> >(1100,tmpV));
        }
    }
    if(time>1100 && time <=1130){
        if(tmpM.count(1130)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1130,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1130);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1130);
            tmpM.insert(pair <int, std::vector<location> >(1130,tmpV));
        }
    }
    if(time>1130 && time <=1200){
        if(tmpM.count(1200)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1200,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1200);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1200);
            tmpM.insert(pair <int, std::vector<location> >(1200,tmpV));
        }
    }
    if(time>1200 && time <=1230){
        if(tmpM.count(1230)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1230,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1230);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1230);
            tmpM.insert(pair <int, std::vector<location> >(1230,tmpV));
        }
    }
    if(time>1230 && time <=1300){
        if(tmpM.count(1300)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1300,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1300);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1300);
            tmpM.insert(pair <int, std::vector<location> >(1300,tmpV));
        }
    }
    if(time>1300 && time <=1330){
        if(tmpM.count(1330)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1330,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1330);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1330);
            tmpM.insert(pair <int, std::vector<location> >(1330,tmpV));
        }
    }
    if(time>1330 && time <=1400){
        if(tmpM.count(1400)==0){
            // we are sure not having samples so far, so we create the vector, include the new sample and that's it
            std::vector<location> tmpV;
            tmpV.push_back(loc);
            // include it in the map
            tmpM.insert(pair <int, std::vector<location> >(1400,tmpV));
        }else{
            // here we have first to retrieve the correct vector from the map and then push_back the new value
            std::map<int, std::vector<location> >::iterator it_tmpM;
            it_tmpM=tmpM.find(1400);
            std::vector<location> tmpV=(*it_tmpM).second;
            
            tmpV.push_back(loc);
            
            // first erase the record then re-insert it again
            tmpM.erase(1400);
            tmpM.insert(pair <int, std::vector<location> >(1400,tmpV));
        }
    }
}

void printAreaSamplesPerTime(std::map<int, std::vector<location> >  &tmpM){
    std::map<int, std::vector<location> >::iterator it_tmpM;
    // Printing Map with samples grouped per time sample
    for (it_tmpM = tmpM.begin(); it_tmpM != tmpM.end(); ++it_tmpM) {
        cout<<(*it_tmpM).first<<": ";
        std::vector<location> tmpV=(*it_tmpM).second;
        std::vector<location>::iterator it_tmpv;
        
        /*for (it_tmpv = tmpV.begin(); it_tmpv != tmpV.end(); ++it_tmpv) {
         cout<<(*it_tmpv).xval<<" "<<(*it_tmpv).yval<<" - ";
         }
         */
        cout<<tmpV.size();
        cout<<endl;
        
    }
    
}

// Computing the metric of Sample Distribution Coefficient

std::map<int, float > area_one_metricM;
std::map<int, float > area_two_metricM;
std::map<int, float > area_three_metricM;
std::map<int, float > area_four_metricM;
std::map<int, float > area_five_metricM;

std::map<int, float >::iterator it_aremetric;

float computeSampleDistributionCoefficient(std::map<int, std::vector<location> >  &tmpM, int time_c){
    std::map<int, std::vector<location> >::iterator it_areasm;
    it_areasm=tmpM.find(time_c);
    std::vector<location> tmpV=(*it_areasm).second;
    
    std::vector<location>::iterator it_tmpvi;
    std::vector<location>::iterator it_tmpvj;
    
    // counter to be sure remove zeros from the list
    int k=0;
    double dist=0.0;
    for (it_tmpvi = tmpV.begin(); it_tmpvi != tmpV.end()-1; ++it_tmpvi) {
        for (it_tmpvj = it_tmpvi+1; it_tmpvj != tmpV.end(); ++it_tmpvj) {
            // distance in meters
            double dist_ij=havdist((*it_tmpvi).xval,(*it_tmpvi).yval,(*it_tmpvj).xval,(*it_tmpvj).yval)*1000;
            if(dist_ij>0){// just to avoid distances = 0
                dist+=dist_ij;
                k++;
            }
        }
    }
    float avg_sample_dist=dist/k; // average distance between any two samples in the area - similar to average path length in graph theory
    float sample_dist_coeff=(*it_areasm).second.size()/avg_sample_dist;
    
    return sample_dist_coeff;
    //cout<<"Average Sample Distance: "<<avg_sample_dist<<" Sample Distribution Coefficient: "<<sample_dist_coeff<<" sample/m"<<endl;
}
   int main(){
     
   ifstream infile( "../heatmapendsim.dat");

   // in this file there's no dummy line as first line
   
   for (string row; getline(infile, row, row_delim); ) {

      // 1: value, 2: latitude, 3: longitude, 4: timestamp
       
      rows.push_back(Rows::value_type());
      istringstream ss(row);
      for (string field; getline(ss, field, field_delim); ) {
          rows.back().push_back(field);
      }
       
       int time=atoi(rows.back().at(3).c_str());
       //cout<<rows.back().at(1)<<"\t"<<rows.back().at(2)<<endl;

       double lt=atof(rows.back().at(1).c_str());
       double ln=atof(rows.back().at(2).c_str());
       
       location loc(lt,ln);
       
       if(coordinate_is_inside_area(lt,ln,area_oneLAT,area_oneLON,area_one_sides)==true){
           // statistics per time
           collectStatsTimeMap(area_one_samplesM,time);
           
           // collect statistics to compute the metric "sample distribution"
           collectSamplesTimeVector(area_one_sM, time,loc);
       }else if(coordinate_is_inside_area(lt,ln,area_twoLAT,area_twoLON,area_two_sides)==true){
           // statistics per time
           collectStatsTimeMap(area_two_samplesM,time);
           
           // collect statistics to compute the metric "sample distribution"
           collectSamplesTimeVector(area_two_sM, time,loc);
       }else if(coordinate_is_inside_area(lt,ln,area_threeLAT,area_threeLON,area_three_sides)==true){
           // statistics per time
           collectStatsTimeMap(area_three_samplesM,time);
           // collect statistics to compute the metric "sample distribution"
           collectSamplesTimeVector(area_three_sM, time,loc);

       }else if(coordinate_is_inside_area(lt,ln,area_fourLAT,area_fourLON,area_four_sides)==true){
           // statistics per time
           collectStatsTimeMap(area_four_samplesM,time);
           // collect statistics to compute the metric "sample distribution"
           collectSamplesTimeVector(area_four_sM, time,loc);

       }else if(coordinate_is_inside_area(lt,ln,area_fiveLAT,area_fiveLON,area_five_sides)==true){
           // statistics per time
           collectStatsTimeMap(area_five_samplesM,time);
           // collect statistics to compute the metric "sample distribution"
           collectSamplesTimeVector(area_five_sM, time,loc);
       }

   }
   if (!infile.eof()){
   cerr << "Error in reading the file 'map-points-center-lux.txt'!"<<endl;
   }
   infile.close();
       
   /*
   // Printing Map with samples grouped per time sample
   cout<<" Area 1"<<endl;
   printAreaSamplesPerTime(area_one_sM);
   cout<<" Area 2"<<endl;
   printAreaSamplesPerTime(area_two_sM);
   cout<<" Area 3"<<endl;
   printAreaSamplesPerTime(area_three_sM);
   cout<<" Area 4"<<endl;
   printAreaSamplesPerTime(area_four_sM);
   cout<<" Area 5"<<endl;
   printAreaSamplesPerTime(area_five_sM);
   */
       
       
   // printing areas values + finding the maximum # of samples generated for comparison and map coloring later
   float max_over_interval=0.0;
   cout<<endl;
   for (it_areaonesamples = area_one_samplesM.begin(); it_areaonesamples != area_one_samplesM.end(); ++it_areaonesamples) {
       if((*it_areaonesamples).second>max_over_interval){
           max_over_interval=(*it_areaonesamples).second;
       }
       //cout<<(*it_areaonesamples).first<<" "<<(*it_areaonesamples).second<<endl;
   }
    //cout<<endl;
    for (it_areatwosamples = area_two_samplesM.begin(); it_areatwosamples != area_two_samplesM.end(); ++it_areatwosamples) {
        if((*it_areatwosamples).second>max_over_interval){
            max_over_interval=(*it_areatwosamples).second;
        }
        //cout<<(*it_areatwosamples).first<<" "<<(*it_areatwosamples).second<<endl;
    }
    //cout<<endl;
    for (it_areathreesamples = area_three_samplesM.begin(); it_areathreesamples != area_three_samplesM.end(); ++it_areathreesamples) {
        if((*it_areathreesamples).second>max_over_interval){
            max_over_interval=(*it_areathreesamples).second;
        }
        //cout<<(*it_areathreesamples).first<<" "<<(*it_areathreesamples).second<<endl;
    }
    //cout<<endl;
    for (it_areafoursamples = area_four_samplesM.begin(); it_areafoursamples != area_four_samplesM.end(); ++it_areafoursamples) {
        if((*it_areafoursamples).second>max_over_interval){
            max_over_interval=(*it_areafoursamples).second;
        }
        //cout<<(*it_areafoursamples).first<<" "<<(*it_areafoursamples).second<<endl;
    }
    //cout<<endl;
    for (it_areafivesamples = area_five_samplesM.begin(); it_areafivesamples != area_five_samplesM.end(); ++it_areafivesamples) {
        if((*it_areafivesamples).second>max_over_interval){
            max_over_interval=(*it_areafivesamples).second;
        }
        //cout<<(*it_areafivesamples).first<<" "<<(*it_areafivesamples).second<<endl;
    }
   cout<<"Max over period: "<<max_over_interval<<endl<<endl;
       
   //computeSampleDistributionCoefficient();
   
   // * * * * * * * * *
   // PRINTING OVER FILE
       
   ofstream outfile,coordfile;
   outfile.open("map.tex");
   coordfile.open("coordinates-map.tex");
       
   outfile<<"\\documentclass[tikz,border=10pt]{standalone}"<<endl<<endl;
    // color specifications
   outfile<<"\\definecolor{mapblue}{HTML}{0000FF}%"<<endl;
   outfile<<"\\definecolor{maproyalblue}{HTML}{0088FF}%"<<endl;
   outfile<<"\\definecolor{mapskyblue}{HTML}{00CCFF}%"<<endl;
   outfile<<"\\definecolor{mapgreen}{HTML}{00CC00}%"<<endl;
   outfile<<"\\definecolor{mapdandelion}{HTML}{FFCC00}%"<<endl;
   outfile<<"\\definecolor{maporange}{HTML}{FF8800}%"<<endl;
   outfile<<"\\definecolor{mapred}{HTML}{FF0000}%"<<endl<<endl;
       
       
   outfile<<"\\usepackage{pgfplots}"<<endl;
   outfile<<"\\usetikzlibrary{backgrounds,shadings}"<<endl;
       
   outfile<<"% #1=name, #2 shading width, #3=list of colors"<<endl;
   outfile<<"\\tikzset{relevance map/.code={"<<endl;
   outfile<<"\\pgfdeclareverticalshading{relevancemapshading}{2cm}{"<<endl;
   outfile<<"color(0cm)=(mapblue);"<<endl;
   outfile<<"color(0.675cm)=(mapskyblue);"<<endl;
   outfile<<"color(1cm)=(mapgreen);"<<endl;
   outfile<<"color(1.325cm)=(mapdandelion);"<<endl;
   outfile<<"color(1.9cm)=(mapred)"<<endl;
   outfile<<"}"<<endl;
   outfile<<"\\pgfkeysalso{/tikz/shading=relevancemapshading}"<<endl;
   outfile<<"}"<<endl;
   outfile<<"}"<<endl;
       
   outfile<<"\\pgfplotsset{"<<endl;
   outfile<<"colormap={rel-map}{[1cm] color(0cm)=(mapblue) color(1cm)=(mapskyblue) color(2cm)=(mapgreen) color(3cm)=(mapdandelion) color(5cm)=(mapred)}"<<endl;
   outfile<<"}"<<endl;
       
   outfile<<"\\begin{document}"<<endl;
   // here we should create a "foreach" that loops how many intervals of time we do have
   
   string old_hours="8";
   string old_minutes="00";
   for (int i=0; i<=11; i++) {
       int time_considered=time_intervals[i];
       it_areaonesamples = area_one_samplesM.find(time_considered);
       it_areatwosamples = area_two_samplesM.find(time_considered);
       it_areathreesamples = area_three_samplesM.find(time_considered);
       it_areafoursamples = area_four_samplesM.find(time_considered);
       it_areafivesamples = area_five_samplesM.find(time_considered);
       
       float val_area_one=(*it_areaonesamples).second/max_over_interval;
       float val_area_two=(*it_areatwosamples).second/max_over_interval;
       float val_area_three=(*it_areathreesamples).second/max_over_interval;
       float val_area_four=(*it_areafoursamples).second/max_over_interval;
       float val_area_five=(*it_areafivesamples).second/max_over_interval;
       
       /*
       cout<<"Time: "<<time_considered<<"\t "<<val_area_one<<" "
                                           <<val_area_two<<" "
                                           <<val_area_three<<" "
                                           <<val_area_four<<" "
                                           <<val_area_five<<endl;*/
       
       // ** Computing the metric Sample Distribution Coefficient
       
       //cout<<" Area 1"<<endl;
       float area_one_sdc=computeSampleDistributionCoefficient(area_one_sM,time_considered);
       //cout<<" Area 2"<<endl;
       float area_two_sdc=computeSampleDistributionCoefficient(area_two_sM,time_considered);
       //cout<<" Area 3"<<endl;
       float area_three_sdc=computeSampleDistributionCoefficient(area_three_sM,time_considered);
       //cout<<" Area 4"<<endl;
       float area_four_sdc=computeSampleDistributionCoefficient(area_four_sM,time_considered);
       //cout<<" Area 5"<<endl;
       float area_five_sdc=computeSampleDistributionCoefficient(area_five_sM,time_considered);
       
       
       outfile<<"\\begin{tikzpicture}[yscale=-1,rotate=-90,scale=5]"<<endl;
       outfile<<"\\input{city.tex}"<<endl;
       outfile<<"\\begin{pgfonlayer}{background}"<<endl;
       // area 1
       outfile<<"\\pgfplotscolormapaccess[0:1]{"<<val_area_one<<"}{rel-map}"<<endl;
       outfile<<"\\def\\TEMP{\\definecolor{my color}{rgb}}"<<endl;
       outfile<<"\\expandafter\\TEMP\\expandafter{\\pgfmathresult}"<<endl;
       outfile<<"\\fill[my color!40]"<<endl;
       printArea(outfile,area_one_sides,area_oneLAT,area_oneLON);
       outfile<<"\\fill[my color!40,ultra thin](99.5,-15.9)rectangle(99.45,-15.8);"<<endl;
       outfile<<"\\node [right,text=black,font=\\tiny]at(99.475,-15.8){Area 1: "<<area_one_sdc<<" sample/m};"<<endl;
       // area 2
       outfile<<"\\pgfplotscolormapaccess[0:1]{"<<val_area_two<<"}{rel-map}"<<endl;
       outfile<<"\\def\\TEMP{\\definecolor{my color}{rgb}}"<<endl;
       outfile<<"\\expandafter\\TEMP\\expandafter{\\pgfmathresult}"<<endl;
       outfile<<"\\fill[my color!40]"<<endl;
       printArea(outfile,area_two_sides,area_twoLAT,area_twoLON);
       outfile<<"\\fill[my color!40,ultra thin](99.4,-15.9)rectangle(99.35,-15.8);"<<endl;
       outfile<<"\\node [right,text=black,font=\\tiny]at(99.375,-15.8){Area 2: "<<area_two_sdc<<" sample/m};"<<endl;
       // area 3
       outfile<<"\\pgfplotscolormapaccess[0:1]{"<<val_area_three<<"}{rel-map}"<<endl;
       outfile<<"\\def\\TEMP{\\definecolor{my color}{rgb}}"<<endl;
       outfile<<"\\expandafter\\TEMP\\expandafter{\\pgfmathresult}"<<endl;
       outfile<<"\\fill[my color!40]"<<endl;
       printArea(outfile,area_three_sides,area_threeLAT,area_threeLON);
       outfile<<endl;
       outfile<<"\\fill[my color!40,ultra thin](99.3,-15.9)rectangle(99.25,-15.8);"<<endl;
       outfile<<"\\node [right,text=black,font=\\tiny]at(99.275,-15.8){Area 3: "<<area_three_sdc<<" sample/m};"<<endl;
       // area 4
       outfile<<"\\pgfplotscolormapaccess[0:1]{"<<val_area_four<<"}{rel-map}"<<endl;
       outfile<<"\\def\\TEMP{\\definecolor{my color}{rgb}}"<<endl;
       outfile<<"\\expandafter\\TEMP\\expandafter{\\pgfmathresult}"<<endl;
       outfile<<"\\fill[my color!40]"<<endl;
       printArea(outfile,area_four_sides,area_fourLAT,area_fourLON);
       outfile<<"\\fill[my color!40,ultra thin](99.2,-15.9)rectangle(99.15,-15.8);"<<endl;
       outfile<<"\\node [right,text=black,font=\\tiny]at(99.175,-15.8){Area 4: "<<area_four_sdc<<" sample/m};"<<endl;
       // area 5
       outfile<<"\\pgfplotscolormapaccess[0:1]{"<<val_area_five<<"}{rel-map}"<<endl;
       outfile<<"\\def\\TEMP{\\definecolor{my color}{rgb}}"<<endl;
       outfile<<"\\expandafter\\TEMP\\expandafter{\\pgfmathresult}"<<endl;
       outfile<<"\\fill[my color!40]"<<endl;
       printArea(outfile,area_five_sides,area_fiveLAT,area_fiveLON);
       outfile<<"\\fill[my color!40,ultra thin](99.1,-15.9)rectangle(99.05,-15.8);"<<endl;
       outfile<<"\\node [right,text=black,font=\\tiny]at(99.075,-15.8){Area 5: "<<area_five_sdc<<" sample/m};"<<endl;
       // end areas
       outfile<<endl<<"\\end{pgfonlayer}"<<endl;
       outfile<<"\\draw[relevance map,ultra thin](99.5,-16.5)rectangle(99.2,-16.4);"<<endl;
       outfile<<"\\node[right,font=\\tiny, scale=0.8] at (99.5,-16.4) {"<<max_over_interval<<"};"<<endl;
       outfile<<"\\node[right,font=\\tiny, scale=0.8] at (99.2,-16.4) {0};"<<endl;
       
       // transform to time
       string hours,minutes;
       int time=time_considered;
       std::stringstream out;
       out << time;
       string s_time=out.str();
       
       if(s_time.size()==3){
           hours.append(s_time.substr(0,1));
           minutes.append(s_time.substr(1,2));
       }
       if(s_time.size()==4){
           hours.append(s_time.substr(0,2));
           minutes.append(s_time.substr(2,2));
       }
       

       
       outfile<<"\\node[below,font=\\footnotesize, scale=0.8,minimum width=2cm] at (99.1,-16.4) {Time: "<<old_hours<<":"<<old_minutes<<"-"<<hours<<":"<<minutes<<"};"<<endl;
       outfile<<"\\end{tikzpicture}"<<endl;
       
       old_hours=hours;
       old_minutes=minutes;
       
   }
   /*
   outfile<<"\\begin{tikzpicture}[yscale=-1,rotate=-90,scale=5]"<<endl;
   outfile<<"\\input{city.tex}"<<endl;
       
   outfile<<"\\begin{pgfonlayer}{background}"<<endl;
    for (it_stat= statM.begin(); it_stat != statM.end(); ++it_stat) {
        float val=(float) (*it_stat).second/max;
        
        outfile<<"\\pgfplotscolormapaccess[0:1]{"<<val<<"}{rel-map}"<<endl;
        outfile<<"\\def\\TEMP{\\definecolor{my color}{rgb}}"<<endl;
        outfile<<"\\expandafter\\TEMP\\expandafter{\\pgfmathresult}"<<endl;
        outfile<<"\\fill[my color!40]"<<endl;
        
        if((*it_stat).first==1){
            printArea(outfile,area_one_sides,area_oneLAT,area_oneLON);
        }else if((*it_stat).first==2){
            printArea(outfile,area_two_sides,area_twoLAT,area_twoLON);
        }else if((*it_stat).first==3){
            printArea(outfile,area_three_sides,area_threeLAT,area_threeLON);
        }else if((*it_stat).first==4){
            printArea(outfile,area_four_sides,area_fourLAT,area_fourLON);
        }else if((*it_stat).first==5){
            printArea(outfile,area_five_sides,area_fiveLAT,area_fiveLON);
        }
    }
       
    outfile<<endl<<"\\end{pgfonlayer}"<<endl;
       
    outfile<<"\\draw[relevance map,ultra thin](99.5,-16.5)rectangle(99.2,-16.4);"<<endl;
    outfile<<"\\node[right,font=\\tiny, scale=0.8] at (99.5,-16.4) {"<<sum/5<<"};"<<endl;
    outfile<<"\\node[right,font=\\tiny, scale=0.8] at (99.2,-16.4) {0};"<<endl;
    
    outfile<<"\\end{tikzpicture}"<<endl;
    */
    outfile<<"\\end{document}"<<endl;
    outfile.close();
       
       area_one_samplesM.clear();
       area_two_samplesM.clear();
       area_three_samplesM.clear();
       area_four_samplesM.clear();
       area_five_samplesM.clear();
    
   return 0;
   }

