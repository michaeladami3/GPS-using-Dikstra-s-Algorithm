/*main.cpp*/

//
// Original Author:Prof. Joe Hummel
// Revised Author: Michael Adami
// U. of Illinois, Chicago
// CS 251: Spring 2020
// Project #07: open street maps, graphs, and Dijkstra's alg
// 
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:  
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include "graph.h"
#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
using namespace std;
using namespace tinyxml2;
//
// Dijkstra:
//
// Performs Dijkstra's shortest weighted path algorithm from
// the given start vertex.  Returns a vector of vertices in
// the order they were visited, along with a map of (string,int)
// pairs where the string is a vertex V and the int is the 
// distance from the start vertex to V; if no such path exists,
// the distance is INF (defined in util.h).
//
 class prioritize{ //For the priority_queue
    public:
    bool operator()(const pair<long long,double> &p1, const pair<long long,double> &p2) const{
       if(p1.second<p2.second) //Check for the weight
         return false;
       else if(p1.second>p2.second) //Check for the weight
         return true;
       else
         return p2.first<p1.first; //If it was smaller weight
    }
 };

     
//Solving the if the paths of nodes can be reached using Dijkstra's algorithm, returns a map of nodes by IDs
stack<long long> Dijkstra(graph<long long,double>& G, long long startV, long long endV, map<long long, double>& distances,double& tDistance,map<long long,long long>& predecessors,double INF)
{   
   double edgeWeight=0.0; //Intialize
   double alternativePathDistance; 
   set<long long> discoveredSet;
   stack<long long> path; //For the smallest path
  priority_queue<pair<long long,double>,vector<pair<long long,double>>,prioritize> pq;
   for (long long v:G.getVertices())
   {
      distances[v]=INF; //Every distance=INF
      predecessors[v]=0; // Every predecessors=0
      pq.push(make_pair(v,INF)); //Save those pairs
    }
   distances[startV]=0;
   pq.push(make_pair(startV,0));
   while(!pq.empty()){ 
      auto elem = pq.top(); 
      pq.pop();
      if(elem.second==INF) //Check that it wasnt the end
         break;
      else if(discoveredSet.count(elem.first)) //Check the count in the set
         continue;
      else
         discoveredSet.insert(elem.first);  
      for(long long v:G.neighbors(elem.first)){
         G.getWeight(elem.first,v,edgeWeight);
         alternativePathDistance=distances[elem.first]+edgeWeight;
         if(alternativePathDistance<distances[v]){
            distances[v]=alternativePathDistance;  // Save the distance to it
            predecessors[v]=elem.first; //Save the first element
            tDistance=distances[endV]-distances[startV];
            pq.push(make_pair(v,alternativePathDistance));
         }
      }
   }
    //Start the stack 
   path.push(endV); //Push the first
   while(predecessors[endV]){ //Push elements in while it isnt ENDV
       endV=predecessors[endV];
       path.push(endV);
   }
   return path;
}
/* Check: Uses input the insure that the user has put in valid input and will search for that building in the vector of  buildings
 * then will search for that building and save the information needed to output into main. 
 */ 
 void check(vector<BuildingInfo> Buildings,string startName,string endName,string &sName,string &eName,double &sLat,double &sLong,double &desLat,double &desLong,bool& check1, bool& check2){
          for(BuildingInfo v:Buildings){
             if(v.Abbrev==startName){ 
              sName=v.Fullname; //Save the name
              sLat=v.Coords.Lat; //Save the lat
              sLong=v.Coords.Lon; //Save the long
              check1=true;
             }
          if(!check1){
              if(v.Fullname.find(startName)!=string::npos){
              sName=v.Fullname; //Save the name
              sLat=v.Coords.Lat; //Save the lat
              sLong=v.Coords.Lon; // Save the Long
              check1=true;
            }
          }
        }
     if(!check1) //If it was an abbreviation and it didn't get any information
         return;
      for(BuildingInfo v:Buildings){
        if(v.Abbrev==endName){ // If an abbreviation
            eName=v.Fullname; //Save the name
            desLat=v.Coords.Lat; //Set the coords
            desLong=v.Coords.Lon;//Set the coords
            check2=true;
        }
       if(!check2) //If it was an abbreviation
           if(v.Fullname.find(endName)!=string::npos){
            eName=v.Fullname; //Save the name
            desLat=v.Coords.Lat; //Set the coords
            desLong=v.Coords.Lon;//Set the coords
            check2=true;
     }
   }
 }

//AddAllVerticesandEdges:Create the graph
void AddAllVerticesandEdges(vector<FootwayInfo> Footways,map<long long, Coordinates> Nodes,graph<long long, double>& G){
    for(auto v: Nodes){ //ADD ALL THE NODES
           G.addVertex(v.first);
    }
    for(auto v:Footways){ //Look through Footways
        for(size_t i=0;i<v.Nodes.size()-1;i++){ // Get the map and look through that
            double distance = distBetween2Points(Nodes[v.Nodes[i]].Lat,Nodes[v.Nodes[i]].Lon,
                                                Nodes[v.Nodes[i+1]].Lat,Nodes[v.Nodes[i+1]].Lon); //This is getting the distance between the nodes
            G.addEdge(v.Nodes[i],v.Nodes[i+1],distance); // Add the first one
            G.addEdge(v.Nodes[i+1],v.Nodes[i],distance); //On the way back
        }
        
    }
    
}

/*NodesNearby:Searches through all the Nodes and checks if it the closest distance because the buildings have no node coordinated to it.*/
void NodesNearby(vector<FootwayInfo> Footways,map<long long, Coordinates> Nodes,double startLat,double startLong,double desLat,double desLong,long long& sID,long long& dID){
        double prevSDist,prevDDist,sLong,sLat,dLat,dLong,sDistance,dDistance; //The elements for the function
        for(auto v:Footways){ //Loop through each footway
            prevSDist=distBetween2Points(startLat,startLong,
                                         Nodes.at(v.Nodes.at(0)).Lat,Nodes.at(v.Nodes.at(0)).Lon); //Get the distance of each node from the start. We need something to initialize compare
            prevDDist=distBetween2Points(desLat,desLong, 
                                        Nodes.at(v.Nodes.at(0)).Lat,Nodes.at(v.Nodes.at(0)).Lon);//Get the distance of each node from the end. We need something to initialize compare
            } 
        for(auto v: Footways){ //Get each Node
            for(size_t i = 0; i<v.Nodes.size();i++){ //Get the vector
            sDistance=distBetween2Points(startLat,startLong,Nodes.at(v.Nodes.at(i)).Lat,Nodes.at(v.Nodes.at(i)).Lon); //Search the distance between start and node
            dDistance=distBetween2Points(desLat,desLong,Nodes.at(v.Nodes.at(i)).Lat,Nodes.at(v.Nodes.at(i)).Lon); // Search the distance between destination and node again
            if(sDistance<prevSDist){ // if they have a distance less than before store is as the new placeholder
                sLong=Nodes.at(v.Nodes.at(i)).Lon;
                sLat=Nodes.at(v.Nodes.at(i)).Lat;
                prevSDist=sDistance; //Set the distance
                sID=Nodes.at(v.Nodes.at(i)).ID; //Start ID
            }
            if(dDistance<prevDDist){ //If they have a distance less than before make that the new placeholder
                dLong=Nodes.at(v.Nodes.at(i)).Lon;
                dLat=Nodes.at(v.Nodes.at(i)).Lat;
                prevDDist=dDistance; //set the distance
                dID=Nodes.at(v.Nodes.at(i)).ID; //END ID
            }
         } 
        }
        cout<<"Nearest start node:"<<endl;
        cout<<" "<<sID<<endl; //Start ID
        cout<<" ("<<sLat<<", "<<sLong<<")"<<endl;
        cout<<"Nearest destination node:"<<endl;
        cout<<" "<<dID<<endl; //END ID
        cout<<" ("<<dLat<<", "<<dLong<<")"<<endl;
}
//pathUsingDijkstra:Function that outputs the path to the other node using Dijkstra
void pathUsingDijkstra(graph<long long,double>& G,long long startV, long long endV){
    map<long long, long long> predecessors; 
    map<long long,double> distancesBetweenNodes; 
    stack<long long> path; //I use this to print path
    const double INF = numeric_limits<double>::max(); //INF
    double tDistance=0; //Total Distance
    long long visitedV; //Saved for the top of the stack
   stack<long long> visited=Dijkstra(G,startV,endV,distancesBetweenNodes,tDistance,predecessors,INF); //Call the functions and save to stack
    if(startV!=endV){ //If it is the same Building
        if( distancesBetweenNodes[endV]!=INF && !visited.empty()){ //If reachable 
          cout<<"Distance to dest: "<<tDistance<<" miles"<<endl;
          cout<<"Path: ";
           while(!visited.empty()){ //Get ready to output
               visitedV=visited.top();
               visited.pop();
               cout<<visitedV;
               if(!visited.empty()) // Make sure we take the arrow from out front
                cout<<"->";
            else
                cout<<" ";
           }
            cout<<endl;
               }
            else
               cout<<"Sorry, destination unreachable"<<endl;
        predecessors.clear(); //Clear the maps
        distancesBetweenNodes.clear(); //Clear the maps
    }
    else{ //It is the same thing
        cout<<"Distance to dest: "<< tDistance <<" miles"<<endl;
        cout<<"Path: ";
        cout<<startV<<endl;
    }
}

//////////////////////////////////////////////////////////////////
//
// main
//
int main()
{ 
  string startBuilding, destBuilding,sName,eName;
  map<long long, Coordinates>  Nodes;     // maps a Node ID to it's coordinates (lat, lon)
  vector<FootwayInfo>          Footways;  // info about each footway, in no particular order
  vector<BuildingInfo>         Buildings; // info about each building, in no particular order
  XMLDocument                  xmldoc;
  bool                         checker1; //Check startBuilding
  bool                         checker2; //check destBuilding
  graph<long long, double> G; // New Graph
  vector<long long> visited;
  map<long long, double> distances; 
  double startLat,startLong,desLat,desLong; //LAT/LONG
  long long sID,dID; // STARTBUILDING ID/ENDBUILDING ID
  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "")
  {
    filename = def_filename;
  }

  //
  // Load XML-based map file 
  //
  if (!LoadOpenStreetMap(filename, xmldoc))
  {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }
  
  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);
  
  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);
  set<long long>  S;
  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == Nodes.size());
  assert(footwayCount == Footways.size());
  assert(buildingCount == Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  AddAllVerticesandEdges(Footways,Nodes,G);
  cout << "# of buildings: " << Buildings.size() << endl;
  cout<<  "# of vertices: " << G.NumVertices()<<endl;
  cout<<  "# of edges: " << G.NumEdges()<<endl;
  cout << endl;

  //
  // Navigation from building to building
  //

/*Start the input for the buildings. If input isn't in the vector. Checks the abbreviation first
 * and then checks the substrings. If they are neither then continue on to the next two inputs.
 */
  cout << "Enter start (partial name or abbreviation), or #> ";
  getline(cin, startBuilding);

  while (startBuilding != "#")
  {   
    checker1=false; //New check
    checker2=false; // New check
    cout << "Enter destination (partial name or abbreviation)> ";

    getline(cin, destBuilding);
    check(Buildings,startBuilding,destBuilding,sName,eName,startLat,startLong,desLat,desLong,checker1,checker2);//Get all the information
    if(!checker1|| !checker2){ //IF WE HAVE A PROBLEM 
        if(!checker1) //If this passed Start Building Inspection
            cout<<"Start building not found"<<endl;
        else if(!checker2) //If this passed End Building Inspection
            cout<<"Destination building not found"<<endl;
    }
    else{
    cout<<"Starting point:"<<endl; 
    cout<<" "<<sName<<endl;
    cout<<" ("<<startLat<<", "<<startLong<<")"<<endl;
    cout<<"Destination point:"<<endl;
    cout<<" "<<eName<<endl; //from the function "check"
    cout<<" ("<<desLat<<", "<<desLong<<")"<<endl<<endl;
    NodesNearby(Footways,Nodes,startLat,startLong,desLat,desLong,sID,dID); //Nodes nearby output just in a function
    cout<<endl<<"Navigating with Dijkstra..."<<endl;     
    pathUsingDijkstra(G, sID,dID); //Call to output just in function
    }

    //
    // another navigation?
    //
    cout << endl;
    cout << "Enter start (partial name or abbreviation), or #> ";
    getline(cin, startBuilding);
  }
  //
  // done:
  //
  cout << "** Done **" << endl;

  return 0;
}