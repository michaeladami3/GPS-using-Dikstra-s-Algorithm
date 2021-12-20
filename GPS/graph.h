/*graph.h*/

//
// Michael Adami
//
// Basic graph class using adjacency list representation using a map<pair<VertexT,VertexT>,weightT>.
//  We can store vectors and weights dynamically. We store the vertices inside the set and the edges
//  inside of the map.
//
// original author: Prof. Joe Hummel
// U. of Illinois, Chicago
// CS 251: Spring 2020
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <limits>
#include <set>
#include <map>
#include <queue>
#include <iterator>
using namespace std;
//Start of the graph set. 
template<typename VertexT, typename WeightT>
class graph
{
private:
    int count; //Get the count of elements while we input the edges. No need to make a function
    map<pair<VertexT,VertexT>,WeightT> adjacencyList; //Edges data
    set<VertexT> vertices; //Vertices 
    graph<long long, double> G();
    bool _LookupVertex(VertexT v){ //Check if they exist
        if(vertices.count(v)){ //Check the count inside of the set
            return true;
        }
        return false;
    }
public:
  //
  // constructor:
  //
  // Makes count zero to count the number of elements inputed.
  //
  graph(){
      count=0;//Initialize them to zero
  }
    
  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const
  {
    return vertices.size(); //Get the size of the set
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const
  {
    return count; //Get the count from private
  }

  //
  // addVertex
  //
  //Adds vertices to the list. If the edge exists, update the list and return true. If the vertices
  //exist and edge does not exist return true and add to list. Else return false because of non-existant
  //verticies. FIXED DYNAMICALLY
  //
  bool addVertex(VertexT v)
  {
    //
    // is the vertex already in the graph?  If so, we do not 
    // insert again otherwise Vertices may fill with duplicates:
    //
    if (_LookupVertex(v))
      return false;
    //
    // if we get here, vertex does not exist so insert.  Where
    // we insert this vertex in the adjacency list.
    //
    vertices.insert(v);
    return true;
  }

  //
  // addEdge
  //
  //Adds edges to the list. If the edge exists, update the list and return true. If the vertices
  //exist and edge does not exist return true and add to list. Else return false because of non-existant
  //verticies. FIXED DYNAMICALLY
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight)
  {
    if(!_LookupVertex(from)||!_LookupVertex(to)) //If they dont exist return
        return false;
    if(adjacencyList.find({from,to})!=adjacencyList.end()){ //If they already exist update and return
        adjacencyList[{from,to}]=weight;
        return true;
        }
        adjacencyList.insert(pair<pair<VertexT,VertexT>,WeightT>({from,to},weight));//Else they are new and they exist
        count++;
        return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If 
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not 
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight)const 
  {
    if(adjacencyList.find({from,to})!=adjacencyList.end()){
        weight=adjacencyList.find({from,to})->second;
        return true;
    }
    return false;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const
  {
    set<VertexT>  S; //New set
    //Searches a previous set and checks if it was inside the map and if it is in the map. We got the neighbors
    for(VertexT k:vertices) //Search using for each
       if(adjacencyList.find({v,k})!=adjacencyList.end()){ //If we haven't found it in the map
         S.insert(k); //insert it
   }
    return S; //Return the set from the neighbors
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const
  {
    vector<VertexT> verts; //New vector
    for(VertexT v:vertices){
        verts.push_back(v); //push everything from the set
    }
    return verts;  // returns a copy:
  }

  //
  // dump
  // 
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const
  {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    int i=0;
    for(VertexT v: vertices)
    {
      output << " " << i << ". " << v << endl;
      i++;
    }

    output << endl;
    output << "**Edges:" << endl;
    for(auto it = adjacencyList.cbegin(); it != adjacencyList.cend(); ++it)
    {
        cout<<"("<<it->first.first<<", "<<it->first.second<<", "<<it->second<<")"<<endl;
    }
    output << "**************************************************" << endl;
  }

};