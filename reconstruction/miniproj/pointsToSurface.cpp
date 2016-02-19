#include "pointsToSurface.h"
using namespace std;

PointsToSurface::PointsToSurface(const QString &filename) 
  : _acm(0) {
  setlocale(LC_ALL,"C");

  if(!readInputFile(filename)) {
    cout << "Error while loading" << filename.toStdString() << endl;
    exit(0);
  }
}

PointsToSurface::~PointsToSurface() {
  _points.clear();
}

//Return the center of gravity of the points contained in the vector pts
Point3D PointsToSurface::center_of_gravity(v_Point3D pts){
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    for(int i=0;i<pts.size();i++){
        x += pts[i].x;
        y += pts[i].y;
        z += pts[i].z;
    }
    x = x/pts.size();
    y = y/pts.size();
    z = z/pts.size();

    return Point3D(x,y,z);
}

//Return the scalar product of the vectors u and v , both having the same size "size"
double scalar_product(vector<double> u, vector<double> v, int size){

    double res= 0.0;
    for(int i=0;i<size;i++){
        res = res + u[i]*v[i];
    }
    return res;
}

void PointsToSurface::computeNonOrientedNormals() {

    nb_neighbors = 6;
    double A11,A12,A13,A22,A23,A33;
    for(int i=0;i<_points.size();i++){

        Point3D p = _points[i];


        //compute the k-neighbors : pts
        v_Point3D pts = kneighborhoodPoints(p,_points,nb_neighbors);

        //compute center of gravity B of pts
        Point3D b = center_of_gravity(pts);

        //compute the matrix A = pts - B

        v_Point3D A(nb_neighbors);
        for(int i=0;i<pts.size();i++){
            A[i] = pts[i]-b;
        }
        //A[i] is the i-th column of A

        vector<double> a1(nb_neighbors) ;
        vector<double> a2(nb_neighbors) ;
        vector<double> a3(nb_neighbors) ;
        for(int i=0;i<pts.size();i++){
            a1[i] = A[i].x;
            a2[i] = A[i].y;
            a3[i] = A[i].z;
        }
        //a1 is the first column of A , or the first row of A'
        //a2    second
        //a3    third

        //compute A' * A = [A11 A12 A13; A12 A22 A23; A13 A23 A33]
        //A12 == A21 , A23 == A32, A13 == A31 because scalar product is symetric

        A11 = scalar_product(a1,a1,nb_neighbors);
        A12 = scalar_product(a1,a2,nb_neighbors);
        A13 = scalar_product(a1,a3,nb_neighbors);

        A22 = scalar_product(a2,a2,nb_neighbors);
        A23 = scalar_product(a2,a3,nb_neighbors);

        A33 = scalar_product(a3,a3,nb_neighbors);

        //compute the eigen vectors of A'*A
        Point3D u;
        Point3D v;
        Point3D n;

        calcul_repere_vecteurs_propres(A11,A12,A13,A22,A23,A33,u,v,n);


        //return the non oriented normals _noNormals
        _noNormals.push_back(n);
    }
}

double PointsToSurface::compute_radius(){

    double max = -1;
    for(int i=0;i<_points.size();i++){
        v_Point3D neighbors= kneighborhoodPoints(_points[i],_points,2);
        Point3D neighbor = neighbors[1];
        double dist = distance_(neighbor,_points[i]);
        if(dist > max){
            max = dist;
        }
    }
    //We round it up as a double (example if r = 0.17 we have r = 0.2)
    return (ceil(max*10))/10;
}

void PointsToSurface::computeMinimalSpanningTree() {
    
    Point3D pi;
    Point3D pj;
    Graphe tree = Graphe(_points.size());

    //Choice of the radius
    double r = compute_radius();

    //Computing of proximity graph
    for(int i=0; i < _points.size();i++){
        pi = _points[i];
        for(int j=i+1; j< _points.size();j++){
            pj = _points[j];
            if( (i!=j) && distance_(pi,pj) < r){
                double weight =  1 - abs( produit_scalaire(_noNormals[i],_noNormals[j]) ) ;
                tree.ajouter_arc(i,j,weight);
            }
        }
    }

    //Computing of minimal spanning tree
   // _acm = tree;
    _acm = tree.arbre_couvrant_minimal();
}

void PointsToSurface::normalsRedirection(Graphe g,Noeud node,int s,QVector<bool> seen_nodes){

    l_int arcs_list = node.la;
    //for all the arcs going from/to node
    for (l_int::iterator k=arcs_list.begin();
     k!=arcs_list.end(); k++)
  {
    Arc a = g.arc(*k);
    unsigned int j = (a.n1 == s) ? a.n2 : a.n1;
    //Look at scalar product son/father
    //if <0 redirect the normals
    if(!seen_nodes.at(j)){

        seen_nodes[j] = true;

        Point3D normal_father = _oNormals[s];
        Point3D normal_son = _oNormals[j];

        if(produit_scalaire(normal_father,normal_son) < 0){
            _oNormals[j] = -normal_son;
        }else{
            _oNormals[j] = normal_son;
        }

        normalsRedirection(g,g.noeud(j),j,seen_nodes);
    }
  }


}

void PointsToSurface::computeOrientedNormals() {
  // a remplir : _oNormals
    QVector<bool> seen_nodes(_points.size(),false);
    seen_nodes[0] = true;
    _oNormals = _noNormals;
    normalsRedirection(_acm,_acm.noeud(0),0,seen_nodes);
}

double PointsToSurface::computeImplicitFunc(double x,double y,double z) {
  // a faire : determiner la fonction implicite (MLS)
  double sigma = 0.1;

  double sum = 0.0;
  double sum_weight = 0.0;
  double weight = 0.0;
  Point3D X = Point3D(x,y,z);

  for(int i=0;i<_points.size();i++){
      weight = exp(-pow( norme(X-_points[i])/sigma ,2 ));
      sum +=  ( produit_scalaire( _oNormals[i] , (X-_points[i]) ) * weight );
      sum_weight += weight;
  }

  return (sum/sum_weight);
}

void PointsToSurface::computeNormalsFromImplicitFunc() {

    _surfacen = _surfacep;
    double nx,ny,nz,x,y,z;
    for(int i=0;i<_surfacep.size();i++){
        Point3D s0 = _surfacep[i].S0;
        Point3D s1 = _surfacep[i].S1;
        Point3D s2 = _surfacep[i].S2;
        x = s0.x;
        y = s0.y;
        z = s0.z;

        nx = computeImplicitFunc(x-0.01,y,z) - computeImplicitFunc(x+0.01,y,z);
        ny = computeImplicitFunc(x,y-0.01,z) - computeImplicitFunc(x,y+0.01,z);
        nz = computeImplicitFunc(x,y,z-0.01) - computeImplicitFunc(x,y,z+0.01);
        Point3D ns0 = normalise(Point3D(nx,ny,nz));

        x = s1.x;
        y = s1.y;
        z = s1.z;

        nx = computeImplicitFunc(x-0.01,y,z) - computeImplicitFunc(x+0.01,y,z);
        ny = computeImplicitFunc(x,y-0.01,z) - computeImplicitFunc(x,y+0.01,z);
        nz = computeImplicitFunc(x,y,z-0.01) - computeImplicitFunc(x,y,z+0.01);
        Point3D ns1 = normalise(Point3D(nx,ny,nz));

        x = s2.x;
        y = s2.y;
        z = s2.z;

        nx = computeImplicitFunc(x-0.01,y,z) - computeImplicitFunc(x+0.01,y,z);
        ny = computeImplicitFunc(x,y-0.01,z) - computeImplicitFunc(x,y+0.01,z);
        nz = computeImplicitFunc(x,y,z-0.01) - computeImplicitFunc(x,y,z+0.01);
        Point3D ns2 = normalise(Point3D(nx,ny,nz));

        Triangle3D triangle;
        triangle.S0 = ns0;
        triangle.S1 = ns1;
        triangle.S2 = ns2;
        _surfacen[i] = triangle;
    }
}

void PointsToSurface::computeMesh() {
  // a remplir : _surfacep

    // Create a 3DGrid G
    //Covering box and step in the 3 dimensions

    Point3D min = boundingMin();
    Point3D max = boundingMax();
    unsigned int n_x = 9;
    unsigned int n_y = 9;
    unsigned int n_z = 9;
    Grille3D grille = Grille3D(min.x,min.y,min.z,max.x,max.y,max.z,n_x,n_y,n_z);

    int DIM_X = n_x+1;
    int DIM_Y = n_y+1;
    int DIM_Z = n_z+1;


  //Create an array v with the values of the implicit function
    //same size as the grid
    double v[DIM_X*DIM_Y*DIM_Z];
    //for each position of the 3DGrid compute the value of the implicit function

    for(int i=0;i<n_x+1;i++){
        for(int j=0;j<n_y+1;j++){
            for(int k=0;k<n_z+1;k++){
                v[i+DIM_X*(j+DIM_Y*k)] = computeImplicitFunc(grille.x(i),grille.y(j),grille.z(k));
            }
        }
    }


  //Compute the isovalue surface
    SurfaceIsovaleurGrille* isovalue_surface = new SurfaceIsovaleurGrille();
    isovalue_surface->surface_isovaleur(_surfacep,grille,v,0.0);
}






void PointsToSurface::computeSurface() {
  // appelle toutes les fonctions pour reconstruire la surface 
  cout << "computing non-oriented normals..." << endl;
  computeNonOrientedNormals();
  cout << "computing minimal spanning tree..." << endl;
  computeMinimalSpanningTree();
  cout << "computing oriented normals..." << endl;
  computeOrientedNormals();
  cout << "computing mesh surface..." << endl;
  computeMesh();
  cout << "computing final normals..." << endl;
  computeNormalsFromImplicitFunc();
}





v_Point3D  PointsToSurface::kneighborhoodPoints(const Point3D &p,const v_Point3D &pts,unsigned int k) const {
  // renvoie les k points les plus proches de p dans la liste de points pts

  v_Point3D neighbors;

  if(k==0) return neighbors;

  if(pts.size()<k) k=pts.size();
  neighbors.resize(k);
  
  std::vector<double> dist(pts.size());
  unsigned int imax = 0;

  for(unsigned int i=0;i<pts.size();++i) {
    dist[i]=distance_(p,pts[i]);
    
    if(dist[i]>dist[imax]) 
      imax=i;
  }


  for(unsigned int i=0;i<k;++i) {
    unsigned int index=0;
    
    for(unsigned int j=0;j<pts.size();++j) {
      if(dist[j]<dist[index]) 
	index=j;
    }

    neighbors[i]=pts[index];
    dist[index] = dist[imax]+1;
  }

  return neighbors;
}

bool PointsToSurface::readInputFile(const QString &filename) {

  _points.clear();
  
  FILE *f;
  unsigned int nS;
  double r;
  
  // ouverture du fichier
  if ((f = fopen(filename.toStdString().c_str(),"r"))==(FILE *)NULL)
    return false;
	
  // lecture du nombre de points et du rayon d'�chantillonnage
  if(fscanf(f, "%i %lf\n", &nS, &r)==EOF) return false;
	
  // lecture des points
  _points.resize(nS);
  for (unsigned int i=0; i<nS; i++) {
    if(fscanf(f, " %lf %lf %lf\n", &_points[i].x, &_points[i].y, &_points[i].z)==EOF) return false;
  }

  boite_englobante(_points,_boundingBox[0],_boundingBox[1]);

  // mean of the min distance between 2 points
  _meanDist = 0.0;
  for(unsigned int i=0;i<_points.size();++i) {
    v_Point3D pts = kneighborhoodPoints(_points[i],_points,2);
    _meanDist = _meanDist+distance_(pts[0],pts[1]);
  }
  _meanDist = _meanDist/(double)_points.size();

  // fermeture du fichier
  fclose(f);
	
  return true;
}
