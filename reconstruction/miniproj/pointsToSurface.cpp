#include "pointsToSurface.h"
#include <iostream>

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


void PointsToSurface::computeNonOrientedNormals() {

    int k = 3;
    double A11,A12,A13,A22,A23,A33;
    for(int i=0;i<_points.size();i++){

        Point3D p = _points[i];


        //compute the k-neighbors : pts
        v_Point3D pts = kneighborhoodPoints(p,_points,k);

        //compute center of gravity B of pts
        Point3D b = center_of_gravity(pts);

        //compute the matrix A = pts - B
        Point3D a1 = pts[0]-b;
        Point3D a2 = pts[1]-b;
        Point3D a3 = pts[2]-b;

        //compute A' * A = [A11 A12 A13; A12 A22 A23; A13 A23 A33]
        //A12 == A21 , A23 == A32, A13 == A31 because scalar product is symetric

        A11 = produit_scalaire(a1,a1);
        A12 = produit_scalaire(a1,a2);
        A13 = produit_scalaire(a1,a3);

        A22 = produit_scalaire(a2,a2);
        A23 = produit_scalaire(a2,a3);

        A33 = produit_scalaire(a3,a3);

        //compute the eigen vectors of A'*A
        Point3D u;
        Point3D v;
        Point3D n;


        //Problem here ?
        calcul_repere_vecteurs_propres(A11,A12,A13,A22,A23,A33,u,v,n);


        //return the non oriented normals _noNormals

        Point3D norm = produit_vectoriel(v,u);
        _noNormals.push_back( n);
    }
}

/*Si l'arc (i,j) est present dans le graphe g alors retourne vrai
* sinon retourne faux*/
bool presentDansArbre(Graphe g,int i){

    int k=0;
    while(k < g.nb_arcs() ){
        if(g.arc(k).n1 == i || g.arc(k).n2 == i){

            return true;
        }
        k++;
    }
    return false;
}

bool estLie(Graphe g,int i){
    return false;
}

void PointsToSurface::computeMinimalSpanningTree() {
    //Pour chaque point
    //d(Pi,Pj) <r et i!=j
    float r = 0.5;
    Point3D pi;
    Point3D pj;
    Graphe arbre = Graphe(_points.size());
    for(int i=0; i< _points.size();i++){
        pi = _points[i];
        //Si on a déjà regardé pi on fait rien
        //Sinon
        //Pour tous les points à une distance r de pi pj
        //on ajoute l'edge pi pj
        for(int j=0;j< _points.size();j++){
            pj = _points[j];
            //Si i = j ou d(pi,pj) >= r ou pi,pj deja present dans l'arbre
            if((i == j) || (distance_(pi,pj) >= r) || presentDansArbre(arbre,j)){
                //On passe
            }else{
                //On ajoute l'edge (pi,pj)
                arbre.ajouter_arc(i,j);
            }
        }
    }
    _acm = arbre.arbre_couvrant_minimal();
}

void PointsToSurface::computeOrientedNormals() {
  // a remplir : _oNormals
}

double PointsToSurface::computeImplicitFunc(double x,double y,double z) {
  // a faire : déterminer la fonction implicite (MLS)

  return x+y+z;
}

void PointsToSurface::computeNormalsFromImplicitFunc() {
  // a remplir : _surfacen
}

void PointsToSurface::computeMesh() {
  // a remplir : _surfacep
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
	
  // lecture du nombre de points et du rayon d'échantillonnage
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
