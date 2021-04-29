// [[Rcpp::depends(lidR)]]

#include <Rcpp.h>
#include <SpatialIndex.h>
using namespace Rcpp;
using namespace lidR;

// [[Rcpp::export]]
Rcpp::XPtr<lidR::QuadTree> quadtree(S4 las) {
  QuadTree* qt = new QuadTree(las);
  Rcpp::XPtr<QuadTree> p(qt, true);
  return p;
}

// [[Rcpp::export]]
IntegerVector filter_orectangle_with_index(SEXP xptr, double xmin, double ymin, double xmax, double ymax, double angle) {

  Rcpp::XPtr<QuadTree> tree(xptr);
  OrientedRectangle orect(xmin, ymin, xmax, ymax, angle);
  std::vector<PointXYZ> pts;
  tree->lookup(orect, pts);

  IntegerVector ids(pts.size());
  for(int i = 0 ; i < pts.size(); i++)
    ids[i] = pts[i].id;

  return ids+1;
}
