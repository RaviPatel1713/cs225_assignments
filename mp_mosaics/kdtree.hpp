/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <cmath>
using namespace std;

template<int Dim> 
void KDTree<Dim>::swap(Point<Dim>& point1, Point<Dim>& point2)
{
  Point<Dim> temp = point1;
  point1 = point2;
  point2 = temp;
}

template<int Dim> 
void KDTree<Dim>::partition(vector<Point<Dim>>& points, 
                            size_t left, size_t right, 
                            int curr_dim, 
                            size_t partition_idx)
{
  size_t pivot_idx = right;
  size_t i = left;
  for (size_t j = i; j < right; ++j) {
    if (smallerDimVal(points[j], points[pivot_idx], curr_dim)) {
      swap(points[i], points[j]); 
      ++i;
    }
  }
  swap(points[pivot_idx], points[i]);
  pivot_idx = i;

  if (partition_idx < pivot_idx) {
    partition(points, left, pivot_idx-1, curr_dim, partition_idx);
  }
  else if (partition_idx > pivot_idx) {
    partition(points, pivot_idx+1, right, curr_dim, partition_idx);
  }
}

template <int Dim> 
void KDTree<Dim>::build_KDTree(vector<Point<Dim>>&points, 
                               KDTreeNode*& subroot, 
                               long int left, long int right, int d) 
{
  if (left > right) { return; } 

  int curr_dim = d % Dim;
  size_t median = (left + right)/2;
  partition(points, left, right, curr_dim, /*point to partition around*/median); 
  
  subroot = new KDTreeNode(points[median]);

  build_KDTree(points, subroot->left, left, median - 1, d + 1);
  build_KDTree(points, subroot->right, median + 1, right, d + 1);
}

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
    if (std::fabs(first[curDim] - second[curDim]) < 0.0000001f) {
      // breaks tie if the value along a dimension is same 
      return first < second;
    }
    return first[curDim] < second[curDim];
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    // euclidean distance approach
    double curr_best_dist = 0;
    for (size_t i = 0; i < Dim; ++i) 
    { curr_best_dist += std::pow(target[i] - currentBest[i], 2); }

    double potential_best_dist = 0;
    for (size_t i = 0; i < Dim; ++i) 
    { potential_best_dist += std::pow(target[i] - potential[i], 2); }

    // ties broken with Point::operator<() method 

    if (std::fabs(curr_best_dist - potential_best_dist) < 0.0000001f) {
      return potential < currentBest;
    }

    // found potential to be closer than current best neighbor 
    return potential_best_dist < curr_best_dist;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
    root = NULL;
    size = newPoints.size();
    if (!newPoints.empty()) {
      // copy newPoints contents into points
      vector<Point<Dim>> points(newPoints);
      size_t left = 0;
      size_t right = points.size() - 1;
      int d = 0;
      build_KDTree(points, root, left, right, d);
    }
}

template <int Dim> 
void KDTree<Dim>::_copy(KDTreeNode*& this_tree_subroot, const KDTreeNode*& other_tree_subroot) {
  if (other_tree_subroot != NULL) {
    this_tree_subroot = new KDTreeNode(other_tree_subroot->point);
    _copy(this_tree_subroot->left, other_tree_subroot->left);
    _copy(this_tree_subroot->right, other_tree_subroot->right);
  }
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
  root = NULL;
  size = 0;
  _copy(root, other->root);
  this->size = other->size;
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */ 
  if (this != &rhs) {
    _destroy(root);
    _copy(root, rhs->root);
    this->size = rhs->size;
  }
  return *this;
}


template <int Dim> 
void KDTree<Dim>::_destroy(KDTreeNode*& subroot) {
  if (subroot != NULL) {
    _destroy(subroot->left);
    _destroy(subroot->right);
    delete subroot;
    subroot = NULL;
  }
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  _destroy(root);
  size = 0;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    if (root != NULL) {
      // initiallize with some point (for it to be passed as reference)
      Point<Dim> closest = root->point;
      closest = findNearestNeighbor(root, query, /*depth*/0, closest);
      return closest;
    }
    return Point<Dim>();
}

template <int Dim>
const Point<Dim>& KDTree<Dim>::findNearestNeighbor(KDTreeNode* subroot, const Point<Dim>& target, int d, 
                                                   Point<Dim>& closest) const
{
  if (subroot == NULL) { return closest; }

  int curr_dim = d % Dim;
  KDTreeNode*  subtree_dir = subroot->left;
  bool went_left = true;

  if (smallerDimVal(subroot->point, target, curr_dim)) {
    subtree_dir = subroot->right;
    went_left = false;
  }

  // check if current subroot node point is closer than current best found 
  if (shouldReplace(target, findNearestNeighbor(subtree_dir, target, d+1, closest), subroot->point)) {
    closest = subroot->point;
  }

  double curr_best_dist = 0;
  for (size_t i = 0; i < Dim; ++i) { curr_best_dist += std::pow(target[i] - closest[i], 2); }

  if (std::pow(target[curr_dim] - subroot->point[curr_dim], 2) < curr_best_dist
  || std::fabs(std::pow(target[curr_dim] - subroot->point[curr_dim], 2) - curr_best_dist) < 0.0000001f) {
    subtree_dir = went_left ? subroot->right : subroot->left;
    Point<Dim> potential_closest = findNearestNeighbor(subtree_dir, target, d+1, closest);    

    if (shouldReplace(target, closest, potential_closest)) {
      closest = potential_closest;
    }
  }
  
  return closest;
}
