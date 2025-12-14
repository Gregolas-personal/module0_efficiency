#ifndef CLUSTER_UTILS_H
#define CLUSTER_UTILS_H

#include <vector>
#include "cluster.h"
#include "layer.h"

namespace cluster_utils {

struct LayerClusters {
    std::vector<Cluster> eta1;
    std::vector<Cluster> eta2;
};

std::vector<Cluster> buildClusters(const std::vector<Layer::Hit>& hits);

}  // namespace cluster_utils

#endif  // CLUSTER_UTILS_H
