#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <utility>

// Simple container describing a time-based cluster built from adjacent channels.
class Cluster {
public:
    std::vector<int> channels;       // Contributing channels (relative indices)
    double time = 0.0;               // Average leading time across hits in the cluster
    double timeOverThreshold = 0.0;  // Sum of TOTs for all hits
    int centerChannel = -1;          // Channel having the highest TOT inside the cluster

    Cluster() = default;
    Cluster(std::vector<int> channelsIn, double clusterTime, double tot, int center)
        : channels(std::move(channelsIn)), time(clusterTime), timeOverThreshold(tot), centerChannel(center) {}
};

#endif  // CLUSTER_H
