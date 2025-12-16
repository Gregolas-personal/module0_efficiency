#include "cluster_utils.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace cluster_utils {

std::vector<Cluster> buildClusters(const std::vector<Layer::Hit>& hits, int layerIndex, int etaSide) {
    std::vector<Cluster> clusters;
    if (hits.empty()) {
        return clusters;
    }

    std::vector<Layer::Hit> sortedHits = hits;
    std::sort(sortedHits.begin(), sortedHits.end(), [](const Layer::Hit& a, const Layer::Hit& b) {
        if (a.channel == b.channel) {
            return a.leadingTime < b.leadingTime;
        }
        return a.channel < b.channel;
    });

    Cluster currentCluster;
    int hitCount = 0;
    int lastChannel = 0;
    float lastTime = 0.0F;
    float highestTot = -1.0F;
    int centerChannel = -1;
    float earliestTime = std::numeric_limits<float>::max();
    int earliestChannel = -1;
    float earliestTot = 0.0F;

    auto flushCluster = [&]() {
        if (hitCount == 0) {
            return;
        }
        currentCluster.time = earliestTime;
        currentCluster.timeOverThreshold = earliestTot;
        currentCluster.centerChannel = centerChannel;
        currentCluster.layer = layerIndex;
        currentCluster.etaSide = etaSide;
        if (earliestChannel >= 0 && centerChannel >= 0 && earliestChannel != centerChannel) {
            std::cerr << "Warning: cluster earliest channel " << earliestChannel
                      << " differs from highest TOT channel " << centerChannel << std::endl;
        }
        clusters.push_back(currentCluster);
        currentCluster = Cluster();
        hitCount = 0;
        highestTot = -1.0F;
        centerChannel = -1;
        earliestTime = std::numeric_limits<float>::max();
        earliestChannel = -1;
        earliestTot = 0.0F;
    };

    for (const auto& hit : sortedHits) {
        if (hit.timeOverThreshold <= 0.0F || hit.leadingTime < 0.0F) {
            continue;
        }

        const bool channelClose = hitCount > 0 && std::abs(hit.channel - lastChannel) <= 2;
        const bool timeClose = hitCount > 0 && std::fabs(hit.leadingTime - lastTime) <= 15.0;

        if (hitCount == 0 || (!channelClose || !timeClose)) {
            flushCluster();
            currentCluster.channels.push_back(hit.channel);
            hitCount = 1;
            highestTot = hit.timeOverThreshold;
            centerChannel = hit.channel;
            earliestTime = hit.leadingTime;
            earliestChannel = hit.channel;
            earliestTot = hit.timeOverThreshold;
        } else {
            currentCluster.channels.push_back(hit.channel);
            ++hitCount;
            if (hit.timeOverThreshold > highestTot) {
                highestTot = hit.timeOverThreshold;
                centerChannel = hit.channel;
            }
            if (hit.leadingTime < earliestTime) {
                earliestTime = hit.leadingTime;
                earliestChannel = hit.channel;
                earliestTot = hit.timeOverThreshold;
            }
        }

        lastChannel = hit.channel;
        lastTime = hit.leadingTime;
    }

    flushCluster();
    return clusters;
}

}  // namespace cluster_utils
