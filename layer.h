#ifndef LAYER_H
#define LAYER_H

#include <vector>
#include <algorithm>

// Layer stores all hits for eta1/eta2 sides of a detector layer.
class Layer {
public:
    // Hit keeps the first leading/trailing edges and the derived TOT.
    struct Hit {
        int channel = 0;
        float leadingTime = -1.0F;
        float trailingTime = -1.0F;
        float timeOverThreshold = 0.0F;
    };

    std::vector<Hit> eta1Hits;  // Hits belonging to eta1 side
    std::vector<Hit> eta2Hits;  // Hits belonging to eta2 side

    // Adds a leading or trailing edge to the channel, computing TOT when possible.
    void addHit(int channel, bool isEta1, bool isLeadingEdge, float time) {
        std::vector<Hit>& hits = isEta1 ? eta1Hits : eta2Hits;
        auto it = std::find_if(hits.begin(), hits.end(), [channel](const Hit& h) {
            return h.channel == channel;
        });
        if (it == hits.end()) {
            hits.push_back(Hit{channel});
            it = hits.end() - 1;
        }

        if (isLeadingEdge) {
            if (it->leadingTime >= 0.0F) {
                return;  // first leading edge already stored
            }
            it->leadingTime = time;
        } else {
            it->trailingTime = time;
        }

        if (it->leadingTime >= 0.0F && it->trailingTime >= 0.0F) {
            float tot = it->trailingTime - it->leadingTime;
            if (tot != 9999.0F) {
                constexpr float period = 256.0F * 30.0F;
                constexpr float halfPeriod = period / 2.0F;
                if (tot < -halfPeriod) {
                    tot += period;
                }
                if (tot > halfPeriod) {
                    tot -= period;
                }
                // If trailing arrives before leading (even after wrap correction), discard this
                // trailing edge and keep the hit as “missing trailing” with TOT sentinel.
                if (tot <= 0.0F) {
                    it->trailingTime = -1.0F;
                    it->timeOverThreshold = 9999.0F;
                    return;
                }
            }
            it->timeOverThreshold = tot;
        } else if (it->leadingTime >= 0.0F && it->trailingTime < 0.0F) {
            // Missing trailing edge: keep the hit and mark TOT with sentinel.
            it->timeOverThreshold = 9999.0F;
        }
    }
};

#endif  // LAYER_H
