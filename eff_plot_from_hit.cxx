#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1I.h>
#include <TH1F.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <regex>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include "cluster.h"
#include "layer.h"
#include "cluster_utils.h"

using cluster_utils::LayerClusters;
using cluster_utils::buildClusters;

// Entry point for efficiency studies: consumes hit data, builds layers, and prints clusters.

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: ./plot <root_file1> ..." << std::endl;
        return 1;
    }

    // Open input ROOT file containing the hit tree.
    TFile* file = TFile::Open(argv[1]);
    if (!file || file->IsZombie()) {
        std::cout << "Could not open file: " << argv[1] << std::endl;
        return 1;
    }

    // Grab the tree that stores per-event hit information.
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cout << "Tree 'tree' not found in the file." << std::endl;
        file->Close();
        return 1;
    }
        
        int nevents = tree->GetEntries();
        std::cout << "Entries: " << nevents << std::endl;

        int nHits = 0;
        std::vector<int>* hit_channel = nullptr;
        std::vector<float>* hit_time1 = nullptr;
        std::vector<float>* hit_time2 = nullptr;
        std::vector<int>* hit_rise = nullptr;
        std::vector<int>* hit_bcid = nullptr;

        
        tree->SetBranchAddress("nHits", &nHits);
        tree->SetBranchAddress("hit_channel", &hit_channel);
        tree->SetBranchAddress("hit_time1", &hit_time1);
        tree->SetBranchAddress("hit_time2", &hit_time2);
        tree->SetBranchAddress("hit_rise", &hit_rise);
        tree->SetBranchAddress("hit_bcid", &hit_bcid);

        Long64_t nEntries = tree->GetEntries();
        // Analysis buffers per event: layers, resulting clusters, hit channels, and trigger time.
        std::vector<std::array<Layer, 3>> layersPerEvent;
        std::vector<std::array<LayerClusters, 3>> clustersPerEvent;
        std::vector<std::vector<int>> hitChannelsPerEvent;
        std::vector<float> triggerTimes;
        layersPerEvent.reserve(nEntries);
        clustersPerEvent.reserve(nEntries);
        hitChannelsPerEvent.reserve(nEntries);
        triggerTimes.reserve(nEntries);

        // Event loop: fill layers, compute clusters, capture trigger info.
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            Layer ly0, ly1, ly2;
            float triggerTime = -1.0F;
            std::vector<int> hitChannelsThisEvent;

            const auto hitsInEvent = hit_channel ? hit_channel->size() : 0;
            // Helper converting raw TDC edges into global time using BCID correction.
            auto correctedTime = [&](size_t idx, bool useTime1) -> float {
                float rawTime = useTime1 ? hit_time1->at(idx) : hit_time2->at(idx);
                int bcid = hit_bcid ? hit_bcid->at(idx) : 0;
                return static_cast<float>(bcid % 256) * 30.0F + (rawTime - 1.0F);
            };
            for (size_t i = 0; i < hitsInEvent; ++i) {
                int currentChannel = hit_channel->at(i);
                bool belongsToEta1 = hit_time1->at(i) > 0 && hit_time2->at(i) == 0;
                bool belongsToEta2 = hit_time1->at(i) == 0 && hit_time2->at(i) > 0;
                if (!belongsToEta1 && !belongsToEta2) {
                    continue;  // hit must be clearly assigned to exactly one eta
                }

                bool isLeadingEdge = hit_rise && hit_rise->at(i) == 1;
                float timeValue = correctedTime(i, belongsToEta1);
                if (currentChannel == 143) {
                    if (isLeadingEdge && triggerTime < 0.0F) {
                        triggerTime = timeValue;  // trigger defined by ch. 143 rising edge
                    }
                    continue;  // do not cluster the trigger channel
                }

                Layer* currentLayer = nullptr;
                if (currentChannel < 24) {
                    currentLayer = &ly0;
                } else if (currentChannel >= 24 && currentChannel < 48) {
                    currentLayer = &ly1;
                } else if (currentChannel >= 48 && currentChannel < 72) {
                    currentLayer = &ly2;
                } else {
                    continue;
                }

                currentLayer->addHit(currentChannel % 24 + 1, belongsToEta1, isLeadingEdge, timeValue);
                hitChannelsThisEvent.push_back(currentChannel);
            } // End loop on hits

            // Persist per-event context for later analysis/printing.
            layersPerEvent.push_back({ly0, ly1, ly2});
            hitChannelsPerEvent.push_back(std::move(hitChannelsThisEvent));
            triggerTimes.push_back(triggerTime);

            auto buildEtaClusters = [&](const std::vector<Layer::Hit>& hits) {
                std::vector<Layer::Hit> validHits;
                validHits.reserve(hits.size());
                for (const auto& hit : hits) {
                    if (hit.leadingTime >= 0.0F && hit.trailingTime >= 0.0F) {
                        validHits.push_back(hit);
                    }
                }
                return buildClusters(validHits);
            };

            std::array<LayerClusters, 3> eventClusters;
            eventClusters[0].eta1 = buildEtaClusters(ly0.eta1Hits);
            eventClusters[0].eta2 = buildEtaClusters(ly0.eta2Hits);
            eventClusters[1].eta1 = buildEtaClusters(ly1.eta1Hits);
            eventClusters[1].eta2 = buildEtaClusters(ly1.eta2Hits);
            eventClusters[2].eta1 = buildEtaClusters(ly2.eta1Hits);
            eventClusters[2].eta2 = buildEtaClusters(ly2.eta2Hits);
            clustersPerEvent.push_back(std::move(eventClusters));
        } // End loop over entries

        std::cout << "Stored layer info for " << layersPerEvent.size() << " events." << std::endl;
        std::cout << "Built clusters for " << clustersPerEvent.size() << " events (separated per layer)." << std::endl;
        std::cout << "Stored trigger times for " << triggerTimes.size() << " events." << std::endl;

        // Compute efficiency per layer/eta: (triggered events with ≥1 cluster) / (events with trigger).
        int triggeredEvents = 0;
        std::array<std::array<int, 2>, 3> eventsWithClusters{};
        for (size_t evt = 0; evt < triggerTimes.size(); ++evt) {
            if (triggerTimes[evt] < 0.0F) {
                continue;
            }
            ++triggeredEvents;
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                if (!clustersPerEvent[evt][layerIdx].eta1.empty()) {
                    ++eventsWithClusters[layerIdx][0];
                }
                if (!clustersPerEvent[evt][layerIdx].eta2.empty()) {
                    ++eventsWithClusters[layerIdx][1];
                }
            }
        }
        if (triggeredEvents > 0) {
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                for (size_t side = 0; side < 2; ++side) {
                    double eff = static_cast<double>(eventsWithClusters[layerIdx][side]) / triggeredEvents;
                    std::cout << "Layer " << layerIdx << " eta" << (side + 1) << " efficiency: "
                              << eff << " (" << eventsWithClusters[layerIdx][side] << "/" << triggeredEvents << ")"
                              << std::endl;
                }
            }
        } else {
            std::cout << "No triggered events found; efficiencies unavailable." << std::endl;
        }

        // Report per-event hit channels and resulting clusters for debugging/validation.
        for (size_t evt = 0; evt < clustersPerEvent.size(); ++evt) {
            std::cout << "Event " << evt << ":\n";
            std::cout << "  Hit channels: ";
            const auto& eventChannels = hitChannelsPerEvent[evt];
            for (size_t i = 0; i < eventChannels.size(); ++i) {
                std::cout << eventChannels[i];
                if (i + 1 != eventChannels.size()) {
                        std::cout << ", ";
                }
            }
            std::cout << std::endl;
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                const auto& layerClusters = clustersPerEvent[evt][layerIdx];
                // Print helper to display center channel and multiplicity per eta.
                auto printClusters = [&](const std::vector<Cluster>& clusters, const char* etaLabel) {
                    std::cout << "  Layer " << layerIdx << " " << etaLabel << " clusters: " << clusters.size() << std::endl;
                    for (size_t i = 0; i < clusters.size(); ++i) {
                        const auto& cluster = clusters[i];
                        std::cout << "    Cluster " << i << " center=" << cluster.centerChannel
                                  << " channels=" << cluster.channels.size() << std::endl;
                    }
                };
                printClusters(layerClusters.eta1, "eta1");
                printClusters(layerClusters.eta2, "eta2");
            }
        }
        
        file->Close();
        
    return 0;
}
