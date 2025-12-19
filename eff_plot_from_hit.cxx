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
        std::vector<std::vector<float>> hitTimesPerEvent;
        std::vector<float> triggerTimes;
        size_t correctTriggerCount = 0;
        layersPerEvent.reserve(nEntries);
        clustersPerEvent.reserve(nEntries);
        hitChannelsPerEvent.reserve(nEntries);
        hitTimesPerEvent.reserve(nEntries);
        triggerTimes.reserve(nEntries);

        // Event loop: fill layers, compute clusters, capture trigger info.
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            Layer ly0, ly1, ly2;
            float triggerTime = -1.0F;
            std::vector<int> hitChannelsThisEvent;
            std::vector<float> hitTimesThisEvent;
            bool hasValidTrigger = false;
            struct PendingHit {
                int layerIndex = -1;
                int relativeChannel = -1;
                int originalChannel = -1;
                bool isEta1 = false;
                bool isLeadingEdge = false;
                float time = 0.0F;
            };
            std::vector<PendingHit> pendingHits;
            bool skipEvent = false;

            const auto hitsInEvent = hit_channel ? hit_channel->size() : 0;
            // Helper converting raw TDC edges into global time using BCID correction.
            auto correctedTime = [&](size_t idx, bool useTime1) -> float {
                float rawTime = useTime1 ? hit_time1->at(idx) : hit_time2->at(idx);
                int bcid = hit_bcid ? hit_bcid->at(idx) : 0;
                return static_cast<float>(bcid % 256) * 30.0F + (rawTime - 1.0F);
            };
            for (size_t i = 0; i < hitsInEvent; ++i) {
                int currentChannel = hit_channel->at(i);
                const float rawTime1 = hit_time1->at(i);
                const float rawTime2 = hit_time2->at(i);
                bool hasEta1Time = rawTime1 > 0.0F;
                bool hasEta2Time = rawTime2 > 0.0F;
                //std::cout << "Processing hit " << i << ": channel=" << currentChannel
                //          << "  Total hits in event: " << hitsInEvent
                //          << " time1=" << rawTime1
                //          << " time2=" << rawTime2
                //          << " hasEta1Time=" << hasEta1Time
                //          << " hasEta2Time=" << hasEta2Time
                //          << std::endl;
                if (!hasEta1Time && !hasEta2Time) {
                    continue;  // no usable timing information
                }

                bool isLeadingEdge = hit_rise && hit_rise->at(i) == 1;

                if (currentChannel == 143) {
                    // Trigger is defined only on eta2 rising edge of channel 143.
                    if (isLeadingEdge && hasEta2Time && triggerTime < 0.0F) {
                        float triggerCandidateTime = correctedTime(i, /*useTime1=*/false);
                        if (triggerCandidateTime < 120.0F) {
                            skipEvent = true;  // ignore events with early trigger time
                            break;
                        }
                        triggerTime = triggerCandidateTime;  // trigger defined by ch. 143 rising edge
                        hasValidTrigger = true;
                    }
                    continue;  // do not cluster the trigger channel
                }
                
                int layerIndex = (currentChannel % 24) / 8;
                
                //std::cout << "For entry " << entry << ", Hit channel: " << currentChannel << ", Layer: " << layerIndex << std::endl;
                int connector = currentChannel / 24;
                int strip = 8 * connector + currentChannel % 8;
                //std::cout << "entry " << entry << ": Hit channel " << currentChannel
                //          << " mapped to layer " << layerIndex
                //          << ", connector " << connector
                //          << ", strip " << strip 
                //          << ", eta " << (belongsToEta1 ? "1" : "2")
                //          << ", edge " << (isLeadingEdge ? "Leading" : "Trailing")
                //          << ", time " << timeValue
                //          << std::endl;
                if (hasEta1Time) {
                    pendingHits.push_back(PendingHit{
                        layerIndex,
                        strip,
                        currentChannel,
                        /*isEta1=*/true,
                        isLeadingEdge,
                        correctedTime(i, /*useTime1=*/true)});
                }
                if (hasEta2Time) {
                    pendingHits.push_back(PendingHit{
                        layerIndex,
                        strip,
                        currentChannel,
                        /*isEta1=*/false,
                        isLeadingEdge,
                        correctedTime(i, /*useTime1=*/false)});
                }
            } // End loop on hits
            if (skipEvent) {
                continue;  // discard event when trigger time is below threshold
            }
            if (hasValidTrigger) {
                ++correctTriggerCount;  // track how many times a valid trigger is found
            }
            //std::cout << "Event " << entry << ": found " << pendingHits.size() << " valid hits over a total of hits available " << hit_channel->size() << " channels, trigger time = " << triggerTime << std::endl;
            //std::cout << "  Hit details (layer, channel, eta, edge, time):" << std::endl;
            //for (const auto& hit : pendingHits) {
              //  std::cout << "    Layer " << hit.layerIndex
                //          << ", Ch " << hit.originalChannel
                  //        << ", Eta " << (hit.isEta1 ? "1" : "2")
                    //      << ", " << (hit.isLeadingEdge ? "Leading" : "Trailing")
                      //    << ", Time " << hit.time << std::endl;
            //}
            // Keep only hits close in time to the trigger before forming clusters.
            for (const auto& hit : pendingHits) {
                if (hit.isLeadingEdge && triggerTime >= 0.0F) {
                    //std::cout << "Processing hit at time " << hit.time << " with trigger time " << triggerTime << " delta " << (hit.time - triggerTime) << std::endl;
                    float delta = hit.time - triggerTime;
                    //if (delta < -200.0F || delta > -130.0F) {
                    //    continue;
                    //}
                }

                Layer* currentLayer = nullptr;
                if (hit.layerIndex == 0) {
                    currentLayer = &ly0;
                } else if (hit.layerIndex == 1) {
                    currentLayer = &ly1;
                } else if (hit.layerIndex == 2) {
                    currentLayer = &ly2;
                }
                if (!currentLayer) {
                    continue;
                }

                currentLayer->addHit(hit.relativeChannel, hit.isEta1, hit.isLeadingEdge, hit.time);
                hitChannelsThisEvent.push_back(hit.originalChannel);
                hitTimesThisEvent.push_back(hit.time);
            }

            // Persist per-event context for later analysis/printing.
            layersPerEvent.push_back({ly0, ly1, ly2});
            hitChannelsPerEvent.push_back(std::move(hitChannelsThisEvent));
            hitTimesPerEvent.push_back(std::move(hitTimesThisEvent));
            triggerTimes.push_back(triggerTime);

            auto buildEtaClusters = [&](const std::vector<Layer::Hit>& hits, int layerIndex, int etaSide) {
                std::vector<Layer::Hit> validHits;
                validHits.reserve(hits.size());
                for (const auto& hit : hits) {
                    if (hit.leadingTime >= 0.0F &&
                        (hit.trailingTime >= 0.0F || hit.timeOverThreshold == 9999.0F)) {
                        validHits.push_back(hit);
                    }
                }
                return buildClusters(validHits, layerIndex, etaSide);
            };

            std::array<LayerClusters, 3> eventClusters;
            eventClusters[0].eta1 = buildEtaClusters(ly0.eta1Hits, 0, 1);
            eventClusters[0].eta2 = buildEtaClusters(ly0.eta2Hits, 0, 2);
            eventClusters[1].eta1 = buildEtaClusters(ly1.eta1Hits, 1, 1);
            eventClusters[1].eta2 = buildEtaClusters(ly1.eta2Hits, 1, 2);
            eventClusters[2].eta1 = buildEtaClusters(ly2.eta1Hits, 2, 1);
            eventClusters[2].eta2 = buildEtaClusters(ly2.eta2Hits, 2, 2);
            clustersPerEvent.push_back(std::move(eventClusters));
        } // End loop over entries

        std::cout << "Stored layer info for " << layersPerEvent.size() << " events." << std::endl;
        std::cout << "Built clusters for " << clustersPerEvent.size() << " events (separated per layer)." << std::endl;
        std::cout << "Stored trigger times for " << triggerTimes.size() << " events." << std::endl;
        std::cout << "Valid triggers found: " << correctTriggerCount << std::endl;

        // Hit presence efficiency per layer/side: events with hits over triggered events.
        std::array<size_t, 3> eventsWithEta1Hits{};
        std::array<size_t, 3> eventsWithEta2Hits{};
        for (size_t evt = 0; evt < triggerTimes.size(); ++evt) {
            if (triggerTimes[evt] < 0.0F) {
                continue;
            }
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                const auto& ly = layersPerEvent[evt][layerIdx];
                if (!ly.eta1Hits.empty()) {
                    ++eventsWithEta1Hits[layerIdx];
                }
                if (!ly.eta2Hits.empty()) {
                    ++eventsWithEta2Hits[layerIdx];
                }
            }
        }
        if (correctTriggerCount > 0) {
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                double effEta1Hits = static_cast<double>(eventsWithEta1Hits[layerIdx]) /
                                     static_cast<double>(correctTriggerCount);
                double effEta2Hits = static_cast<double>(eventsWithEta2Hits[layerIdx]) /
                                     static_cast<double>(correctTriggerCount);
                std::cout << "Layer " << layerIdx << " hit efficiency: "
                          << "eta1=" << effEta1Hits << " (" << eventsWithEta1Hits[layerIdx]
                          << "/" << correctTriggerCount << "), "
                          << "eta2=" << effEta2Hits << " (" << eventsWithEta2Hits[layerIdx]
                          << "/" << correctTriggerCount << ")"
                          << std::endl;
            }
        } else {
            std::cout << "No valid triggers found; hit efficiencies unavailable." << std::endl;
        }

        // Compute efficiency per layer: events with ≥1 eta1 cluster matched to an eta2 cluster (|Δchannel| ≤ 1)
        // divided by the total number of processed events.
        const size_t totalEvents = clustersPerEvent.size();
        std::array<size_t, 3> eventsWithMatchedClusters{};
        std::array<size_t, 3> eventsWithEta1Clusters{};
        std::array<size_t, 3> eventsWithEta2Clusters{};
        for (size_t evt = 0; evt < totalEvents; ++evt) {
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                const auto& eta1Clusters = clustersPerEvent[evt][layerIdx].eta1;
                const auto& eta2Clusters = clustersPerEvent[evt][layerIdx].eta2;
                if (!eta1Clusters.empty()) {
                    ++eventsWithEta1Clusters[layerIdx];
                }
                if (!eta2Clusters.empty()) {
                    ++eventsWithEta2Clusters[layerIdx];
                }
                bool matched = false;
                for (const auto& c1 : eta1Clusters) {
                    for (const auto& c2 : eta2Clusters) {
                        if (std::abs(c1.centerChannel - c2.centerChannel) <= 2) {
                            matched = true;
                            break;
                        }
                    }
                    if (matched) {
                        break;
                    }
                }
                if (matched) {
                    ++eventsWithMatchedClusters[layerIdx];
                }
            }
        }
        if (totalEvents > 0) {
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                double eff = static_cast<double>(eventsWithMatchedClusters[layerIdx]) /
                             static_cast<double>(totalEvents);
                std::cout << "Layer " << layerIdx << " matched efficiency: "
                          << eff << " (" << eventsWithMatchedClusters[layerIdx] << "/" << totalEvents << ")"
                          << std::endl;
                double effEta1 = static_cast<double>(eventsWithEta1Clusters[layerIdx]) /
                                 static_cast<double>(totalEvents);
                double effEta2 = static_cast<double>(eventsWithEta2Clusters[layerIdx]) /
                                 static_cast<double>(totalEvents);
                std::cout << "  Eta1 efficiency: " << effEta1 << " (" << eventsWithEta1Clusters[layerIdx]
                          << "/" << totalEvents << ")" << std::endl;
                std::cout << "  Eta2 efficiency: " << effEta2 << " (" << eventsWithEta2Clusters[layerIdx]
                          << "/" << totalEvents << ")" << std::endl;
            }
        } else {
            std::cout << "No events processed; efficiencies unavailable." << std::endl;
        }

        // Report per-event hit channels and resulting clusters for debugging/validation.
        const bool debug = true;
        if(debug){
        for (size_t evt = 0; evt < clustersPerEvent.size(); ++evt) {
            std::cout << "Event " << evt << ":\n";
            std::cout << "  Hit channels (time): ";
            const auto& eventChannels = hitChannelsPerEvent[evt];
            const auto& eventTimes = hitTimesPerEvent[evt];
            for (size_t i = 0; i < eventChannels.size(); ++i) {
                std::cout << eventChannels[i] << "(" << eventTimes[i] << ")";
                if (i + 1 != eventChannels.size()) {
                        std::cout << ", ";
                }
            }
            std::cout << std::endl;
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                const auto& layerClusters = clustersPerEvent[evt][layerIdx];
                const auto& layerHits = layersPerEvent[evt][layerIdx];
                std::cout << "  Layer " << layerIdx << " hits: eta1=" << layerHits.eta1Hits.size()
                          << " eta2=" << layerHits.eta2Hits.size() << std::endl;
                size_t totalHitsLayer = layerHits.eta1Hits.size() + layerHits.eta2Hits.size();
                size_t totalClustersLayer = layerClusters.eta1.size() + layerClusters.eta2.size();
                std::cout << "    Total hits (eta1+eta2): " << totalHitsLayer << std::endl;
                std::cout << "    Total clusters (eta1+eta2): " << totalClustersLayer << std::endl;
                // Print helper to display center channel and multiplicity per eta.
                auto printClusters = [&](const std::vector<Cluster>& clusters, const char* etaLabel) {
                    std::cout << "  Layer " << layerIdx << " " << etaLabel << " clusters: " << clusters.size() << std::endl;
                    for (size_t i = 0; i < clusters.size(); ++i) {
                        const auto& cluster = clusters[i];
                        std::cout << "    Cluster " << i << " center=" << cluster.centerChannel
                                  << " size=" << cluster.channels.size()
                                  << " time=" << cluster.time << std::endl;
                    }
                };
                printClusters(layerClusters.eta1, "eta1");
                printClusters(layerClusters.eta2, "eta2");
            }
        }
    }
        file->Close();
        
    return 0;
}
