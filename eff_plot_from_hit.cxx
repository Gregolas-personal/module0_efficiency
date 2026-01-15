// ROOT framework includes for GUI and data structures
#include <TApplication.h>  // ROOT application for event loop
#include <TCanvas.h>        // Canvas for drawing plots
#include <TFile.h>          // ROOT file I/O
#include <TH2D.h>           // 2D histogram
#include <TROOT.h>          // Global ROOT object
#include <TString.h>        // ROOT string utilities
#include <TTree.h>          // ROOT tree for structured data

// Standard C++ library includes
#include <memory>     // Smart pointers (unique_ptr)
#include <string>     // String handling
#include <vector>     // Dynamic arrays
#include <regex>      // Regular expressions (not used in current code)
#include <array>      // Fixed-size arrays
#include <algorithm>  // Algorithms like std::abs
#include <map>        // Associative containers
#include <cmath>      // Mathematical functions
#include <limits>     // Numeric limits (not used in current code)

// Custom headers for detector data structures
#include "cluster.h"         // Cluster class definition
#include "layer.h"           // Layer class definition
#include "cluster_utils.h"   // Cluster building utilities
#include "event_display.h"   // Event visualization functions

// Import cluster utility types into current namespace
// Import cluster utility types into current namespace
using cluster_utils::LayerClusters;
using cluster_utils::buildClusters;

// Forward declarations of utility print functions (defined elsewhere)
void printUsage();  // Print command-line usage instructions
void printCouldNotOpenFile(const std::string& inputFile);  // Error message for file opening
void printTreeNotFound();  // Error when TTree is not found in file
void printEntries(int nevents);  // Print number of events in dataset
void printStoredInfo(size_t layersCount, size_t clustersCount, size_t triggersCount, size_t correctTriggerCount);  // Summary of collected data
void printLayerHitEfficiency(size_t layerIdx,
                             double effEta1, size_t numEta1,
                             double effEta2, size_t numEta2,
                             size_t denom);  // Print hit efficiency per layer and eta side
void printNoValidTriggers();  // Print warning when no valid triggers found
void printLayerMatchedEfficiency(size_t layerIdx, double eff, size_t numerator, size_t denom);  // Print matched cluster efficiency
void printEtaEfficiency(const char* label, double eff, size_t numerator, size_t denom);  // Print efficiency for a specific eta side
void printNoEventsProcessed();  // Print warning when no events were processed
void printEventDebug(size_t evt,
                     const std::vector<int>& eventChannels,
                     const std::vector<float>& eventTimes,
                     const std::array<Layer, 3>& layers,
                     const std::array<LayerClusters, 3>& clusters);  // Debug output for each event

// Main entry point for efficiency analysis program
// Processes detector hit data, builds clusters, and calculates detection efficiencies
int main(int argc, char* argv[]) {
    // Check for minimum required command-line arguments
    if (argc < 2) {
        printUsage();  // Show usage if no input file specified
        return 1;      // Exit with error code
    }
    // Get input ROOT file path from command line
    std::string inputFile = argv[1];
    // By default, show interactive event display
    bool showDisplay = true;
    // Check for optional --skip-display flag
    if (argc >= 3 && std::string(argv[2]) == "--skip-display") {
        showDisplay = false;  // Disable GUI display
    }

    // Initialize ROOT application for interactive graphics (if enabled)
    std::unique_ptr<TApplication> appPtr;  // Smart pointer manages lifetime
    if (showDisplay) {
        // Create ROOT application with command-line args
        appPtr = std::make_unique<TApplication>("app", &argc, argv);
        gROOT->SetBatch(kFALSE);  // Enable interactive mode
    } else {
        gROOT->SetBatch(kTRUE);  // Batch mode (no GUI windows)
    }

    // Open the ROOT file containing detector hit data
    TFile* file = TFile::Open(inputFile.c_str());
    if (!file || file->IsZombie()) {  // Check if file opened successfully
        printCouldNotOpenFile(inputFile);
        return 1;  // Exit with error
    }

    // Retrieve the TTree named "tree" containing hit information
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {  // Check if tree exists
        printTreeNotFound();
        file->Close();  // Clean up file handle
        return 1;       // Exit with error
    }
        
        // Get total number of events in the tree
        int nevents = tree->GetEntries();
        printEntries(nevents);  // Print to console

        // Declare variables to store hit data from tree branches
        int nHits = 0;                      // Number of hits in current event
        std::vector<int>* hit_channel = nullptr;     // Channel numbers (detector strips)
        std::vector<float>* hit_time1 = nullptr;     // Eta1 timing information
        std::vector<float>* hit_time2 = nullptr;     // Eta2 timing information
        std::vector<int>* hit_rise = nullptr;        // Edge type (1=rising/leading, 0=falling/trailing)
        std::vector<int>* hit_bcid = nullptr;        // Bunch crossing ID for time correction

        // Link tree branches to local variables for reading data
        tree->SetBranchAddress("nHits", &nHits);           // Number of hits
        tree->SetBranchAddress("hit_channel", &hit_channel);  // Channel IDs
        tree->SetBranchAddress("hit_time1", &hit_time1);   // Eta1 times
        tree->SetBranchAddress("hit_time2", &hit_time2);   // Eta2 times
        tree->SetBranchAddress("hit_rise", &hit_rise);     // Rising/falling edge
        tree->SetBranchAddress("hit_bcid", &hit_bcid);     // Bunch crossing ID

        Long64_t nEntries = tree->GetEntries();  // Total entries in tree
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

        // Histogram to capture clusters-per-event vs cluster size (axes intentionally swapped).
        const int maxClusterSizeBins = 20;
        const int maxClustersPerEventBins = 20;
        TH2D* hClustersPerEventVsSize = new TH2D(
            "hClustersPerEventVsSize",
            "Clusters per event vs cluster size;Clusters per event;Cluster size (#channels)",
            maxClustersPerEventBins, -0.5, maxClustersPerEventBins - 0.5,
            maxClusterSizeBins, 0.5, maxClusterSizeBins + 0.5);

        // Event loop: fill layers, compute clusters, capture trigger info.
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            Layer ly0, ly1, ly2;
            float triggerTime = -1.0F;
            std::vector<int> hitChannelsThisEvent;
            std::vector<float> hitTimesThisEvent;
            bool hasValidTrigger = false;
            // Structure to hold pending hit information before filtering
            struct PendingHit {
                int layerIndex = -1;       // Which layer (0, 1, or 2)
                int relativeChannel = -1;  // Strip number within layer
                int originalChannel = -1;  // Original global channel number
                bool isEta1 = false;       // True if eta1, false if eta2
                bool isLeadingEdge = false;  // True for rising edge, false for falling
                float time = 0.0F;         // Corrected time value
            };
            std::vector<PendingHit> pendingHits;  // Temporary storage for all hits
            bool skipEvent = false;  // Flag to skip events with bad trigger timing

            const auto hitsInEvent = hit_channel ? hit_channel->size() : 0;  // Number of hits in this event
            // Lambda function: convert raw DCT time to global time using BCID correction
            auto correctedTime = [&](size_t idx, bool useTime1) -> float {
                // Select either eta1 (time1) or eta2 (time2) raw time
                float rawTime = useTime1 ? hit_time1->at(idx) : hit_time2->at(idx);
                int bcid = hit_bcid ? hit_bcid->at(idx) : 0;  // Bunch crossing ID
                // BCID correction: (bcid % 256) * 30ns + (rawTime - 1)
                return static_cast<float>(bcid % 256) * 30.0F + (rawTime - 1.0F);
            };
            // Loop over all hits in this event
            for (size_t i = 0; i < hitsInEvent; ++i) {
                int currentChannel = hit_channel->at(i);  // Global channel number
                const float rawTime1 = hit_time1->at(i);  // Eta1 raw time
                const float rawTime2 = hit_time2->at(i);  // Eta2 raw time
                bool hasEta1Time = rawTime1 > 0.0F;  // Check if eta1 has valid time
                bool hasEta2Time = rawTime2 > 0.0F;  // Check if eta2 has valid time
                //std::cout << "Processing hit " << i << ": channel=" << currentChannel
                //          << "  Total hits in event: " << hitsInEvent
                //          << " time1=" << rawTime1
                //          << " time2=" << rawTime2
                //          << " hasEta1Time=" << hasEta1Time
                //          << " hasEta2Time=" << hasEta2Time
                //          << std::endl;
                if (!hasEta1Time && !hasEta2Time) {
                    continue;  // Skip hits with no timing information
                }

                bool isLeadingEdge = hit_rise && hit_rise->at(i) == 1;  // 1 = rising edge

                // Special handling for trigger channel (143)
                if (currentChannel == 143) {
                    // Trigger defined only on eta2 rising edge of channel 143
                    if (isLeadingEdge && hasEta2Time && triggerTime < 0.0F) {
                        float triggerCandidateTime = correctedTime(i, /*useTime1=*/false);
                        if (triggerCandidateTime < 120.0F) {
                            skipEvent = true;  // Reject events with early trigger (< 120ns)
                            break;
                        }
                        triggerTime = triggerCandidateTime;  // Store trigger time
                        hasValidTrigger = true;
                    }
                    continue;  // Don't include trigger channel in clustering
                }
                
                // Map global channel to layer index (0, 1, or 2)
                // Formula: layerIndex = (channel % 24) / 8
                int layerIndex = (currentChannel % 24) / 8;
                
                //std::cout << "For entry " << entry << ", Hit channel: " << currentChannel << ", Layer: " << layerIndex << std::endl;
                // Calculate connector number (0 to 5) based on channel
                int connector = currentChannel / 24;
                // Calculate strip number within layer (0 to 39)
                int strip = 8 * connector + currentChannel % 8;
                //std::cout << "entry " << entry << ": Hit channel " << currentChannel
                //          << " mapped to layer " << layerIndex
                //          << ", connector " << connector
                //          << ", strip " << strip 
                //          << ", eta " << (belongsToEta1 ? "1" : "2")
                //          << ", edge " << (isLeadingEdge ? "Leading" : "Trailing")
                //          << ", time " << timeValue
                //          << std::endl;
                // Add eta1 hit if valid time exists
                if (hasEta1Time) {
                    pendingHits.push_back(PendingHit{
                        layerIndex,                        // Layer 0, 1, or 2
                        strip,                             // Strip number
                        currentChannel,                    // Original channel ID
                        /*isEta1=*/true,                   // Eta1 side
                        isLeadingEdge,                     // Rising/falling edge
                        correctedTime(i, /*useTime1=*/true)});  // Corrected time
                }
                // Add eta2 hit if valid time exists
                if (hasEta2Time) {
                    pendingHits.push_back(PendingHit{
                        layerIndex,                         // Layer 0, 1, or 2
                        strip,                              // Strip number
                        currentChannel,                     // Original channel ID
                        /*isEta1=*/false,                   // Eta2 side
                        isLeadingEdge,                      // Rising/falling edge
                        correctedTime(i, /*useTime1=*/false)});  // Corrected time
                }
            } // End loop over hits in this event
            if (skipEvent) {
                continue;  // Discard event if trigger time was too early
            }
            if (hasValidTrigger) {
                ++correctTriggerCount;  // Increment counter for valid triggers
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
            // Filter hits: keep only leading edges within time window relative to trigger
            for (const auto& hit : pendingHits) {
                // Only apply time filter to leading edges with valid trigger
                if (hit.isLeadingEdge && triggerTime >= 0.0F) {
                    //std::cout << "Processing hit at time " << hit.time << " with trigger time " << triggerTime << " delta " << (hit.time - triggerTime) << std::endl;
                    float delta = hit.time - triggerTime;  // Time difference from trigger
                    // Keep hits between -200ns and -130ns relative to trigger
                    if (delta < -200.0F || delta > -130.0F) {
                        continue;  // Skip hits outside time window
                    }
                }

                // Select the correct layer object based on layer index
                Layer* currentLayer = nullptr;
                if (hit.layerIndex == 0) {
                    currentLayer = &ly0;  // Layer 0
                } else if (hit.layerIndex == 1) {
                    currentLayer = &ly1;  // Layer 1
                } else if (hit.layerIndex == 2) {
                    currentLayer = &ly2;  // Layer 2
                }
                if (!currentLayer) {
                    continue;  // Skip if layer index is invalid
                }

                // Add this hit to the appropriate layer
                currentLayer->addHit(hit.relativeChannel, hit.isEta1, hit.isLeadingEdge, hit.time);
                // Store hit information for this event
                hitChannelsThisEvent.push_back(hit.originalChannel);
                hitTimesThisEvent.push_back(hit.time);
            }

            // Store per-event data in vectors for later analysis
            layersPerEvent.push_back({ly0, ly1, ly2});  // Save all 3 layers
            hitChannelsPerEvent.push_back(std::move(hitChannelsThisEvent));  // Save hit channels
            hitTimesPerEvent.push_back(std::move(hitTimesThisEvent));  // Save hit times
            triggerTimes.push_back(triggerTime);  // Save trigger time

            // Lambda function: build clusters from hits with valid timing
            auto buildEtaClusters = [&](const std::vector<Layer::Hit>& hits, int layerIndex, int etaSide) {
                std::vector<Layer::Hit> validHits;  // Filter for valid hits only
                validHits.reserve(hits.size());
                for (const auto& hit : hits) {
                    // Require valid leading time and either trailing time or TOT=9999
                    if (hit.leadingTime >= 0.0F &&
                        (hit.trailingTime >= 0.0F || hit.timeOverThreshold == 9999.0F)) {
                        validHits.push_back(hit);
                    }
                }
                // Call cluster building algorithm
                return buildClusters(validHits, layerIndex, etaSide);
            };

            // Build clusters for all layers and both eta sides
            std::array<LayerClusters, 3> eventClusters;
            eventClusters[0].eta1 = buildEtaClusters(ly0.eta1Hits, 0, 1);  // Layer 0, eta1
            eventClusters[0].eta2 = buildEtaClusters(ly0.eta2Hits, 0, 2);  // Layer 0, eta2
            eventClusters[1].eta1 = buildEtaClusters(ly1.eta1Hits, 1, 1);  // Layer 1, eta1
            eventClusters[1].eta2 = buildEtaClusters(ly1.eta2Hits, 1, 2);  // Layer 1, eta2
            eventClusters[2].eta1 = buildEtaClusters(ly2.eta1Hits, 2, 1);  // Layer 2, eta1
            eventClusters[2].eta2 = buildEtaClusters(ly2.eta2Hits, 2, 2);  // Layer 2, eta2
            clustersPerEvent.push_back(std::move(eventClusters));  // Store clusters for this event
        } // End loop over all events

        /********************************** END OF CLUSTER ALGORITHM **********************************/
        // Analysis and plotting section

        // Fill 2D histogram: count clusters by size for each event
        for (size_t evt = 0; evt < clustersPerEvent.size(); ++evt) {
            std::map<int, int> clustersPerSize;  // Map: cluster size -> count
            const auto& eventClusters = clustersPerEvent[evt];
            // Lambda to accumulate cluster sizes across all clusters
            auto accumulateClusterSizes = [&clustersPerSize](const std::vector<Cluster>& clusters) {
                for (const auto& cluster : clusters) {
                    int size = static_cast<int>(cluster.channels.size());  // Number of channels in cluster
                    ++clustersPerSize[size];  // Increment count for this size
                }
            };
            // Process clusters from all layers and both eta sides
            for (const auto& layerClusters : eventClusters) {
                accumulateClusterSizes(layerClusters.eta1);  // Eta1 clusters
                accumulateClusterSizes(layerClusters.eta2);  // Eta2 clusters
            }
            // Fill histogram for each cluster size found in this event
            for (const auto& [clusterSize, count] : clustersPerSize) {
                hClustersPerEventVsSize->Fill(count, clusterSize);  // X=count, Y=size
            }
        }


        // Print summary of collected data
        printStoredInfo(layersPerEvent.size(), clustersPerEvent.size(), triggerTimes.size(), correctTriggerCount);

        // === HIT PRESENCE EFFICIENCY CALCULATION ===
        // Count events with hits in each layer/eta side (only for triggered events)
        std::array<size_t, 3> eventsWithEta1Hits{};  // Initialize to zero
        std::array<size_t, 3> eventsWithEta2Hits{};  // Initialize to zero
        for (size_t evt = 0; evt < triggerTimes.size(); ++evt) {
            if (triggerTimes[evt] < 0.0F) {
                continue;  // Skip events without valid trigger
            }
            // Check each of the 3 layers
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                const auto& ly = layersPerEvent[evt][layerIdx];
                if (!ly.eta1Hits.empty()) {
                    ++eventsWithEta1Hits[layerIdx];  // Count event if eta1 has hits
                }
                if (!ly.eta2Hits.empty()) {
                    ++eventsWithEta2Hits[layerIdx];  // Count event if eta2 has hits
                }
            }
        }
        // Calculate and print hit efficiency per layer
        if (correctTriggerCount > 0) {
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                // Efficiency = (events with hits) / (total triggered events)
                double effEta1Hits = static_cast<double>(eventsWithEta1Hits[layerIdx]) /
                                     static_cast<double>(correctTriggerCount);
                double effEta2Hits = static_cast<double>(eventsWithEta2Hits[layerIdx]) /
                                     static_cast<double>(correctTriggerCount);
                printLayerHitEfficiency(layerIdx,
                                        effEta1Hits, eventsWithEta1Hits[layerIdx],
                                        effEta2Hits, eventsWithEta2Hits[layerIdx],
                                        correctTriggerCount);
            }
        } else {
            printNoValidTriggers();  // Warning if no triggers found
        }

        // === MATCHED CLUSTER EFFICIENCY CALCULATION ===
        // Count events with matched clusters (|Δchannel| ≤ 2 between eta1 and eta2)
        const size_t totalEvents = clustersPerEvent.size();
        std::array<size_t, 3> eventsWithMatchedClusters{};  // Events with matched clusters
        std::array<size_t, 3> eventsWithEta1Clusters{};     // Events with any eta1 cluster
        std::array<size_t, 3> eventsWithEta2Clusters{};     // Events with any eta2 cluster
        for (size_t evt = 0; evt < totalEvents; ++evt) {
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                const auto& eta1Clusters = clustersPerEvent[evt][layerIdx].eta1;
                const auto& eta2Clusters = clustersPerEvent[evt][layerIdx].eta2;
                if (!eta1Clusters.empty()) {
                    ++eventsWithEta1Clusters[layerIdx];  // Count event
                }
                if (!eta2Clusters.empty()) {
                    ++eventsWithEta2Clusters[layerIdx];  // Count event
                }
                // Check for matched cluster pairs
                bool matched = false;
                for (const auto& c1 : eta1Clusters) {
                    for (const auto& c2 : eta2Clusters) {
                        // Clusters match if center channels are within 2 of each other
                        if (std::abs(c1.centerChannel - c2.centerChannel) <= 2) {
                            matched = true;
                            break;  // Found a match, stop searching
                        }
                    }
                    if (matched) {
                        break;  // Found a match, stop searching
                    }
                }
                if (matched) {
                    ++eventsWithMatchedClusters[layerIdx];  // Count event with match
                }
            }
        }
        
        // Run interactive event display (if enabled)
        const double layerSpacingCm = 3.0;  // 3 cm between layers
        runPerEventDisplay(clustersPerEvent, showDisplay, appPtr.get(), layerSpacingCm);
        
        // Calculate and print efficiencies
        if (totalEvents > 0) {
            for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
                // Matched cluster efficiency = (events with matched clusters) / (total events)
                double eff = static_cast<double>(eventsWithMatchedClusters[layerIdx]) /
                             static_cast<double>(totalEvents);
                printLayerMatchedEfficiency(layerIdx, eff, eventsWithMatchedClusters[layerIdx], totalEvents);
                // Eta1 cluster efficiency = (events with eta1 clusters) / (total events)
                double effEta1 = static_cast<double>(eventsWithEta1Clusters[layerIdx]) /
                                 static_cast<double>(totalEvents);
                // Eta2 cluster efficiency = (events with eta2 clusters) / (total events)
                double effEta2 = static_cast<double>(eventsWithEta2Clusters[layerIdx]) /
                                 static_cast<double>(totalEvents);
                printEtaEfficiency("Eta1", effEta1, eventsWithEta1Clusters[layerIdx], totalEvents);
                printEtaEfficiency("Eta2", effEta2, eventsWithEta2Clusters[layerIdx], totalEvents);
            }
        } else {
            printNoEventsProcessed();  // Warning if no events were processed
        }

        // === DEBUG OUTPUT ===
        // Print detailed information for each event (channels, times, clusters)
        const bool debug = true;  // Enable debug output
        if(debug){
        for (size_t evt = 0; evt < clustersPerEvent.size(); ++evt) {
            printEventDebug(evt,
                            hitChannelsPerEvent[evt],  // Hit channels in this event
                            hitTimesPerEvent[evt],     // Hit times in this event
                            layersPerEvent[evt],       // Layer data for this event
                            clustersPerEvent[evt]);    // Cluster data for this event
        }
    }

        // === SAVE HISTOGRAM ===
        // Create canvas and draw clusters-per-event vs cluster-size histogram
        TCanvas* cClustersVsSize = new TCanvas("cClustersVsSize", "Clusters per event vs cluster size", 800, 600);
        hClustersPerEventVsSize->SetStats(false);  // Don't show statistics box
        hClustersPerEventVsSize->Draw("COLZ");     // Draw as color plot (Z-axis shows counts)
        cClustersVsSize->SaveAs("clusters_per_event_vs_size.png");  // Save to file
        file->Close();  // Close the ROOT input file
        
    return 0;  // Exit successfully
}  // End of main function
