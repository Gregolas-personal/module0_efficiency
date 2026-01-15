#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "cluster_utils.h"
#include "layer.h"

void printUsage() {
    std::cout << "Usage: ./plot <root_file1> [--skip-display]" << std::endl;
}

void printCouldNotOpenFile(const std::string& inputFile) {
    std::cout << "Could not open file: " << inputFile << std::endl;
}

void printTreeNotFound() {
    std::cout << "Tree 'tree' not found in the file." << std::endl;
}

void printEntries(int nevents) {
    std::cout << "Entries: " << nevents << std::endl;
}

void printStoredInfo(size_t layersCount, size_t clustersCount, size_t triggersCount, size_t correctTriggerCount) {
    std::cout << "Stored layer info for " << layersCount << " events." << std::endl;
    std::cout << "Built clusters for " << clustersCount << " events (separated per layer)." << std::endl;
    std::cout << "Stored trigger times for " << triggersCount << " events." << std::endl;
    std::cout << "Valid triggers found: " << correctTriggerCount << std::endl;
}

void printLayerHitEfficiency(size_t layerIdx,
                             double effEta1, size_t numEta1,
                             double effEta2, size_t numEta2,
                             size_t denom) {
    std::cout << "Layer " << layerIdx << " hit efficiency: "
              << "eta1=" << effEta1 << " (" << numEta1 << "/" << denom << "), "
              << "eta2=" << effEta2 << " (" << numEta2 << "/" << denom << ")"
              << std::endl;
}

void printNoValidTriggers() {
    std::cout << "No valid triggers found; hit efficiencies unavailable." << std::endl;
}

void printLayerMatchedEfficiency(size_t layerIdx, double eff, size_t numerator, size_t denom) {
    std::cout << "Layer " << layerIdx << " matched efficiency: "
              << eff << " (" << numerator << "/" << denom << ")" << std::endl;
}

void printEtaEfficiency(const char* label, double eff, size_t numerator, size_t denom) {
    std::cout << "  " << label << " efficiency: " << eff
              << " (" << numerator << "/" << denom << ")" << std::endl;
}

void printNoEventsProcessed() {
    std::cout << "No events processed; efficiencies unavailable." << std::endl;
}

void printEventDisplayPrompt(size_t evt) {
    std::cout << "Displayed event " << evt << ". Press Enter to continue." << std::endl;
}

void print3DLineFit(size_t evt,
                    double etaIntercept, double etaSlope, double etaChi2,
                    double dtIntercept, double dtSlope, double dtChi2) {
    std::cout << "Event " << evt << " 3D line fit (eta vs z): y="
              << etaIntercept << "+" << etaSlope << "x, chi2=" << etaChi2 << std::endl;
    std::cout << "Event " << evt << " 3D line fit (dt vs z): y="
              << dtIntercept << "+" << dtSlope << "x, chi2=" << dtChi2 << std::endl;
}

static void printClusters(size_t layerIdx, const std::vector<Cluster>& clusters, const char* etaLabel) {
    std::cout << "  Layer " << layerIdx << " " << etaLabel << " clusters: " << clusters.size() << std::endl;
    for (size_t i = 0; i < clusters.size(); ++i) {
        const auto& cluster = clusters[i];
        std::cout << "    Cluster " << i << " center=" << cluster.centerChannel
                  << " size=" << cluster.channels.size()
                  << " time=" << cluster.time << std::endl;
    }
}

void printEventDebug(size_t evt,
                     const std::vector<int>& eventChannels,
                     const std::vector<float>& eventTimes,
                     const std::array<Layer, 3>& layers,
                     const std::array<cluster_utils::LayerClusters, 3>& clusters) {
    std::cout << "Event " << evt << ":\n";
    std::cout << "  Hit channels (time): ";
    for (size_t i = 0; i < eventChannels.size(); ++i) {
        std::cout << eventChannels[i] << "(" << eventTimes[i] << ")";
        if (i + 1 != eventChannels.size()) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
    for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
        const auto& layerClusters = clusters[layerIdx];
        const auto& layerHits = layers[layerIdx];
        std::cout << "  Layer " << layerIdx << " hits: eta1=" << layerHits.eta1Hits.size()
                  << " eta2=" << layerHits.eta2Hits.size() << std::endl;
        size_t totalHitsLayer = layerHits.eta1Hits.size() + layerHits.eta2Hits.size();
        size_t totalClustersLayer = layerClusters.eta1.size() + layerClusters.eta2.size();
        std::cout << "    Total hits (eta1+eta2): " << totalHitsLayer << std::endl;
        std::cout << "    Total clusters (eta1+eta2): " << totalClustersLayer << std::endl;
        printClusters(layerIdx, layerClusters.eta1, "eta1");
        printClusters(layerIdx, layerClusters.eta2, "eta2");
    }
}
