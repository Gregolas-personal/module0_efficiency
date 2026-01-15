// Include header file with function declarations
#include "event_display.h"

// ROOT framework includes for GUI and plotting
#include <TApplication.h>   // ROOT application for GUI event loop
#include <TCanvas.h>         // Canvas for drawing plots
#include <TF1.h>             // 1D function for fitting
#include <TGraph.h>          // 2D graph for scatter plots
#include <TGraph2D.h>        // 3D graph for 3D scatter plots
#include <TH1F.h>            // 1D histogram
#include <TH2F.h>            // 2D histogram for frame/axes
#include <TH3F.h>            // 3D histogram for frame/axes
#include <TLegend.h>         // Legend for plots
#include <TPolyLine3D.h>     // 3D polyline for drawing tracks and planes
#include <TString.h>         // ROOT string utilities
#include <TSystem.h>         // System interface for event processing

// Standard C++ library includes
#include <algorithm>  // For std::min, std::max
#include <cmath>      // For mathematical functions
#include <iostream>   // For input/output
#include <memory>     // For smart pointers (unique_ptr)
#include <string>     // For std::string
#include <vector>     // For std::vector containers

// Custom header for cluster utilities
#include "cluster_utils.h"

// Forward declaration: prints prompt for user interaction between events
void printEventDisplayPrompt(size_t evt);
// Forward declaration: prints 3D line fit parameters (eta and time difference)
void print3DLineFit(size_t evt,
                    double etaIntercept, double etaSlope, double etaChi2,
                    double dtIntercept, double dtSlope, double dtChi2);

// Main function: runs event-by-event display and performs linear fits
// Parameters:
//   clustersPerEvent: vector of cluster data for each event (3 layers per event)
//   showDisplay: whether to show interactive GUI plots
//   appPtr: pointer to ROOT application for GUI
//   layerSpacingCm: spacing between detector layers in cm
void runPerEventDisplay(
    const std::vector<std::array<cluster_utils::LayerClusters, 3>>& clustersPerEvent,
    bool showDisplay,
    TApplication* appPtr,
    double layerSpacingCm) {
    // Create histogram to track chi-squared values from cluster center fits
    TH1F* hChi2CentersFit = new TH1F("hChi2CentersFit",
                                     "#chi^{2} of center linear fit;#chi^{2};Entries",
                                     100, 0.0, 20.0);  // 100 bins from 0 to 20
    // Create histogram to track chi-squared values from time difference fits
    TH1F* hChi2DtFit = new TH1F("hChi2DtFit",
                                "#chi^{2} of #Delta t linear fit;#chi^{2};Entries",
                                100, 0.0, 20.0);  // 100 bins from 0 to 20
    // Smart pointer to canvas (auto-deleted when out of scope)
    std::unique_ptr<TCanvas> cEvt;
    if (showDisplay) {  // Only create canvas if display is enabled
        // Create canvas with 3 pads: 1500px wide, 600px tall
        cEvt = std::make_unique<TCanvas>("cEvt", "Event display", 1500, 600);
        cEvt->Divide(3, 1);  // Divide into 3 columns, 1 row
    }
    // Loop over all events in the dataset
    for (size_t evt = 0; evt < clustersPerEvent.size(); ++evt) {
        // Initialize vectors for display (contains all clusters, no filtering)
        std::vector<double> zEta1Display;           // Z positions of eta1 clusters
        std::vector<double> centersEta1Display;     // Center channels of eta1 clusters
        std::vector<double> zEta2Display;           // Z positions of eta2 clusters
        std::vector<double> centersEta2Display;     // Center channels of eta2 clusters
        std::vector<double> zMatchedDisplay;        // Z positions of matched cluster pairs
        std::vector<double> timeDiffMatchedDisplay; // Time difference (eta1 - eta2) for matched pairs
        std::vector<double> etaMatchedDisplay;      // Average center channel for matched pairs
        
        // Initialize vectors for fitting (filtered by layer multiplicity)
        std::vector<double> zEta1Fit;           // Z positions of eta1 clusters for fitting
        std::vector<double> centersEta1Fit;     // Center channels of eta1 clusters for fitting
        std::vector<double> zEta2Fit;           // Z positions of eta2 clusters for fitting
        std::vector<double> centersEta2Fit;     // Center channels of eta2 clusters for fitting
        std::vector<double> zMatchedFit;        // Z positions of matched pairs for fitting
        std::vector<double> timeDiffMatchedFit; // Time differences for matched pairs (fitting)
        std::vector<double> etaMatchedFit;      // Average centers for matched pairs (fitting)
        
        size_t eligibleFitLayers = 0;  // Count layers suitable for fitting (1-3 total clusters)

        // Loop over the 3 detector layers
        for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
            // Calculate Z position of eta2 layer (layers are equally spaced)
            double zEta2Layer = layerIdx * layerSpacingCm;
            // Eta1 layer is 2 mm (0.2 cm) above eta2 layer
            double zEta1Layer = zEta2Layer + 0.2;
            // Get references to cluster collections for this layer
            const auto& eta1Clusters = clustersPerEvent[evt][layerIdx].eta1;
            const auto& eta2Clusters = clustersPerEvent[evt][layerIdx].eta2;
            // Calculate total number of clusters in this layer
            const size_t totalClusters = eta1Clusters.size() + eta2Clusters.size();
            // Only use layers with 1-3 total clusters for fitting (reduces noise)
            const bool includeLayerForFit = (totalClusters > 0 && totalClusters <= 3);
            if (includeLayerForFit) {
                ++eligibleFitLayers;  // Increment count of layers suitable for fitting
            }
            // Process all eta1 clusters in this layer
            for (const auto& c1 : eta1Clusters) {
                // Add Z position to display vector
                zEta1Display.push_back(zEta1Layer);
                // Add center channel position to display vector
                centersEta1Display.push_back(static_cast<double>(c1.centerChannel));
                // If this layer is suitable for fitting, also add to fit vectors
                if (includeLayerForFit) {
                    zEta1Fit.push_back(zEta1Layer);
                    centersEta1Fit.push_back(static_cast<double>(c1.centerChannel));
                }
            }
            // Process all eta2 clusters in this layer
            for (const auto& c2 : eta2Clusters) {
                // Add Z position to display vector
                zEta2Display.push_back(zEta2Layer);
                // Add center channel position to display vector
                centersEta2Display.push_back(static_cast<double>(c2.centerChannel));
                // If this layer is suitable for fitting, also add to fit vectors
                if (includeLayerForFit) {
                    zEta2Fit.push_back(zEta2Layer);
                    centersEta2Fit.push_back(static_cast<double>(c2.centerChannel));
                }
            }

            // Find matched cluster pairs between eta1 and eta2
            for (const auto& c1 : eta1Clusters) {
                for (const auto& c2 : eta2Clusters) {
                    // Calculate spatial separation between cluster centers
                    int delta = std::abs(c1.centerChannel - c2.centerChannel);
                    // Clusters are considered matched if within 2 channels of each other
                    if (delta <= 2) {
                        // Use eta1 Z position for matched points (includes vertical offset)
                        double zMatch = zEta2Layer + 0.2;
                        // Store Z position of match
                        zMatchedDisplay.push_back(zMatch);
                        // Calculate and store time difference (eta1 - eta2)
                        timeDiffMatchedDisplay.push_back(c1.time - c2.time);
                        // Store average center channel of the matched pair
                        etaMatchedDisplay.push_back(0.5 * (c1.centerChannel + c2.centerChannel));
                        // If this layer qualifies for fitting, add to fit vectors
                        if (includeLayerForFit) {
                            zMatchedFit.push_back(zMatch);
                            timeDiffMatchedFit.push_back(c1.time - c2.time);
                            etaMatchedFit.push_back(0.5 * (c1.centerChannel + c2.centerChannel));
                        }
                    }
                }
            }
        }  // End loop over layers

        // Prepare combined dataset for linear fit (both eta1 and eta2 points)
        std::vector<double> zAll;         // Combined Z positions
        std::vector<double> centersAll;   // Combined center channels
        // Pre-allocate memory for efficiency
        zAll.reserve(zEta1Fit.size() + zEta2Fit.size());
        centersAll.reserve(centersEta1Fit.size() + centersEta2Fit.size());
        // Insert all eta1 Z positions
        zAll.insert(zAll.end(), zEta1Fit.begin(), zEta1Fit.end());
        // Insert all eta2 Z positions
        zAll.insert(zAll.end(), zEta2Fit.begin(), zEta2Fit.end());
        // Insert all eta1 center channels
        centersAll.insert(centersAll.end(), centersEta1Fit.begin(), centersEta1Fit.end());
        // Insert all eta2 center channels
        centersAll.insert(centersAll.end(), centersEta2Fit.begin(), centersEta2Fit.end());

        // Perform linear fit on cluster centers vs Z position
        std::unique_ptr<TGraph> gCenters;  // Graph for center positions
        std::unique_ptr<TF1> fCenters;     // Fit function (linear)
        // Only fit if we have at least 3 eligible layers and 3+ data points
        if (eligibleFitLayers >= 3 && zAll.size() >= 3) {
            // Create graph with all combined points (X=Z, Y=center channel)
            gCenters = std::make_unique<TGraph>(static_cast<int>(zAll.size()), zAll.data(), centersAll.data());
            // Set unique name for this event's graph
            gCenters->SetName(Form("gAll_evt%zu", evt));
            // Set graph appearance (black line and markers)
            gCenters->SetLineColor(kBlack);
            gCenters->SetLineWidth(2);
            gCenters->SetMarkerColor(kBlack);
            // Create linear fit function (pol1 = a + b*x)
            fCenters = std::make_unique<TF1>(Form("f_eta_evt%zu", evt), "pol1",
                                             -1.0, layerSpacingCm * 3 + 1.0);
            // Perform fit (Q=quiet, 0=don't draw automatically)
            gCenters->Fit(fCenters.get(), "Q0");
            // Record chi-squared value in histogram
            hChi2CentersFit->Fill(fCenters->GetChisquare());
        }

        // Perform linear fit on time difference vs Z position for matched pairs
        std::unique_ptr<TGraph> gDt;  // Graph for time differences
        std::unique_ptr<TF1> fDt;     // Fit function (linear)
        // Only fit if we have at least 3 eligible layers and 3+ matched pairs
        if (eligibleFitLayers >= 3 && zMatchedFit.size() >= 3) {
            // Create graph with matched pairs (X=Z, Y=time difference)
            gDt = std::make_unique<TGraph>(static_cast<int>(zMatchedFit.size()),
                                           zMatchedFit.data(), timeDiffMatchedFit.data());
            // Set marker style (triangle)
            gDt->SetMarkerStyle(22);
            // Set marker color (green)
            gDt->SetMarkerColor(kGreen + 2);
            // Calculate data range for fit function limits
            double xmin, ymin, xmax, ymax;
            gDt->ComputeRange(xmin, ymin, xmax, ymax);
            // Create linear fit function over data range
            fDt = std::make_unique<TF1>(Form("f_dt_evt%zu", evt), "pol1", xmin, xmax);
            // Set fit line color (green)
            fDt->SetLineColor(kGreen + 2);
            // Perform fit (Q=quiet, 0=don't draw)
            gDt->Fit(fDt.get(), "Q0");
            // Record chi-squared value in histogram
            hChi2DtFit->Fill(fDt->GetChisquare());
        }

        // Perform linear fit on average eta position vs Z for matched pairs
        std::unique_ptr<TGraph> gEtaMatched;  // Graph for matched eta positions
        std::unique_ptr<TF1> fEtaMatched;     // Fit function (linear)
        // Only fit if we have at least 3 eligible layers and 2+ matched pairs
        if (eligibleFitLayers >= 3 && zMatchedFit.size() >= 2) {
            // Create graph with matched pairs (X=Z, Y=average center channel)
            gEtaMatched = std::make_unique<TGraph>(static_cast<int>(zMatchedFit.size()),
                                                   zMatchedFit.data(), etaMatchedFit.data());
            // Set marker style (open triangle)
            gEtaMatched->SetMarkerStyle(24);
            // Set marker color (magenta)
            gEtaMatched->SetMarkerColor(kMagenta + 2);
            // Calculate data range for fit function limits
            double xmin, ymin, xmax, ymax;
            gEtaMatched->ComputeRange(xmin, ymin, xmax, ymax);
            // Create linear fit function over data range
            fEtaMatched = std::make_unique<TF1>(Form("f_etaMatched_evt%zu", evt), "pol1", xmin, xmax);
            // Set fit line color (magenta)
            fEtaMatched->SetLineColor(kMagenta + 2);
            // Perform fit (Q=quiet, 0=don't draw)
            gEtaMatched->Fit(fEtaMatched.get(), "Q0");
        }

        // If both fits exist, print the 3D line fit parameters
        if (fEtaMatched && fDt) {
            print3DLineFit(evt,
                           fEtaMatched->GetParameter(0),  // Eta intercept
                           fEtaMatched->GetParameter(1),  // Eta slope
                           fEtaMatched->GetChisquare(),   // Eta chi-squared
                           fDt->GetParameter(0),          // Time diff intercept
                           fDt->GetParameter(1),          // Time diff slope
                           fDt->GetChisquare());          // Time diff chi-squared
        }

        // If display is disabled, skip visualization and continue to next event
        if (!showDisplay) {
            continue;  // Fits and chi2 histograms are still filled
        }

        // === Pad 1: Cluster centers vs layer Z position ===
        cEvt->cd(1);  // Select first pad in canvas
        gPad->Clear();  // Clear any previous content
        // Create 2D frame histogram for axes (X=Z position, Y=center channel)
        TH2F* frameCenters = new TH2F(Form("frameCenters_%zu", evt),
                                      "Cluster centers;Layer z (cm);Center channel",
                                      10, -1.0, layerSpacingCm * 3 + 1.0,  // X axis: -1 to ~15 cm
                                      41, -0.5, 40.5);  // Y axis: -0.5 to 40.5 channels
        frameCenters->SetStats(false);  // Don't show statistics box
        frameCenters->Draw();  // Draw frame to establish axes
        // Initialize graph pointers for eta1 and eta2 clusters
        TGraph* gEta1 = nullptr;
        TGraph* gEta2 = nullptr;
        // Draw eta1 clusters if any exist
        if (!zEta1Display.empty()) {
            // Create graph from eta1 data (X=Z, Y=center channel)
            gEta1 = new TGraph(static_cast<int>(zEta1Display.size()), zEta1Display.data(), centersEta1Display.data());
            gEta1->SetName(Form("gEta1_evt%zu", evt));  // Unique name
            gEta1->SetMarkerStyle(20);  // Filled circle
            gEta1->SetMarkerColor(kRed);  // Red color
            gEta1->Draw("P SAME");  // Draw points on same pad
        }
        // Draw eta2 clusters if any exist
        if (!zEta2Display.empty()) {
            // Create graph from eta2 data (X=Z, Y=center channel)
            gEta2 = new TGraph(static_cast<int>(zEta2Display.size()), zEta2Display.data(), centersEta2Display.data());
            gEta2->SetName(Form("gEta2_evt%zu", evt));  // Unique name
            gEta2->SetMarkerStyle(21);  // Filled square
            gEta2->SetMarkerColor(kBlue);  // Blue color
            gEta2->Draw("P SAME");  // Draw points on same pad
        }
        // Create legend in top-right corner (NDC coordinates: 0.60-0.88 in X, 0.65-0.88 in Y)
        TLegend* legCenters = new TLegend(0.60, 0.65, 0.88, 0.88);
        // Add event number as legend header (no symbol)
        legCenters->AddEntry((TObject*)nullptr, Form("Event %zu", evt), "");
        // Add eta1 entry to legend if graph exists
        if (gEta1) {
            legCenters->AddEntry(gEta1, "Eta1 center", "p");  // "p" = point marker
        }
        // Add eta2 entry to legend if graph exists
        if (gEta2) {
            legCenters->AddEntry(gEta2, "Eta2 center", "p");  // "p" = point marker
        }
        // Draw combined fit if it exists
        if (gCenters && fCenters) {
            gCenters->Draw("P SAME");  // Draw fit graph points
            fCenters->SetLineColor(kBlack);  // Black fit line
            fCenters->SetLineWidth(2);  // Thick line
            fCenters->Draw("SAME");  // Draw fit line
            // Add fit parameters to legend
            legCenters->AddEntry(fCenters.get(),
                                 Form("Combined fit: y=%.2f+%.2fx, #chi^{2}=%.2f",
                                      fCenters->GetParameter(0),    // Intercept
                                      fCenters->GetParameter(1),    // Slope
                                      fCenters->GetChisquare()),    // Chi-squared
                                 "l");  // "l" = line
        }
        legCenters->Draw();  // Draw the legend

        // === Pad 2: Time difference of matched clusters vs layer Z ===
        cEvt->cd(2);  // Select second pad in canvas
        gPad->Clear();  // Clear any previous content
        // Create 2D frame histogram for axes (X=Z position, Y=time difference)
        TH2F* frameDt = new TH2F(Form("frameDt_%zu", evt),
                                 "Matched cluster #Delta t (eta1 - eta2);Layer z (cm);#Delta t (ns)",
                                 10, -1.0, layerSpacingCm * 3 + 1.0,  // X axis: -1 to ~15 cm
                                 200, -100.0, 100.0);  // Y axis: -100 to 100 ns
        frameDt->SetStats(false);  // Don't show statistics box
        frameDt->Draw();  // Draw frame to establish axes
        // Draw matched cluster time differences if any exist
        if (!zMatchedDisplay.empty()) {
            // Create graph from matched pairs (X=Z, Y=time difference)
            TGraph* gDtDisplay = new TGraph(static_cast<int>(zMatchedDisplay.size()),
                                            zMatchedDisplay.data(), timeDiffMatchedDisplay.data());
            gDtDisplay->SetMarkerStyle(22);  // Triangle marker
            gDtDisplay->SetMarkerColor(kGreen + 2);  // Green color
            gDtDisplay->Draw("P SAME");  // Draw points on same pad
        }
        // Draw fit and legend if fit exists
        if (gDt) {
            // Create legend in top-right corner
            TLegend* legDt = new TLegend(0.60, 0.75, 0.88, 0.88);
            // Add event number as header
            legDt->AddEntry((TObject*)nullptr, Form("Event %zu", evt), "");
            // Draw fit function and add to legend if it exists
            if (fDt) {
                fDt->Draw("SAME");  // Draw fit line
                // Add fit parameters to legend
                legDt->AddEntry(fDt.get(),
                                Form("dt fit: y=%.2f+%.2fx, #chi^{2}=%.2f",
                                     fDt->GetParameter(0),    // Intercept
                                     fDt->GetParameter(1),    // Slope
                                     fDt->GetChisquare()),    // Chi-squared
                                "l");  // "l" = line
            }
            legDt->Draw();  // Draw the legend
        }

        // === Pad 3: 3D view (eta vs time difference vs Z) ===
        cEvt->cd(3);  // Select third pad in canvas
        gPad->Clear();  // Clear any previous content
        // Define 3D plot axis limits
        const double etaMin = -0.5;    // Min center channel
        const double etaMax = 40.5;    // Max center channel
        const double dtMin = -100.0;   // Min time difference (ns)
        const double dtMax = 100.0;    // Max time difference (ns)
        const double zMin = -0.5;      // Min Z position (cm)
        const double zMax = layerSpacingCm * 3 + 0.5;  // Max Z position (cm)
        // Create 3D frame histogram for axes
        TH3F* frame3D = new TH3F(Form("frame3D_%zu", evt),
                                 "Matched clusters 3D;#eta (center channel);#Delta t (ns);Layer z (cm)",
                                 10, etaMin, etaMax,  // X axis (eta)
                                 10, dtMin, dtMax,    // Y axis (time diff)
                                 10, zMin, zMax);     // Z axis (layer position)
        frame3D->SetStats(false);  // Don't show statistics box
        frame3D->Draw("BOX");  // Draw 3D box frame
        gPad->SetTheta(25);  // Set viewing angle theta (vertical)
        gPad->SetPhi(35);    // Set viewing angle phi (horizontal)

        // Draw matched cluster points in 3D if any exist
        if (!zMatchedDisplay.empty()) {
            // Create 3D graph with space for all matched points
            TGraph2D* g3d = new TGraph2D(static_cast<int>(zMatchedDisplay.size()));
            // Fill graph with all matched cluster data
            for (size_t i = 0; i < zMatchedDisplay.size(); ++i) {
                // Set point coordinates (X=eta, Y=time diff, Z=layer position)
                g3d->SetPoint(static_cast<int>(i),
                              etaMatchedDisplay[i],        // X: center channel
                              timeDiffMatchedDisplay[i],   // Y: time difference
                              zMatchedDisplay[i]);         // Z: layer position
            }
            g3d->SetMarkerStyle(20);  // Filled circle
            g3d->SetMarkerColor(kGreen + 2);  // Green color
            g3d->Draw("P SAME");  // Draw points in 3D
        }

        // Draw layer planes as rectangles in 3D space
        for (size_t layerIdx = 0; layerIdx < 3; ++layerIdx) {
            // Calculate Z position of this layer (at eta1 position, +0.2 cm offset)
            double zPlane = layerIdx * layerSpacingCm + 0.2;
            // Create closed rectangle with 5 points (last = first to close)
            TPolyLine3D* layerPlane = new TPolyLine3D(5);
            // Bottom-left corner
            layerPlane->SetPoint(0, etaMin, dtMin, zPlane);
            // Bottom-right corner
            layerPlane->SetPoint(1, etaMax, dtMin, zPlane);
            // Top-right corner
            layerPlane->SetPoint(2, etaMax, dtMax, zPlane);
            // Top-left corner
            layerPlane->SetPoint(3, etaMin, dtMax, zPlane);
            // Close the rectangle (back to bottom-left)
            layerPlane->SetPoint(4, etaMin, dtMin, zPlane);
            layerPlane->SetLineColor(kGray + 1);  // Gray color
            layerPlane->SetLineStyle(2);  // Dashed line
            layerPlane->Draw("SAME");  // Draw on same 3D plot
        }

        // Draw fitted 3D track line if both fits exist
        if (fEtaMatched && fDt) {
            // Find min and max Z values from matched fit data
            double zStart = zMatchedFit.front();  // Initialize with first point
            double zEnd = zMatchedFit.front();
            for (double z : zMatchedFit) {
                zStart = std::min(zStart, z);  // Find minimum Z
                zEnd = std::max(zEnd, z);      // Find maximum Z
            }
            // Create 3D line with 2 endpoints
            TPolyLine3D* trackLine = new TPolyLine3D(2);
            // Start point: evaluate both fit functions at zStart
            trackLine->SetPoint(0, fEtaMatched->Eval(zStart),  // X from eta fit
                                   fDt->Eval(zStart),          // Y from time diff fit
                                   zStart);                    // Z coordinate
            // End point: evaluate both fit functions at zEnd
            trackLine->SetPoint(1, fEtaMatched->Eval(zEnd),    // X from eta fit
                                   fDt->Eval(zEnd),            // Y from time diff fit
                                   zEnd);                      // Z coordinate
            trackLine->SetLineColor(kBlack);  // Black color
            trackLine->SetLineWidth(2);       // Thick line
            trackLine->Draw("SAME");  // Draw on same 3D plot
        }
        cEvt->Update();  // Update canvas to show all drawings
        gSystem->ProcessEvents();  // Process GUI events (make window responsive)
        printEventDisplayPrompt(evt);  // Print prompt asking user to continue
        std::string tmp;  // Temporary string for user input
        std::getline(std::cin, tmp);  // Wait for user to press Enter
    }  // End event loop
    
    // Keep the GUI alive after all events if display was shown
    if (showDisplay && appPtr) {
        appPtr->Run(kTRUE);  // Run ROOT event loop (kTRUE = return on termination signal)
    }

    // Save chi-squared distributions to PNG files (filled regardless of display setting)
    // Create canvas for center fit chi-squared histogram
    TCanvas* cChi2Centers = new TCanvas("cChi2Centers", "Chi2 of center fits", 800, 600);
    hChi2CentersFit->Draw();  // Draw the histogram
    cChi2Centers->SaveAs("chi2_centers_fit.png");  // Save to file
    // Create canvas for time difference fit chi-squared histogram
    TCanvas* cChi2Dt = new TCanvas("cChi2Dt", "Chi2 of dt fits", 800, 600);
    hChi2DtFit->Draw();  // Draw the histogram
    cChi2Dt->SaveAs("chi2_dt_fit.png");  // Save to file
}  // End of runPerEventDisplay function
