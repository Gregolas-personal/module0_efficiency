g++ -o plot_eff eff_plot_from_hit.cxx event_display.cxx printer.cxx cluster_utils.cxx `root-config --cflags --libs`
./plot_eff ./ly0/2025-12-11_23-36_triplet520_ly0_HV6000_vth19_Vamp15.root 
