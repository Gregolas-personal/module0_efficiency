g++ -o plot_eff eff_plot_from_hit.cxx cluster_utils.cxx `root-config --cflags --libs`
./plot_eff *.root
