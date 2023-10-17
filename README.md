# scan-to-map-scan-matching-tests
Trilogy of methods for correspondenceless matching of scans from 2D panoramic LIDAR sensors to map-scans

You may run x1, DBH, and FSM  for each dataset with

./s2msm_node 2 100 0 7373  0.05 0.035 0 0 0 0 360 360 SKG aces
./s2msm_node 1 100 0 4933  0.05 0.035 0 0 0 0 360 360 DBH fr079
./s2msm_node 1 100 0 13630 0.05 0.035 0 0 0 0 360 360 FMT intel
./s2msm_node 2 100 0 1987  0.05 0.035 0 0 0 0 360 360 SKG mit_csail
./s2msm_node 1 100 0 17479 0.05 0.035 0 0 0 0 360 360 FMT mit_killian
