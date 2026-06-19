- ! where it should be 
- "freq" or "band" or "spct"
- smwh_modl or smwh_ft
- nvlp, sidepeak etc.. naming
- 162 164 use symbols
- functional, map, functors
- function placement
- current centers and rois
  in unit of pixels
  smwh_roi = (40, 80)
  smwh_essn = (30, 60)
  smwh_core = (20, 40)
  - smwh_roi: determined by the entire dataset (per run-info), roi centered by the entire center
  - smwh_essn: currently unused
  - calculating essenstials: 
    - smwh = smwh_roi, center = smwh_roi .+ 1: choose the same center 
      - used for modulation profile
      - used for dens2d = dens_roi, which is used for profile fitting
    - smwh_strip = smwh = smwh_roi (default), used for prfl_strip
      - which is unused for now
    - smwh_core, cent_core=xy: the center now determined per IB per rep
      - for desn2d_core, used in PCA
