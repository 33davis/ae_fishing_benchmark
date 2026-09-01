This folder contains all the synthetic generate files used in my analysis. 

All files are split into subfolders based on their fragment count values **_n_**, e.g. n10, n25, n50...etc
  -Within each of these subfolders are the the three types .fa files:
    1. Iteratively fragmented (length determined by derived parameters) file with **_n_** fragments per species
    2. A "selected_$number_**_n_**" file, where the original file were subsampled without replacement for deamination damage (3,7,18) with gargammel program
    3. A "non_selected_$number_**_n_**" file that is the counterpart non-subsampled file to be merged back with the deaminated treatment fragments (3,7,18) post-gargammel
