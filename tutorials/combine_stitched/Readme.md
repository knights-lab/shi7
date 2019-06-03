```
#!/usr/bin/env bash

# Run the paired end
!shi7 -i ../tests/testfq/GA01NEX/ -o ./data/pe --adaptor Nextera

# Run the single end
!shi7 -i ../tests/testfq/GA01NEX/ -o ./data/se --adaptor Nextera --flash False

# Combine them together
python ./combine_stitched_se.py --input_se ./data/se/combined_seqs.fna --input_pe data/pe/combined_seqs.fna --output ./data/test.fna
```