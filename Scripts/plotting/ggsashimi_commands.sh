# These are the ggsashimi commands used to plot the splicing outliers in Figure 4
# For each standard and depleted reference, there are 3 representative samples plotted
# https://github.com/guigolab/ggsashimi


#this is for sample with unknown causal mutation
#chr10:102890678:102900593
#AS3MT
./ggsashimi.py \
  -b bam16_samples.tsv \
  -g gencode.v48.chr_patch_hapl_scaff.annotation.gtf \
  -c chr10:102890078-102900993 \
  -j s16_splice.bed \
  -s NONE \
  -M 5 \
  -O 3 \
  -C 3 \
  -L 3 \
  --height 2 \
  --ann-height 2.5 \
  --shrink \
  -o sashimi_chr10_102890678_102900593_multi \
  -P palette.txt \
  --alpha 0.9 \
  --width=8 \
  -F png

#these are for sample 5-4
#chr3:174440532-174551069
#NAALADL2
./ggsashimi.py \
  -b bam_samples.tsv \
  -g gencode.v48.chr_patch_hapl_scaff.annotation.gtf \
  -c chr3:174440532-174551069 \
  -s NONE \
  -M 10 \
  -O 3 \
  -C 3 \
  -L 3 \
  --height 2 \
  --ann-height 2.5 \
  --shrink \
  --fix-y-scale \
  -o sashimi_chr3_174441032_174550569_multi \
  -P palette.txt \
  --alpha 0.9 \
  --width=8 \
  -F png

#these are for sample 3-2 RNUATAC 
#chr21:37816250:37840428
#KCNJ6
./ggsashimi.py \
  -b bam32_samples.tsv \
  -g gencode.v48.chr_patch_hapl_scaff.annotation.gtf \
  -c chr21:37815250-37841428 \
  -s NONE \
  -M 10 \
  -O 3 \
  -C 3 \
  -L 3 \
  --height 2 \
  --ann-height 3 \
  --shrink \
  --fix-y-scale \
  -o sashimi_chr21_37816250_37840428_multi \
  -P palette.txt \
  --alpha 0.9 \
  --width=8 \
  -F png


