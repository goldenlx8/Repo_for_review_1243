## whole genome
cnvkit.py heatmap -d -o outputs/CNVKit/cnv_heatmap_call_cns.pdf outputs/CNVKit/*.call.cns -d
cnvkit.py heatmap -o outputs/CNVKit/cnv_heatmap_spontanous_tumor_call_cns.pdf \
outputs/CNVKit/N1307.call.cns \
outputs/CNVKit/N1637.call.cns \
outputs/CNVKit/N1690.call.cns \
outputs/CNVKit/N1255.call.cns \
outputs/CNVKit/N1369.call.cns \
outputs/CNVKit/N1314.call.cns \
outputs/CNVKit/N1345.call.cns \
outputs/CNVKit/N1010.call.cns \
outputs/CNVKit/N1019.call.cns \
outputs/CNVKit/N1062.call.cns \
outputs/CNVKit/N1108.call.cns \
outputs/CNVKit/N1285.call.cns \
outputs/CNVKit/N1319.call.cns \
outputs/CNVKit/N1692.call.cns

cnvkit.py heatmap -o outputs/CNVKit/cnv_heatmap_celllines.pdf \
outputs/CNVKit/N1011.call.cns \
outputs/CNVKit/N1076.call.cns \
outputs/CNVKit/N1018.call.cns \
outputs/CNVKit/N1128.call.cns


## Pten
cnvkit.py heatmap -d -o outputs/CNVKit/pten_heatmap.png outputs/CNVKit/N[0-9][0-9][0-9][0-9].cnr -c chr19:32734897-32803560

## Trp53
cnvkit.py heatmap -d -o outputs/CNVKit/Trp53_heatmap.png outputs/CNVKit/N[0-9][0-9][0-9][0-9].cnr -c  chr11:69471185-69482698
