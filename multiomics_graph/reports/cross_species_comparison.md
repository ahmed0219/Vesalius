# Cross-Species & Cross-Strain Comparison

> Generated from saved pipeline outputs. Layers: G=genome, T=transcriptome, P=proteome, M=metabolome.
> All 16 strains have real eggNOG-mapper COG annotations and (except E. coli MS_14385) a transcriptome layer.

## 1. Per-strain overview

| Species | Strain | ST | Layers | CDS | RNA genes | Proteins | Aligned | Gene nodes | Protein nodes | Metab. nodes | Annot. nodes | Edges | Hyperedges | AMR (loci) | COG cats |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Escherichia coli | B36 | ST131 | G+T+P+M | 5130 | 5186 | 2335 | 2283 | 5118 | 2283 | 219 | 2833 | 77321 | 3699 | 8 (8) | 21 |
| Escherichia coli | MS_14384 | ST963 | G+T+P+M | 4904 | 4946 | 2396 | 2348 | 4891 | 2348 | 219 | 2833 | 83248 | 3767 | 1 (1) | 20 |
| Escherichia coli | MS_14385 | - | G+P+M | 5322 | 0 | 2393 | 2348 | 5296 | 2348 | 219 | 2828 | 47148 | 3762 | 0 (0) | 21 |
| Escherichia coli | MS_14386 | ST224 | G+T+P+M | 4821 | 4844 | 2356 | 2310 | 4808 | 2310 | 219 | 2843 | 78064 | 3736 | 13 (14) | 20 |
| Escherichia coli | MS_14387 | ST69 | G+T+P+M | 5016 | 5041 | 2283 | 2246 | 5000 | 2246 | 219 | 2831 | 79161 | 3661 | 3 (4) | 20 |
| Klebsiella pneumoniae | AJ218 | ST2121 | G+T+P | 5433 | 5529 | 2378 | 2284 | 5420 | 2284 | 0 | 2158 | 77047 | 3302 | 4 (4) | 20 |
| Klebsiella pneumoniae | KPC2 | ST258 | G+T+P | 5523 | 5622 | 2001 | 1921 | 5509 | 1921 | 0 | 2157 | 76107 | 2934 | 9 (9) | 20 |
| Klebsiella variicola | 0331100710 | - | G+T+P | 5091 | 5175 | 1913 | 1837 | 5089 | 1837 | 0 | 1049 | 63548 | 2383 | 0 (0) | 20 |
| Klebsiella variicola | 04153260899 | - | G+T+P | 5082 | 5164 | 1624 | 1575 | 5081 | 1575 | 0 | 1035 | 59789 | 2117 | 0 (0) | 20 |
| Klebsiella variicola | AJ055 | - | G+T+P | 5393 | 5483 | 2251 | 2165 | 5390 | 2165 | 0 | 1057 | 69244 | 2727 | 0 (0) | 20 |
| Klebsiella variicola | AJ_292 | - | G+T+P | 5112 | 5212 | 2187 | 2104 | 5109 | 2104 | 0 | 1058 | 66230 | 2665 | 0 (0) | 20 |
| Staphylococcus aureus | BPH2760 | ST1 | G+T+P | 2613 | 2735 | 1321 | 1295 | 2612 | 1295 | 0 | 1500 | 39108 | 1961 | 1 (1) | 21 |
| Staphylococcus aureus | BPH2819 | ST5 | G+T+P | 2657 | 2771 | 1306 | 1280 | 2656 | 1280 | 0 | 1509 | 36819 | 1957 | 2 (2) | 21 |
| Staphylococcus aureus | BPH2900 | ST22 | G+T+P | 2732 | 2822 | 1265 | 1238 | 2728 | 1238 | 0 | 1493 | 36886 | 1905 | 3 (3) | 21 |
| Staphylococcus aureus | BPH2947 | ST239 | G+T+P | 3042 | 3180 | 1251 | 1233 | 3040 | 1233 | 0 | 1510 | 39734 | 1911 | 12 (14) | 21 |
| Staphylococcus aureus | BPH2986 | ST8 (USA300) | G+T+P | 2888 | 3019 | 1255 | 1230 | 2887 | 1230 | 0 | 1502 | 39270 | 1904 | 8 (8) | 21 |

## 2. Edge counts by relation type

| Relation | B36 | MS_14384 | MS_14385 | MS_14386 | MS_14387 | AJ218 | KPC2 | 0331100710 | 04153260899 | AJ055 | AJ_292 | BPH2760 | BPH2819 | BPH2900 | BPH2947 | BPH2986 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| gene--genomic_proximity--gene | 4159 | 4017 | 4285 | 3951 | 4082 | 4428 | 4441 | 4152 | 4137 | 4406 | 4178 | 2248 | 2275 | 2324 | 2587 | 2447 |
| gene--transcriptional_correlation--gene | 33718 | 32374 | 0 | 31227 | 32635 | 36364 | 36950 | 33923 | 33768 | 36091 | 34155 | 17224 | 16855 | 17721 | 19831 | 18906 |
| gene--encodes--protein | 2283 | 2348 | 2348 | 2310 | 2246 | 2284 | 1921 | 1837 | 1575 | 2165 | 2104 | 1295 | 1280 | 1238 | 1233 | 1230 |
| protein--abundance_correlation--protein | 15916 | 15734 | 16631 | 16166 | 15547 | 16277 | 13342 | 12989 | 11004 | 14870 | 14248 | 8692 | 8292 | 8200 | 8167 | 8372 |
| protein--ppi--protein | 3711 | 6540 | 4978 | 5237 | 5225 | 2843 | 4015 | 92 | 64 | 81 | 75 | 630 | 265 | 208 | 310 | 422 |
| gene--in_pathway--annotation | 7978 | 8480 | 8212 | 8156 | 8116 | 8654 | 7550 | 7056 | 6322 | 7956 | 7856 | 4412 | 4624 | 4406 | 4378 | 4448 |
| gene--cog_category--annotation | 2370 | 2449 | 2422 | 2402 | 2345 | 2416 | 2016 | 1943 | 1668 | 2276 | 2233 | 1345 | 1331 | 1282 | 1276 | 1271 |
| protein--annotated_by--annotation | 4126 | 7541 | 5070 | 5296 | 5618 | 3208 | 5138 | 1444 | 1164 | 1294 | 1285 | 3007 | 1751 | 1381 | 1762 | 1973 |
| gene--associated_with--annotation | 8 | 1 | 0 | 14 | 4 | 4 | 9 | 0 | 0 | 0 | 0 | 1 | 2 | 3 | 14 | 8 |
| protein--member_of_cluster--annotation | 451 | 710 | 525 | 579 | 564 | 569 | 725 | 112 | 87 | 105 | 96 | 254 | 144 | 123 | 176 | 193 |
| metabolite--in_pathway--annotation | 1673 | 1673 | 1673 | 1673 | 1673 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| metabolite--metabolised_by--protein | 345 | 798 | 421 | 470 | 523 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| metabolite--abundance_correlation--metabolite | 583 | 583 | 583 | 583 | 583 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## 3. Hyperedge counts by family

| Family | B36 | MS_14384 | MS_14385 | MS_14386 | MS_14387 | AJ218 | KPC2 | 0331100710 | 04153260899 | AJ055 | AJ_292 | BPH2760 | BPH2819 | BPH2900 | BPH2947 | BPH2986 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| multi_omics_triple | 2283 | 2348 | 2348 | 2310 | 2246 | 2284 | 1921 | 1837 | 1575 | 2165 | 2104 | 1295 | 1280 | 1238 | 1233 | 1230 |
| kegg_pathway | 372 | 380 | 376 | 380 | 374 | 386 | 380 | 336 | 332 | 352 | 350 | 254 | 264 | 254 | 260 | 258 |
| cog_category | 20 | 20 | 20 | 20 | 20 | 20 | 19 | 19 | 19 | 19 | 20 | 20 | 20 | 20 | 20 | 20 |
| go_term | 779 | 779 | 779 | 779 | 779 | 595 | 595 | 178 | 178 | 178 | 178 | 378 | 378 | 378 | 378 | 378 |
| ppi_cluster | 10 | 10 | 10 | 10 | 10 | 13 | 13 | 13 | 13 | 13 | 13 | 13 | 13 | 13 | 13 | 13 |
| amr_mechanism | 7 | 1 | 0 | 7 | 3 | 4 | 6 | 0 | 0 | 0 | 0 | 1 | 2 | 2 | 7 | 5 |
| metabolic_pathway | 228 | 229 | 229 | 230 | 229 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## 4. AMR determinants per strain

- **B36** (Escherichia coli): tet(A) (tetracycline), dfrA17 (trimethoprim), aadA5 (aminoglycoside), sul1 (sulfonamide), mph(A) (macrolide), blaCTX-M-15 (beta_lactam), blaOXA-1 (beta_lactam), aac(6')-Ib-cr5 (fluoroquinolone)
- **MS_14384** (Escherichia coli): blaCMY-2 (beta_lactam)
- **MS_14385** (Escherichia coli): none
- **MS_14386** (Escherichia coli): blaTEM-1b (beta_lactam), blaTEM-215 (beta_lactam), aac(3)-IId (aminoglycoside), aadA1 (aminoglycoside), aph(3')-Ia (aminoglycoside), strAB (aminoglycoside), sul2 (sulfonamide), sul3 (sulfonamide), dfrA12 (trimethoprim), dfrA14 (trimethoprim), floR (phenicol), fosA4 (fosfomycin), mph(A) (macrolide)
- **MS_14387** (Escherichia coli): blaTEM-1b (beta_lactam), strAB (aminoglycoside), sul2 (sulfonamide)
- **AJ218** (Klebsiella pneumoniae): aadA1 (aminoglycoside), blaSHV-44 (beta_lactam), tet(D) (tetracycline), sul1 (sulfonamide)
- **KPC2** (Klebsiella pneumoniae): aac(6')-Ib (aminoglycoside), aadA2 (aminoglycoside), aph(3')-Ia (aminoglycoside), blaKPC-2 (beta_lactam), blaSHV-12 (beta_lactam), catA1 (phenicol), dfrA12 (trimethoprim), mph(A) (macrolide), sul1 (sulfonamide)
- **0331100710** (Klebsiella variicola): none
- **04153260899** (Klebsiella variicola): none
- **AJ055** (Klebsiella variicola): none
- **AJ_292** (Klebsiella variicola): none
- **BPH2760** (Staphylococcus aureus): blaZ (beta_lactam)
- **BPH2819** (Staphylococcus aureus): blaZ (beta_lactam), fosB (fosfomycin)
- **BPH2900** (Staphylococcus aureus): blaZ (beta_lactam), erm(C) (macrolide), mecA (beta_lactam)
- **BPH2947** (Staphylococcus aureus): aac(6')-Ie-aph(2'')-Ia (aminoglycoside), ant(6)-Ia (aminoglycoside), ant(9)-Ia (aminoglycoside), aph(3')-IIIa (aminoglycoside), blaZ (beta_lactam), dfrG (trimethoprim), erm(A) (macrolide), fosB (fosfomycin), mecA (beta_lactam), sat (streptothricin), tet(K) (tetracycline), tet(M) (tetracycline)
- **BPH2986** (Staphylococcus aureus): ant(6)-Ia (aminoglycoside), aph(3')-IIIa (aminoglycoside), blaZ (beta_lactam), fosB (fosfomycin), mecA (beta_lactam), mph(C) (macrolide), msr(A) (macrolide), sat (streptothricin)

## 5. AMR class coverage by species

| Species | strains | total markers | total loci | classes |
|---|---|---|---|---|
| Escherichia coli | 5 | 25 | 27 | aminoglycoside (6), beta_lactam (6), fluoroquinolone (1), fosfomycin (1), macrolide (2), phenicol (1), sulfonamide (4), tetracycline (1), trimethoprim (3) |
| Klebsiella pneumoniae | 2 | 13 | 13 | aminoglycoside (4), beta_lactam (3), macrolide (1), phenicol (1), sulfonamide (2), tetracycline (1), trimethoprim (1) |
| Klebsiella variicola | 4 | 0 | 0 |  |
| Staphylococcus aureus | 5 | 26 | 28 | aminoglycoside (6), beta_lactam (8), fosfomycin (3), macrolide (4), streptothricin (2), tetracycline (2), trimethoprim (1) |
