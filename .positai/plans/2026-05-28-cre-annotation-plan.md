# Plano: Análise Comparativa de CREs Humano-Camundongo

## Objetivo
Anotar, comparar e alinhar CREs (candidate cis-Regulatory Elements) entre humano (hg38) e camundongo (mm10), calculando identidade de sequência, conservação filogenética (PhastCons) e presença de elementos transponíveis. Manter registros separados para CREs homólogos e não-homólogos para estatística descritiva.

---

## Escopo e Estrutura de Dados

### Arquivos de entrada
- `Genomic/hg38-mm10-Homologous.tsv` → 537,822 pares de CREs homólogos
- `Genomic/GRCh38-cCREs.bed` → ~1M CREs humanos
- `Genomic/mm10-cCREs.bed` → CREs de camundongo

### Saídas esperadas
1. **cres_annotated_homologous** → CREs com correspondências 1:1 ou N:M anotadas
2. **cres_annotated_human_only** → CREs humanos sem homólogo em camundongo
3. **cres_annotated_mouse_only** → CREs de camundongo sem homólogo humano
4. Tabelas de estatística descritiva (quantidades, tipos de CRE, etc.)

---

## Fases da Análise

### **Fase 1: Importação e Estruturação Básica**
- [ ] Importar 3 arquivos (homologous_cres, cres_human, cres_mouse)
- [ ] Adicionar sufixos às colunas (_human, _mouse)
- [ ] Full join para capturar homólogos + não-homólogos
- [ ] Criar flag `is_homologous` e `match_type` (1:1, 1:N, N:1, orphan)
- [ ] **Output:** `cres_annotated` (tabela unificada com ~2M linhas estimadas)

### **Fase 2: Anotação de Genes Associados**
**Ferramentas:** ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, TxDb.Mmusculus.UCSC.mm10.knownGene

#### Para humanos:
- [ ] Converter para GRanges com seqinfo correto
- [ ] `annotatePeak()` com tssRegion = c(-3000, 3000)
- [ ] Extrair: gene_name, gene_id, annotation (promoter/intronic/distal), distance_to_tss
- [ ] Registrar genes para cada CRE (pode ser N:1 se overlap múltiplos genes)

#### Para camundongos:
- [ ] Mesmo pipeline com mm10 TxDb
- [ ] Alinhar nomenclatura com homólogos humanos (quando possível via biomaRt)

### **Fase 3: Extração de Sequências**
**Ferramentas:** BSgenome.Hsapiens.UCSC.hg38, BSgenome.Mmusculus.UCSC.mm10

#### Para CREs:
- [ ] Extrair sequência de cada CRE (humano e camundongo)
- [ ] Armazenar em coluna `sequence_human`, `sequence_mouse`

#### Para regiões gênicas (exons/promotores):
- [ ] Para cada gene associado, extrair range customizado: **2500bp upstream + 300bp downstream**
- [ ] Armazenar como `sequence_gene_region`

### **Fase 4: Alinhamento e Cálculo de Identidade**
**Ferramentas:** Biostrings (pairwiseAlignment com gap penalty: -2/-1, local alignment)

#### Para CREs homólogos:
- [ ] Alinhar sequência human_CRE vs mouse_CRE (Smith-Waterman)
- [ ] Extrair: identity_percent, alignment_length, num_mismatches, num_gaps
- [ ] Registrar local alignment score

#### Para regiões gênicas:
- [ ] Alinhar sequence_human_gene vs sequence_mouse_gene
- [ ] Extrair métricas equivalentes → `gene_identity_percent`

### **Fase 5: Anotação de Conservação (PhastCons)**
**Ferramentas:** AnnotationHub (phastCons BigWig)

**Dataset:** **7 vertebrados** (mais restrito, alto sinal)
- [ ] Buscar em AnnotationHub: phastCons hg38 7-way vertebrates
- [ ] Importar arquivo BigWig para hg38
- [ ] Sobrepor CREs humanos → extrair mean/median phastcons score
- [ ] Para camundongo: usar score de conservação correspondente (mm10 7-way)
- [ ] Armazenar: `phastcons_score_human`, `phastcons_score_mouse`

### **Fase 6: Anotação de Elementos Transponíveis (RepeatMasker)**
**Ferramentas:** AnnotationHub (RepeatMasker GRanges)

#### Para humanos:
- [ ] Buscar RepeatMasker hg38 em AnnotationHub
- [ ] `findOverlaps()` entre CREs e TE regions
- [ ] Extrair para cada CRE:
  - `te_overlaps` (TRUE/FALSE)
  - `te_class` (SINE, LINE, DNA, LTR, etc.)
  - `te_family` (ex: Alu, L1, etc.)
  - `te_pct_coverage` (% do CRE coberto)

#### Para camundongos:
- [ ] Mesmo pipeline com RepeatMasker mm10

### **Fase 7: Organização Final e Estatística Descritiva**
- [ ] Separar em 3 tabelas finais:
  1. **cres_homologous_annotated** → 537K linhas (ou menos após limpeza)
  2. **cres_human_only_annotated** → human-específicos
  3. **cres_mouse_only_annotated** → mouse-específicos
  
- [ ] Gerar sumários:
  - Total CREs por espécie
  - % conservados (com homólogo)
  - % com homólogo e identidade >80%
  - Distribuição de tipos de CRE (PLS/pELS/dELS/CTCF-only)
  - Média PhastCons por grupo
  - % com elementos transponíveis

- [ ] Exportar tabelas como RDS (rápido) e CSV (compartilhável)

---

## Detalhes Técnicos

### Gap Penalty para Alinhamento
```r
pairwiseAlignment(
  seq1, seq2,
  type = "local",  # Smith-Waterman
  gapOpening = -2,
  gapExtension = -1
)
```

### PhastCons Dataset
```r
# No AnnotationHub:
ah <- AnnotationHub()
query(ah, c("phastCons", "Homo sapiens", "hg38", "7way"))
# Buscar por vertebrates 7-way
```

### Estrutura de Coluna (cres_annotated final)
```
human_enhancer, mouse_enhancer, is_homologous, match_type,
chrom_human, start_human, end_human, id_dhs_human, id_ccre_human, type_human,
chrom_mouse, start_mouse, end_mouse, id_dhs_mouse, id_ccre_mouse, type_mouse,
sequence_human, sequence_mouse,
gene_name_human, gene_id_human, annotation_human, distance_to_tss_human,
gene_name_mouse, gene_id_mouse, annotation_mouse, distance_to_tss_mouse,
sequence_gene_human, sequence_gene_mouse,
cre_identity_percent, cre_alignment_length, cre_gaps,
gene_identity_percent,
phastcons_score_human, phastcons_score_mouse,
te_overlap_human, te_class_human, te_family_human, te_coverage_human,
te_overlap_mouse, te_class_mouse, te_family_mouse, te_coverage_mouse
```

---

## Estimativas

| Fase | Complexidade | Tempo Estimado |
|------|---|---|
| 1. Estruturação | Baixa | 5–10 min |
| 2. Gene annotation | Média | 20–40 min (ChIPseeker + TxDb) |
| 3. Sequências | Média | 30–60 min (1M+ regiões) |
| 4. Alinhamento | Alta | 60–180 min (Biostrings, ~1.5M pares) |
| 5. PhastCons | Média | 20–40 min (BigWig overlap) |
| 6. RepeatMasker | Baixa | 10–20 min (overlap rápido) |
| 7. Resumo | Baixa | 10 min |
| **Total** | — | **3–5 horas** |

---

## Decisões Confirmadas ✓

1. **PhastCons:** **7 vertebrados** (mais restrito, alto sinal)
2. **Escopo de genes:** **2500bp upstream + 300bp downstream** (customizado)
3. **Alinhamento:** **Local (Smith-Waterman)** – flexível para regiões curtas
4. **Tratamento de N:M matches:** Manter todas as combinações, marcar `match_type`

---

## Próximas Ações
1. Implementar Fase 1 (estruturação)
2. Testar em subset (ex: chr1) antes de escala total
3. Iterativamente executar Fases 2–7
4. Gerar relatório final de conservação
