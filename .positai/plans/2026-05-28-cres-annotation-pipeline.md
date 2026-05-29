# Plano: CREs Annotation Pipeline (Humano + Camundongo)

## Objetivo
Criar um pipeline completo que:
1. Anota cCREs humanos e murinos com genes associados
2. Identifica pares homólogos entre espécies
3. Extrai e alinha sequências de cCREs e regiões promotoras correspondentes
4. Adiciona anotações evolutivas (PhastCons, RepeatMasker)
5. Calcula identidade de sequência entre pares homólogos

---

## Escopo do Trabalho

### Fase 1: Dados Base e Estrutura (JÁ FEITO)
- [x] Carregar `homologous_cres` (hg38-mm10 mapeamento)
- [x] Carregar `cres_human` (GRCh38-cCREs.bed)
- [x] Carregar `cres_mouse` (mm10-cCREs.bed)
- [x] Join das três tabelas → `cres_annotated`
- [x] Extrair sequências de cCREs humanos (BSgenome)

### Fase 2: Anotação Gênica (PARCIAL - Precisa Organizar)
**Humano:**
- [ ] Converter `cres_human` → GRanges com seqinfo correto
- [ ] Usar ChIPseeker + annotatePeak para ligar a genes
- [ ] Adicionar colunas: gene_name, gene_id, annotation, distance_to_tss
- [ ] Extrair sequências promotoras customizadas (upstream 2500bp, downstream 300bp)

**Camundongo:**
- [ ] Repetir para `cres_mouse` com TxDb.Mmusculus.UCSC.mm10.knownGene

### Fase 3: Alinhamento de Sequências
**Para cada par homólogo (human_enhancer ↔ mouse_enhancer):**
- [ ] Extrair sequência do cCRE humano
- [ ] Extrair sequência do cCRE murino
- [ ] Alinhar localmente (pairwise alignment com Biostrings::pairwiseAlignment)
- [ ] Calcular % identidade (sequence_identity)
- [ ] Calcular similarity score

**Para regiões promotoras:**
- [ ] Extrair sequência promotora humana (gene associado)
- [ ] Extrair sequência promotora de camundongo (homólogo de gene)
- [ ] Alinhar e calcular identidade de forma similar

### Fase 4: Anotações Evolutivas e Estruturais
**Conservação (PhastCons):**
- [ ] Consultar AnnotationHub para scores de conservação (vertebrados)
- [ ] Extrair score médio de PhastCons para cada cCRE
- [ ] Adicionar colunas: phastCons_score_human, phastCons_score_mouse

**Elementos Transponíveis (RepeatMasker):**
- [ ] Consultar AnnotationHub para RepeatMasker (hg38 e mm10)
- [ ] Encontrar overlaps entre cCREs e TEs
- [ ] Adicionar colunas: has_te_human (bool), te_classes_human (comma-separated)
- [ ] Repetir para mouse

### Fase 5: Tabela Final
**Estrutura de `cres_final`:**
```
- Identificadores: human_enhancer, mouse_enhancer, homologous (0/1)
- Localização humana: chrom_human, start_human, end_human, type_human
- Localização mouse: chrom_mouse, start_mouse, end_mouse, type_mouse
- Genes humanos: gene_name_human, gene_id_human, distance_to_tss_human, annotation_human
- Genes mouse: gene_name_mouse, gene_id_mouse, distance_to_tss_mouse, annotation_mouse
- Sequências: sequence_ccre_human, sequence_ccre_mouse
- Sequências promotoras: sequence_promoter_human, sequence_promoter_mouse
- Alinhamentos: ccre_alignment_identity, promoter_alignment_identity
- Conservação: phastCons_score_human, phastCons_score_mouse
- TEs: has_te_human, has_te_mouse, te_classes_human, te_classes_mouse
```

---

## Problemas no Código Atual

1. **IDs do AnnotationHub são placeholders** — `AH99341`, `AH75265` são fictícios; precisam ser descobertos via `query()`
2. **Alinhamento não implementado** — Código não faz alinhamento de sequências, apenas as extrai
3. **Falta de tratamento de NAs** — Join `full_join` pode gerar muitos NAs (cCREs únicos em uma espécie)
4. **Promotores vs CREs confundidos** — O código extrai promotores de genes, mas depois faz overlap com cCREs; clareza sobre qual range usar
5. **Falta de Bioconductor::biomaRt** — Para mapear IDs entre espécies (humano ↔ mouse gene ortólogos)

---

## Ordem de Implementação Recomendada

1. **Setup & Discovery** — Instalar pacotes, descobrir IDs corretos no AnnotationHub
2. **Data Preparation** — Validar cres_annotated, tratar NAs, criar GRanges robustos
3. **Annotation (Human)** — ChIPseeker, sequências, PhastCons, TEs
4. **Annotation (Mouse)** — Repetir com genomas/DBs murinos
5. **Sequence Alignment** — Biostrings para pairwise alignment
6. **Final Assembly** — Juntar tudo em tabela final limpa
7. **QC & Export** — Validação, salvar como .tsv/.parquet

---

## Dependências Bioconductor

```r
BiocManager::install(c(
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Mmusculus.UCSC.mm10",
  "AnnotationHub",
  "rtracklayer",
  "Biostrings"
))
```

---

## Próximos Passos

✅ **Aprovado?** Iniciar Fase 1 (reorganizar código base) → Fase 2 (anotação gênica) → continuar sequencialmente.

❓ **Dúvidas antes de começar:**
1. Para os **pares homólogos sem cCRE em uma espécie** → Filtrar fora ou manter como NA?
2. Para **alinhamento de sequências** → Usar gap penalty padrão ou customizado?
3. Para **ranges promotores** → Usar sempre 2500 upstream / 300 downstream, ou variar?
4. Qual nível de **PhastCons** preferir? (Ex: 100 vertebrados, mamíferos, primatas?)
