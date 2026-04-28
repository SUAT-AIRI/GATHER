"""
Answer prompt templates for GATHER cell-type annotation.

Usage:
    from answer_prompts import KG_ANSWER_PROMPT_TEMPLATE

    prompt = KG_ANSWER_PROMPT_TEMPLATE.format(
        cell_sentence="CD3D CD3E CD4 ...",
        kg_evidence="<retrieved knowledge>"
    )
"""

# ============================================================
# 34 Candidate Cell Types
# ============================================================

CANDIDATE_CELL_TYPES_TEXT = """1. CD16-negative, CD56-bright natural killer cell, human
2. CD16-positive, CD56-dim natural killer cell, human
3. CD4-positive helper T cell
4. CD8-positive, alpha-beta memory T cell
5. CD8-positive, alpha-beta memory T cell, CD45RO-positive
6. T follicular helper cell
7. alpha-beta T cell
8. alveolar macrophage
9. classical monocyte
10. conventional dendritic cell
11. dendritic cell, human
12. effector memory CD4-positive, alpha-beta T cell
13. effector memory CD8-positive, alpha-beta T cell, terminally differentiated
14. erythroid lineage cell
15. gamma-delta T cell
16. germinal center B cell
17. group 3 innate lymphoid cell
18. lymphocyte
19. macrophage
20. mast cell
21. megakaryocyte
22. memory B cell
23. mucosal invariant T cell
24. naive B cell
25. naive thymus-derived CD4-positive, alpha-beta T cell
26. naive thymus-derived CD8-positive, alpha-beta T cell
27. non-classical monocyte
28. plasma cell
29. plasmablast
30. plasmacytoid dendritic cell
31. precursor B cell
32. pro-B cell
33. progenitor cell
34. regulatory T cell"""


# ============================================================
# GATHER Final Answer Prompt Template
# ============================================================
# Placeholders:
#   {cell_sentence} - The marker genes expression profile
#   {kg_evidence}   - Knowledge retrieved by GATHER

KG_ANSWER_PROMPT_TEMPLATE = """You are an expert computational biologist specializing in single-cell transcriptomics and cell type annotation.

## Task
Identify the most specific cell type for a human immune cell based on its gene expression profile and the retrieved knowledge evidence.

## Input: Cell Sentence (ranked by expression specificity)
{cell_sentence}

## Candidate Cell Types (34 immune cell types)
""" + CANDIDATE_CELL_TYPES_TEXT + """

## Reference Knowledge
{kg_evidence}

## RULE
Your answer MUST be EXACTLY one of the 34 cell types listed above. Any other answer is INVALID.
If the knowledge mentions a cell type NOT in the list, you MUST find the closest match from the 34 candidates.

## Output Format
[Predicted_Cell_Type: <your answer>]
"""


# ============================================================
# Lung Dataset: 30 Candidate Cell Types
# ============================================================

LUNG_CANDIDATE_CELL_TYPES_TEXT = """1. macrophage
2. basal cell
3. classical monocyte
4. club cell
5. non-classical monocyte
6. capillary endothelial cell
7. basophil
8. CD8-positive, alpha-beta T cell
9. CD4-positive, alpha-beta T cell
10. vein endothelial cell
11. lung microvascular endothelial cell
12. adventitial cell
13. dendritic cell
14. intermediate monocyte
15. pericyte
16. endothelial cell of artery
17. neutrophil
18. plasma cell
19. natural killer cell
20. B cell
21. bronchial smooth muscle cell
22. vascular associated smooth muscle cell
23. endothelial cell of lymphatic vessel
24. smooth muscle cell
25. pulmonary ionocyte
26. plasmacytoid dendritic cell
27. mesothelial cell
28. serous cell of epithelium of bronchus
29. fibroblast
30. myofibroblast cell"""


# ============================================================
# Lung: GATHER Prompt Template
# ============================================================

KG_ANSWER_PROMPT_TEMPLATE_LUNG = """You are an expert computational biologist specializing in single-cell transcriptomics and cell type annotation.

## Task
Identify the most specific cell type for a human lung cell based on its gene expression profile and the retrieved knowledge evidence.

## Input: Cell Sentence (ranked by expression specificity)
{cell_sentence}

## Candidate Cell Types (30 lung cell types)
""" + LUNG_CANDIDATE_CELL_TYPES_TEXT + """

## Reference Knowledge
{kg_evidence}

## RULE
Your answer MUST be EXACTLY one of the 30 cell types listed above. Any other answer is INVALID.
If the knowledge mentions a cell type NOT in the list, you MUST find the closest match from the 30 candidates.

## Output Format
[Predicted_Cell_Type: <your answer>]
"""

