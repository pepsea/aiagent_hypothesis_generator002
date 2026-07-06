"""Hypothesis generation — builds prompt and calls LLM."""

# ──────────────────────────────────────────────────────────────────────────────
# English prompt
# ──────────────────────────────────────────────────────────────────────────────
HYPOTHESIS_PROMPT_TEMPLATE = """You are a drug discovery scientist. Write a hypothesis report for {gene} in {disease} using ONLY the evidence below.

STRICT RULES:
1. Only use reference tags that ACTUALLY APPEAR in the context (e.g. [Paper 1], [UniProt 1], [GWAS 1]). Never invent a tag.
2. If data is absent for a point, write "No data available" — never fabricate.
3. Use specific values from the context: gene names, drug names, p-values, TPM values, pLI scores, etc.
4. Cite every factual claim with a tag from the context.

--- EVIDENCE ---
{context}
--- END EVIDENCE ---

Write the report below. Use actual data from the evidence. Do not copy instructions.

## 1. Target Validity

### 1a. Genetic Evidence
State the number of GWAS hits and ClinVar variants found. Quote specific traits, p-values, and variant names with citations.
Confidence: Low / Moderate / High / Very High — one sentence reason.

### 1b. Functional Evidence
Summarize what the literature (abstract summaries) says about {gene} in {disease}. Quote specific findings.
Confidence: Low / Moderate / High / Very High — one sentence reason.

### 1c. Clinical Relevance
List existing drugs targeting {gene} with phase and mechanism. State whether any overlap with {disease}.

### 1d. Expression & Network
State the top expressing tissues from GTEx (TPM values) and HPA. Note subcellular location and protein class.
Describe key PPI partners and pathway memberships relevant to {disease}.

### 1e. Overall Validity Score
**Overall: Low / Moderate / High / Very High**
Two-sentence synthesis.

---

## 2. Molecular Mechanism
Describe the molecular mechanism linking {gene} to {disease}:
- **Protein function:** what does {gene} normally do at the molecular level?
- **Dysregulation in disease:** how is {gene} expression/activity/structure altered in {disease}? (overexpressed / loss-of-function / mutation / mislocalised)
- **Direct molecular effects:** what immediate molecular events result (e.g. kinase activity change, protein–protein interaction disruption, transcription factor binding)?
- **Key interactors:** which PPI partners or pathway nodes are most affected? Cite network data.
- **Downstream signalling cascade:** step-by-step from {gene} perturbation → effector proteins → cellular phenotype → tissue pathology → {disease} manifestation.
- **Supporting evidence:** cite specific findings from literature, GWAS, ClinVar, and pathway enrichment.

---

## 3. Disease Mechanism
Step-by-step narrative: how does {gene} dysfunction lead to the tissue-level and clinical pathology of {disease}? Connect the molecular events above to organ function and patient symptoms. Use specific pathway names, effectors, and cite literature.

---

## 4. Therapeutic Hypothesis

### 4a. Treatment Hypothesis (bullet points)
- **Intervention:** proposed modulation of {gene} (inhibit / activate / degrade / replace)
- **Effect on target:** how would this change {gene} activity?
- **Downstream effect:** how does that correct the molecular mechanism above?
- **Clinical outcome:** expected improvement in {disease}
- **Key evidence:** cite the strongest supporting data points
- **Falsifiable prediction:** one testable statement that could confirm or refute this

### 4b. One-sentence hypothesis
"If [intervention on {gene}] then [outcome in {disease} patients] because [mechanism]."

---

## 5. Modality
Best modality (small molecule / antibody / PROTAC / ASO / gene therapy) and why:
- Rationale based on subcellular location and druggability
- Tissue specificity strategy (selectivity vs. delivery advantage)
- AlphaFold pLDDT and structural confidence
- Key advantage over alternatives

---

## 6. Existing Drug Landscape
Summarise known drugs/candidates for {gene} with phases. Repositioning or combination opportunities.

---

## 7. Safety Assessment
- On-target risks from gene function
- Safety-relevant tissue expression (heart/liver/kidney/CNS) with TPM or HPA level
- gnomAD pLI/LOEUF interpretation (high = greater on-target risk)
- Off-target / mechanism-based toxicity signals

---

## 8. Recommended Experiments
| Experiment | Endpoint | Expected Result |
|---|---|---|
| (3–5 specific experiments based on the evidence above) | | |

---

## 9. Key Uncertainties
- Evidence gaps
- Alternative hypotheses
- Major risks

---

## References
List only tags that appear in the evidence context and were cited above:

### Papers
### Disease Databases
### Gene/Protein Databases
### Drug/Safety Databases
"""

# ──────────────────────────────────────────────────────────────────────────────
# Japanese prompt
# ──────────────────────────────────────────────────────────────────────────────
HYPOTHESIS_PROMPT_JA = """あなたは創薬の専門家です。以下のエビデンスのみを使って、{gene}を標的とした{disease}の創薬仮説レポートを日本語で作成してください。

【絶対ルール】
1. コンテキスト内に実際に登場する引用タグのみ使うこと（例: [Paper 1]、[UniProt 1]、[GWAS 1]）。存在しないタグを作らない。
2. データがない項目は「データなし」と書く。情報を推測・捏造しない。
3. コンテキストの具体的な数値・名称・スコアを使う（遺伝子名、薬剤名、p値、TPM値、pLIスコアなど）。
4. 事実的な主張には必ずコンテキスト内のタグを引用する。

--- エビデンス ---
{context}
--- エビデンスここまで ---

以下のレポートを書いてください。エビデンスの実際のデータを使用すること。指示文をそのまま書き写さないこと。

## 1. ターゲット妥当性評価

### 1a. 遺伝的エビデンス
GWASヒット数・ClinVar変異数を明記。具体的な形質名・p値・変異名を引用付きで記述。
エビデンス信頼度：低 / 中 / 高 / 非常に高 — 理由を一文で。

### 1b. 機能的エビデンス
文献アブストラクトが{gene}と{disease}についてどう述べているかを要約。具体的知見を引用。
エビデンス信頼度：低 / 中 / 高 / 非常に高 — 理由を一文で。

### 1c. 臨床的関連性
{gene}を標的とする既存薬を列挙（薬剤名・フェーズ・作用機序）。{disease}との治療領域の重複を述べる。

### 1d. 発現プロファイル・ネットワーク
GTEx（TPM値）・HPAから上位発現組織を具体的に記述。細胞内局在・タンパク質クラスを明記。
主要PPI相互作用パートナーと{disease}関連パスウェイを記述。

### 1e. ターゲット妥当性総合スコア
**総合：低 / 中 / 高 / 非常に高**
2文以内で根拠をまとめる。

---

## 2. 分子メカニズム
{gene}と{disease}をつなぐ分子メカニズムを詳述する：
- **タンパク質の正常機能：** {gene}は分子レベルで何をしているか？
- **疾患における機能異常：** {disease}において{gene}の発現・活性・構造はどう変化しているか？（過剰発現 / 機能喪失 / 変異 / 局在異常）
- **直接的な分子イベント：** それによりどのような即時の分子イベントが起こるか？（キナーゼ活性変化、タンパク質間相互作用の破綻、転写因子結合など）
- **主要な相互作用パートナー：** どのPPIパートナーやパスウェイノードが最も影響を受けるか？ネットワークデータを引用。
- **下流シグナル伝達カスケード：** {gene}の変化 → エフェクタータンパク質 → 細胞表現型 → 組織病理 → {disease}発症という流れをステップごとに記述。
- **支持エビデンス：** 文献・GWAS・ClinVar・パスウェイエンリッチメントの具体的な知見を引用。

---

## 3. 疾患メカニズムの考察
上記の分子イベントが組織・臓器レベルの病態と患者症状にどうつながるかを記述。{gene}機能異常から{disease}の臨床症状に至るステップをパスウェイ名・エフェクター・文献引用とともに説明。

---

## 4. 治療仮説

### 4a. 疾患治療仮説（箇条書き）
- **介入方法：** {gene}への具体的介入（阻害/活性化/分解/補充など）
- **ターゲットへの効果：** 介入による{gene}活性・発現の変化
- **下流効果：** 上記の分子メカニズム異常をどう是正するか
- **期待される臨床効果：** {disease}の症状・進行への改善
- **支持エビデンス：** 最も強いエビデンスを引用タグで
- **検証可能な予測：** 反証可能な1つの予測

### 4b. 仮説一文
「{gene}に対して〔介入〕を行うと、{disease}患者において〔臨床効果〕が得られる。これは〔分子メカニズム〕による。」

---

## 5. モダリティ提案
最適モダリティ（低分子/抗体/PROTAC/ASO/遺伝子療法）と理由：
- 細胞内局在・ドラッガビリティに基づく根拠
- 組織特異性（選択性戦略 or デリバリー優位性）
- AlphaFold pLDDT・構造信頼度
- 他モダリティに対する優位性

---

## 6. 既存薬景観とリポジショニング
{gene}標的薬・候補をフェーズ付きで整理。リポジショニング・併用の可能性。

---

## 7. 安全性・毒性リスク評価
- 遺伝子機能に基づくオンターゲットリスク
- 安全性関連組織における発現（心臓/肝臓/腎臓/CNS）とTPM値またはHPAレベル
- gnomAD pLI/LOEUF（高値＝オンターゲットリスク大）の解釈
- オフターゲット・機序由来毒性シグナル

---

## 8. 推奨次期実験
| 実験種別 | エンドポイント | 期待される結果 |
|---|---|---|
| （上記エビデンスに基づく具体的な実験を3〜5件） | | |

---

## 9. 主要な不確実性・限界
- エビデンスのギャップ
- 代替仮説
- 仮説に対する主要リスク

---

## 参考文献
コンテキスト内に登場し、上記で引用したタグのみ列挙：

### 論文（PubMed）
### 疾患データベース
### 遺伝子・タンパク質データベース
### 薬剤・安全性データベース
"""

def generate_hypothesis(
    gene: str,
    disease: str,
    context: str,
    llm_client,
    temperature: float = 0.3,
    lang: str = "en",
    stream_callback=None,
) -> str:
    template = HYPOTHESIS_PROMPT_JA if lang == "ja" else HYPOTHESIS_PROMPT_TEMPLATE
    prompt = template.format(gene=gene, disease=disease, context=context)
    kwargs = dict(temperature=temperature, max_tokens=4000)
    if stream_callback is not None:
        kwargs["stream_callback"] = stream_callback
    return llm_client.generate(prompt, **kwargs)
