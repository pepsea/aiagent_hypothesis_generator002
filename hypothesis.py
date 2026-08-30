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
If direct evidence is limited, also consider the "Related Literature" sections for pathway-connected genes and describe how those genes link to {gene} in {disease} context.
Additionally, use the "Disease Pathway Analysis" section: if {gene} is a direct member of enriched disease pathways, name those pathways and the co-pathway disease genes (e.g., "APP co-occurs with PSEN1, APOE in the synaptic signaling pathway enriched in {disease}"). If {gene} is not a direct member, describe how it may connect to {disease} biology indirectly through shared pathway partners.
Confidence: Low / Moderate / High / Very High — one sentence reason.

### 1c. Clinical Relevance
List existing drugs targeting {gene} with phase and mechanism. State whether any overlap with {disease}.

### 1d. Expression & Network
State the top expressing tissues from GTEx (TPM values) and HPA. Note subcellular location and protein class.
Describe key PPI partners and pathway memberships relevant to {disease}.
Also reference the "Disease Pathway Analysis" section: state the pathway_overlap_score, name the top disease pathways, and whether {gene} is a direct member.

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

## 6. Existing Drug Landscape & Competitive Analysis
- **Known drugs/candidates:** summarise drugs/candidates for {gene} with phases. Repositioning or combination opportunities.
- **Competitive landscape:** using the "Competitive Landscape" evidence (ClinicalTrials.gov trial count and new-entrant risk rating), state how many trials target this gene/disease pair, list the key sponsors/companies, and quote the new-entrant risk level (HIGH/MODERATE/LOW) given in the evidence.
- **Differentiation strategy:** given that risk level, what would a new entrant need to do differently (e.g. better selectivity, different modality, combination approach, unaddressed patient subgroup) to succeed despite the competitive pressure? If risk is LOW, state why this may indicate limited validation rather than a clear opportunity.

---

## 7. Safety Assessment
- On-target risks from gene function
- Safety-relevant tissue expression (heart/liver/kidney/CNS) with TPM or HPA level
- gnomAD pLI/LOEUF interpretation (high = greater on-target risk)
- Off-target / mechanism-based toxicity signals

---

## 8. Recommended Experiments

First, identify the 2–3 most critical gaps that currently block clinical advancement of {gene} as a target in {disease}. For each gap, state what is unknown and why it is a blocker (one sentence each).

Then propose 3–5 experiments in priority order to address those gaps. At least 2 experiments must directly address the clinical advancement blockers identified above. For each experiment state:
1. **Gap addressed** — which clinical blocker this experiment resolves (one sentence).
2. **Purpose** — the specific scientific question being answered (one sentence).
3. **Experimental system** — the most appropriate model or sample type (e.g., disease-relevant cell line, patient-derived primary cells, ex vivo tissue, organoid). Prefer in vitro and ex vivo over animal models. Name a specific system if mentioned in the evidence.
4. **Key readout** — the main measurable endpoint and how a positive result would de-risk the clinical hypothesis.

Keep each experiment to 4–5 sentences. Do not write full protocols.

Prioritise: (1) target/variant validation in disease-relevant tissue, (2) pharmacological proof-of-concept with a tool compound, (3) mechanism/biomarker confirmation, (4) translational evidence from patient samples.

---

## 9. Key Uncertainties
- Evidence gaps
- Alternative hypotheses
- Major risks

Do NOT write a "## References" section — it is generated automatically
from the evidence context after your response, with working links.
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
直接エビデンスが少ない場合は「Related Literature」セクション（パスウェイ隣接遺伝子の論文）も参照し、その遺伝子が{gene}を通じて{disease}にどう関連するかを記述。
さらに「Disease Pathway Analysis」セクションを活用する：{gene}が疾患エンリッチメントパスウェイの直接メンバーであれば、そのパスウェイ名と共存する疾患関連遺伝子（Co-pathway disease genes）を挙げ、{disease}との機能的関連を論じる。直接メンバーでない場合は、共有パスウェイパートナーを通じた間接的な関与を説明する。
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

## 6. 既存薬景観と競合分析
- **既存薬・候補:** {gene}標的薬・候補をフェーズ付きで整理。リポジショニング・併用の可能性。
- **競合状況:** エビデンス中の「Competitive Landscape」（ClinicalTrials.govの試験数・新規参入リスク評価）を用い、この遺伝子×疾患ペアを狙う試験が何件あるか、主要なスポンサー（企業名）、新規参入リスク（HIGH/MODERATE/LOW）を具体的に引用すること。
- **差別化戦略:** そのリスクレベルを踏まえ、新規参入者が成功するために何で差別化すべきか（選択性向上・異なるモダリティ・併用療法・未対応の患者サブグループ等）を述べる。リスクがLOWの場合は、それが明確な機会を意味するのか、単に検証が不十分なだけかにも触れること。

---

## 7. 安全性・毒性リスク評価
- 遺伝子機能に基づくオンターゲットリスク
- 安全性関連組織における発現（心臓/肝臓/腎臓/CNS）とTPM値またはHPAレベル
- gnomAD pLI/LOEUF（高値＝オンターゲットリスク大）の解釈
- オフターゲット・機序由来毒性シグナル

---

## 8. 推奨次期実験

まず、{gene}を{disease}の治療標的として臨床開発を進める上で現在不足している情報のうち、最も重大なブロッカーを2〜3項目特定する。各ブロッカーについて「何が不明で、なぜ臨床化の障壁になるか」を1文で明記すること。

次に、それらのブロッカーを解消するための実験を優先順位の高い順に3〜5件提案する。うち最低2件は上記の臨床化ブロッカーに直接対応するものとする。各実験について以下を記述する：
1. **解消するギャップ** — この実験がどのブロッカーを解決するか（1文）。
2. **目的** — 答えるべき具体的な科学的問い（1文）。
3. **試験系** — 最適なモデルまたはサンプル種別（例：疾患関連細胞株、患者由来初代細胞、ex vivo組織、オルガノイド）。動物モデルよりin vitro・ex vivoを優先する。エビデンスに具体的な系が記載されていればそれを使う。
4. **主要エンドポイントと判定基準** — 主な測定項目と、陽性結果が臨床仮説をどのようにde-riskするか。

各実験は4〜5文程度。詳細なプロトコルは不要。

優先順位：(1) 疾患関連組織での標的/変異の検証、(2) ツール化合物を用いた薬理学的概念実証、(3) メカニズム/バイオマーカーの確認、(4) 患者サンプルを用いたトランスレーショナルエビデンス。

---

## 9. 主要な不確実性・限界
- エビデンスのギャップ
- 代替仮説
- 仮説に対する主要リスク

「## 参考文献」セクションは書かないこと — この後、エビデンスコンテキストから
リンク付きで自動生成される。
"""

def generate_hypothesis(
    gene: str,
    disease: str,
    context: str,
    llm_client,
    temperature: float = 0.3,
    lang: str = "en",
    stream_callback=None,
    num_ctx: int = 24576,
) -> str:
    template = HYPOTHESIS_PROMPT_JA if lang == "ja" else HYPOTHESIS_PROMPT_TEMPLATE
    prompt = template.format(gene=gene, disease=disease, context=context)
    kwargs = dict(temperature=temperature, max_tokens=4000, num_ctx=num_ctx)
    if stream_callback is not None:
        kwargs["stream_callback"] = stream_callback
    return llm_client.generate(prompt, **kwargs)
