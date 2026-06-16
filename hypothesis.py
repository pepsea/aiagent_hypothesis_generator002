"""Hypothesis generation — builds prompt and calls LLM."""

# ──────────────────────────────────────────────────────────────────────────────
# English prompt
# ──────────────────────────────────────────────────────────────────────────────
HYPOTHESIS_PROMPT_TEMPLATE = """You are an expert medicinal chemist and computational biologist specializing in drug target identification and hypothesis generation.

Based on the structured evidence below, generate a detailed drug discovery hypothesis report for targeting {gene} in {disease}.

{context}

---

Produce a structured report using EXACTLY the sections below.
For every factual claim, cite the source using the reference tag provided in the context (e.g. [Paper 1], [ClinVar 2], [GWAS 1], [OpenTargets], [UniProt], etc.).
Use abstract summaries from the literature section to support mechanistic arguments.
Do NOT fabricate data or reference numbers not present in the evidence.

---

## 1. Target Validity Assessment

### 1a. Genetic Association Evidence
Evaluate genetic evidence linking {gene} to {disease}:
- GWAS findings: significant loci, effect sizes, lead SNPs [cite GWAS refs]
- ClinVar pathogenic/likely-pathogenic variants [cite ClinVar refs]
- Genetic evidence confidence: **Low / Moderate / High / Very High** — justify in one sentence

### 1b. Functional Association Evidence
- How does {gene} function relate to the biology of {disease}?
- Key mechanisms from the literature: cite and summarize relevant abstract findings [cite Paper refs]
- Disease mechanism of action (MoA): known or proposed molecular role of {gene} in disease pathophysiology
- Functional evidence confidence: **Low / Moderate / High / Very High**

### 1c. Clinical Relevance
- Existing drugs/clinical candidates targeting {gene} [cite ChEMBL refs]
- Stage of clinical development and therapeutic area overlap
- Patient population relevance and unmet medical need

### 1d. Gene Function, Expression & Network Context
- Core molecular function of {gene} (kinase/receptor/transcription factor etc.) [cite UniProt]
- Tissue expression profile: which tissues/organs express {gene} most highly, and relevance to {disease} tissue [cite GTEx, HPA]
- Protein subcellular localisation and protein class [cite HPA, UniProt]
- Key protein–protein interactions and their disease relevance [cite IntAct/SIGNOR refs if available]
- Pathway/GO enrichment context from network analysis [cite enrichment results if available]
- Overall network-level role in {disease}-relevant biology

### 1e. Overall Target Validity Score
**Overall: Low / Moderate / High / Very High**
One-paragraph synthesis integrating all four sub-assessments above.

---

## 2. Proposed Disease Mechanism
- Step-by-step mechanistic model: how does {gene} dysregulation → {disease} phenotype?
- Supporting evidence from literature abstracts [cite Paper refs]
- Upstream activators and downstream effectors
- Relevant biological pathways

---

## 3. Therapeutic Hypothesis

### 3a. Molecular Mechanism Hypothesis
Describe the proposed molecular mechanism in bullet points:
- **Dysregulation**: How is {gene} dysregulated or dysfunctional in {disease}? [cite evidence]
- **Key pathway**: Which downstream pathway/effector is most critically affected? [cite evidence]
- **Cellular consequence**: What cell-level phenotype does this produce? [cite evidence]
- **Tissue/organ consequence**: How does that cellular phenotype lead to {disease} pathology? [cite evidence]

### 3b. Disease Treatment Hypothesis
State the core therapeutic hypothesis as a structured argument:
- **Intervention**: What specific intervention on {gene} is proposed? (inhibition / activation / degradation / replacement / etc.)
- **Expected effect on target**: How would this intervention alter {gene} activity or expression?
- **Expected downstream effect**: How would that target-level change correct the molecular mechanism above?
- **Expected clinical outcome**: What improvement in {disease} symptoms or progression would result?
- **Supporting rationale**: Key evidence supporting this logic [cite Paper, GWAS, ClinVar, OpenTargets refs]
- **Testable prediction**: One falsifiable prediction that would confirm or refute this hypothesis

### 3c. Hypothesis Statement (one sentence)
Synthesise 3a–3b into a single testable sentence:
"If [specific intervention on {gene}] then [expected therapeutic outcome in {disease} patients] because [mechanistic rationale]."

---

## 4. Modality Recommendation
Recommend the best-fit drug modality (small molecule inhibitor/activator, antibody, PROTAC, ASO, gene therapy, etc.):
- Rationale aligned to target biology and subcellular localisation [cite UniProt, HPA]
- Tissue specificity considerations: if {gene} is broadly expressed, discuss selectivity strategy; if tissue-restricted, note delivery advantage [cite GTEx, HPA]
- AlphaFold structural confidence and druggability implications [cite AlphaFold if available]
- Key technical considerations
- Advantage over alternative modalities

---

## 5. Existing Drug Landscape & Repositioning
- Summary of drugs/candidates targeting {gene} [cite ChEMBL refs]
- Repositioning opportunities or combination strategies
- Competitive considerations

---

## 6. Safety & Toxicity Risk Assessment
- On-target safety risks from gene function [cite UniProt, PubChem]
- Tissue expression safety: list safety-relevant tissues where {gene} is expressed (heart, liver, kidney, CNS, reproductive organs) and discuss associated risks [cite GTEx, HPA]
- Population constraint: interpret pLI/LOEUF scores — high constraint implies greater on-target risk [cite gnomAD if available]
- Off-target / mechanism-based toxicity signals [cite adverse event data]
- Patient population risk factors

---

## 7. Recommended Next Experiments
List 3–5 prioritised experiments:
| Experiment | Endpoint | Expected Result |
|---|---|---|
| ... | ... | ... |

---

## 8. Key Uncertainties & Limitations
- Missing evidence gaps
- Alternative hypotheses
- Major risks to this hypothesis

---

## References

### Papers (PubMed)
(List every [Paper N] cited above, copied verbatim from the context References section)

### Disease & Genetic Databases
(List every [ClinVar N], [GWAS N], [OpenTargets] cited above)

### Gene & Protein Databases
(List every [UniProt], [IntAct], [SIGNOR] cited above)

### Drug & Safety Databases
(List every [ChEMBL N], [PubChem] cited above)
"""

# ──────────────────────────────────────────────────────────────────────────────
# Japanese prompt
# ──────────────────────────────────────────────────────────────────────────────
HYPOTHESIS_PROMPT_JA = """あなたは創薬ターゲット同定と仮説生成を専門とする薬化学者・計算生物学者です。

以下の構造化エビデンスに基づき、{disease}に対する{gene}をターゲットとした創薬仮説レポートを日本語で作成してください。

{context}

---

以下のセクション構成に従って、厳密にレポートを作成してください。
各事実には、コンテキスト内のリファレンスタグを引用してください（例: [Paper 1]、[ClinVar 2]、[GWAS 1]、[OpenTargets]、[UniProt] など）。
論文のアブストラクト要約はメカニズム考察の根拠として積極的に引用・要約してください。
コンテキストにないデータや参照番号を創作しないでください。

---

## 1. ターゲット妥当性評価（Target Validity Assessment）

### 1a. 遺伝的関連エビデンス（Genetic Association）
{gene}と{disease}を結ぶ遺伝的証拠を評価してください：
- GWAS所見：有意な遺伝子座、効果量、リードSNP [GWASリファレンスを引用]
- ClinVarの病的・病的疑い変異 [ClinVarリファレンスを引用]
- 遺伝的エビデンス信頼度：**低 / 中 / 高 / 非常に高** — 一文で根拠を示す

### 1b. 機能的関連エビデンス（Functional Association）
- {gene}の機能が{disease}の生物学とどう関連するか
- 文献からの主要メカニズム：関連するアブストラクト所見を要約・引用 [Paper リファレンスを引用]
- 疾患作用機序（MoA）：{gene}の疾患病態への既知・推定分子的役割
- 機能的エビデンス信頼度：**低 / 中 / 高 / 非常に高**

### 1c. 臨床的関連性（Clinical Relevance）
- {gene}をターゲットとする既存薬・臨床候補 [ChEMBLリファレンスを引用]
- 臨床開発段階と治療領域との重複
- 患者集団との関連性とアンメットニーズ

### 1d. 遺伝子機能・発現プロファイル・ネットワークコンテキスト
- {gene}の中核的分子機能（キナーゼ/受容体/転写因子等）[UniProtを引用]
- 組織発現プロファイル：{gene}が最も高く発現する組織・臓器と{disease}関連組織との関係 [GTEx・HPAを引用]
- タンパク質細胞内局在およびタンパク質クラス [HPA・UniProtを引用]
- 主要なタンパク質相互作用とその疾患関連性 [IntAct/SIGNORを引用]
- ネットワーク解析のパスウェイ・GOエンリッチメント結果 [エンリッチメント結果を引用]
- {disease}関連生物学におけるネットワークレベルの役割

### 1e. ターゲット妥当性総合評価
**総合：低 / 中 / 高 / 非常に高**
上記4つのサブ評価を統合した一段落の総括。

---

## 2. 疾患メカニズムの考察
- {gene}の機能異常 → {disease}表現型に至るメカニズムモデル（ステップバイステップ）
- 文献アブストラクトによる支持エビデンス [Paperリファレンスを引用]
- 上流活性化因子と下流エフェクター
- 関連する生物学的経路

---

## 3. 治療仮説

### 3a. 分子メカニズム仮説
提案する分子メカニズムを箇条書きで記述してください：
- **機能異常**: {disease}において{gene}はどのように異常を来しているか？ [エビデンスを引用]
- **主要パスウェイ**: 最も重要な下流パスウェイ・エフェクターはどれか？ [エビデンスを引用]
- **細胞レベルの結果**: その異常はどのような細胞表現型を引き起こすか？ [エビデンスを引用]
- **組織・臓器レベルの結果**: その細胞表現型が{disease}の病態にどうつながるか？ [エビデンスを引用]

### 3b. 疾患治療仮説
治療仮説を以下の構造で箇条書きにしてください：
- **介入方法**: {gene}に対してどのような介入を提案するか？（阻害 / 活性化 / 分解 / 補充 / 遺伝子補正 など）
- **ターゲットへの効果**: その介入によって{gene}の活性・発現はどう変わるか？
- **下流への効果**: そのターゲット変化が上記の分子メカニズム異常をどう是正するか？
- **期待される臨床効果**: {disease}の症状・進行においてどのような改善が得られるか？
- **支持するエビデンス**: この論理を支持する主要なエビデンス [Paper・GWAS・ClinVar・OpenTargetsを引用]
- **検証可能な予測**: この仮説を支持または否定できる反証可能な予測を1つ示す

### 3c. 仮説ステートメント（一文）
3a〜3bを一文にまとめてください：
「{gene}に対して〔具体的介入〕を行うと、{disease}患者において〔期待される治療効果〕が得られる。これは〔メカニズムの根拠〕による。」

---

## 4. モダリティ提案
最適な創薬モダリティ（低分子阻害薬/活性化薬・抗体・PROTAC・ASO・遺伝子療法など）を推奨：
- ターゲット生物学・細胞内局在に基づく根拠 [UniProt・HPAを引用]
- 組織特異性の考慮：{gene}が広範に発現する場合は選択性戦略を、組織限局発現の場合はデリバリー上の優位性を論じる [GTEx・HPAを引用]
- AlphaFold構造信頼度とドラッガビリティへの示唆 [AlphaFoldが利用可能な場合引用]
- 主要な技術的考慮事項
- 他のモダリティに対する優位性

---

## 5. 既存薬景観とリポジショニング機会
- {gene}をターゲットとする薬剤・候補化合物の概要 [ChEMBLを引用]
- リポジショニング候補または併用戦略
- 競合状況上の考慮事項

---

## 6. 安全性・毒性リスク評価
- 遺伝子機能に基づくオンターゲット安全性懸念 [UniProt、PubChemを引用]
- 組織発現と安全性：心臓・肝臓・腎臓・CNS・生殖器官における{gene}発現と、それに伴うリスクを具体的に論じる [GTEx・HPAを引用]
- 集団制約スコア：pLI/LOEUFが高い場合はオンターゲットリスクが大きいことを考察 [gnomADが利用可能な場合引用]
- オフターゲット・機序由来毒性シグナル [副作用データを引用]
- 患者集団上のリスク要因

---

## 7. 推奨次期実験
優先度順に3〜5件の実験を列挙：
| 実験種別 | エンドポイント | 期待される結果 |
|---|---|---|
| ... | ... | ... |

---

## 8. 主要な不確実性・限界
- エビデンスのギャップ
- 代替仮説
- 仮説に対する主要リスク

---

## 参考文献

### 論文（PubMed）
（上記で引用した [Paper N] をコンテキストの References セクションからそのまま転記）

### 疾患・遺伝的関連データベース
（上記で引用した [ClinVar N]、[GWAS N]、[OpenTargets] を転記）

### 遺伝子・タンパク質情報データベース
（上記で引用した [UniProt]、[IntAct]、[SIGNOR] を転記）

### 薬剤・安全性データベース
（上記で引用した [ChEMBL N]、[PubChem] を転記）
"""

# ──────────────────────────────────────────────────────────────────────────────
# Presentation evaluation prompts (unchanged structure, updated categories)
# ──────────────────────────────────────────────────────────────────────────────
PRESENTATION_EVAL_PROMPT = """You are a drug discovery expert preparing a concise slide-ready evaluation for a scientific presentation.

Gene: {gene}
Disease: {disease}

Evidence data:
{context}

---

For each evidence category below, provide:
- A rating: ✅ Strong / 🟡 Moderate / 🔴 Weak / ⬜ No data
- One sentence of key finding (max 20 words)

Output ONLY the following JSON (no markdown, no extra text):

{{
  "genetic_association": {{
    "rating": "🟡 Moderate",
    "finding": "GWAS/ClinVar evidence summary"
  }},
  "functional_association": {{
    "rating": "✅ Strong",
    "finding": "Literature/MoA evidence summary"
  }},
  "clinical_relevance": {{
    "rating": "🟡 Moderate",
    "finding": "Existing drugs / clinical stage summary"
  }},
  "network_context": {{
    "rating": "✅ Strong",
    "finding": "PPI network and pathway enrichment summary"
  }},
  "target_validity_overall": {{
    "rating": "✅ Strong",
    "finding": "Overall target validity summary"
  }},
  "repositioning_potential": {{
    "rating": "⬜ No data",
    "finding": "one sentence"
  }},
  "safety_risk": {{
    "rating": "🟡 Moderate",
    "finding": "one sentence"
  }},
  "modality_fit": {{
    "rating": "✅ Strong",
    "finding": "one sentence"
  }},
  "overall_confidence": {{
    "rating": "🟡 Moderate",
    "finding": "one sentence summary of the hypothesis"
  }}
}}
"""

PRESENTATION_EVAL_PROMPT_JA = """あなたは創薬の専門家です。科学プレゼンテーション用に、各エビデンスカテゴリを簡潔に評価してください。

遺伝子: {gene}
疾患: {disease}

エビデンスデータ:
{context}

---

各カテゴリについて以下を日本語で提供してください：
- 評価: ✅ 強い / 🟡 中程度 / 🔴 弱い / ⬜ データなし
- 一文の所見（最大30字）

以下のJSONのみを出力（マークダウン・余分なテキスト不要）:

{{
  "genetic_association": {{
    "rating": "🟡 中程度",
    "finding": "GWAS/ClinVarエビデンスの要約"
  }},
  "functional_association": {{
    "rating": "✅ 強い",
    "finding": "文献・MoAエビデンスの要約"
  }},
  "clinical_relevance": {{
    "rating": "🟡 中程度",
    "finding": "既存薬・臨床段階の要約"
  }},
  "network_context": {{
    "rating": "✅ 強い",
    "finding": "PPIネットワーク・パスウェイ解析の要約"
  }},
  "target_validity_overall": {{
    "rating": "✅ 強い",
    "finding": "ターゲット妥当性の総合要約"
  }},
  "repositioning_potential": {{
    "rating": "⬜ データなし",
    "finding": "一文"
  }},
  "safety_risk": {{
    "rating": "🟡 中程度",
    "finding": "一文"
  }},
  "modality_fit": {{
    "rating": "✅ 強い",
    "finding": "一文"
  }},
  "overall_confidence": {{
    "rating": "🟡 中程度",
    "finding": "仮説の一文要約"
  }}
}}
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


def generate_presentation_eval(
    gene: str,
    disease: str,
    context: str,
    llm_client,
    lang: str = "en",
) -> dict:
    import json, re
    template = PRESENTATION_EVAL_PROMPT_JA if lang == "ja" else PRESENTATION_EVAL_PROMPT
    prompt = template.format(gene=gene, disease=disease, context=context)
    raw = llm_client.generate(prompt, temperature=0.1, max_tokens=1000)
    m = re.search(r"\{.*\}", raw, re.DOTALL)
    if m:
        return json.loads(m.group())
    return {}
