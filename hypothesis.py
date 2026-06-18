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

## 2. Disease Mechanism
Step-by-step model: how does {gene} dysfunction lead to {disease}? Use specific pathway names, effectors, and cite literature.

---

## 3. Therapeutic Hypothesis

### 3a. Molecular Mechanism (bullet points)
- **Dysregulation:** how is {gene} altered in {disease}?
- **Key pathway:** which downstream pathway is most affected?
- **Cellular effect:** what cell-level phenotype results?
- **Tissue/organ effect:** how does this drive {disease} pathology?

### 3b. Treatment Hypothesis (bullet points)
- **Intervention:** proposed modulation of {gene} (inhibit / activate / degrade / replace)
- **Effect on target:** how would this change {gene} activity?
- **Downstream effect:** how does that correct the mechanism above?
- **Clinical outcome:** expected improvement in {disease}
- **Key evidence:** cite the strongest supporting data points
- **Falsifiable prediction:** one testable statement that could confirm or refute this

### 3c. One-sentence hypothesis
"If [intervention on {gene}] then [outcome in {disease} patients] because [mechanism]."

---

## 4. Modality
Best modality (small molecule / antibody / PROTAC / ASO / gene therapy) and why:
- Rationale based on subcellular location and druggability
- Tissue specificity strategy (selectivity vs. delivery advantage)
- AlphaFold pLDDT and structural confidence
- Key advantage over alternatives

---

## 5. Existing Drug Landscape
Summarise known drugs/candidates for {gene} with phases. Repositioning or combination opportunities.

---

## 6. Safety Assessment
- On-target risks from gene function
- Safety-relevant tissue expression (heart/liver/kidney/CNS) with TPM or HPA level
- gnomAD pLI/LOEUF interpretation (high = greater on-target risk)
- Off-target / mechanism-based toxicity signals

---

## 7. Recommended Experiments
| Experiment | Endpoint | Expected Result |
|---|---|---|
| (3–5 specific experiments based on the evidence above) | | |

---

## 8. Key Uncertainties
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

## 2. 疾患メカニズムの考察
{gene}の機能異常が{disease}表現型につながるステップを具体的に記述。パスウェイ名・エフェクター・文献を引用。

---

## 3. 治療仮説

### 3a. 分子メカニズム仮説（箇条書き）
- **機能異常：** {disease}において{gene}はどう変化しているか（具体的に）
- **主要パスウェイ：** 最も重要な下流パスウェイ・エフェクター
- **細胞への影響：** それがもたらす細胞表現型
- **組織・疾患への影響：** {disease}病態へのつながり

### 3b. 疾患治療仮説（箇条書き）
- **介入方法：** {gene}への具体的介入（阻害/活性化/分解/補充など）
- **ターゲットへの効果：** 介入による{gene}活性・発現の変化
- **下流効果：** 上記メカニズム異常をどう是正するか
- **期待される臨床効果：** {disease}の症状・進行への改善
- **支持エビデンス：** 最も強いエビデンスを引用タグで
- **検証可能な予測：** 反証可能な1つの予測

### 3c. 仮説一文
「{gene}に対して〔介入〕を行うと、{disease}患者において〔臨床効果〕が得られる。これは〔メカニズム〕による。」

---

## 4. モダリティ提案
最適モダリティ（低分子/抗体/PROTAC/ASO/遺伝子療法）と理由：
- 細胞内局在・ドラッガビリティに基づく根拠
- 組織特異性（選択性戦略 or デリバリー優位性）
- AlphaFold pLDDT・構造信頼度
- 他モダリティに対する優位性

---

## 5. 既存薬景観とリポジショニング
{gene}標的薬・候補をフェーズ付きで整理。リポジショニング・併用の可能性。

---

## 6. 安全性・毒性リスク評価
- 遺伝子機能に基づくオンターゲットリスク
- 安全性関連組織における発現（心臓/肝臓/腎臓/CNS）とTPM値またはHPAレベル
- gnomAD pLI/LOEUF（高値＝オンターゲットリスク大）の解釈
- オフターゲット・機序由来毒性シグナル

---

## 7. 推奨次期実験
| 実験種別 | エンドポイント | 期待される結果 |
|---|---|---|
| （上記エビデンスに基づく具体的な実験を3〜5件） | | |

---

## 8. 主要な不確実性・限界
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

# ──────────────────────────────────────────────────────────────────────────────
# Presentation evaluation prompts
# ──────────────────────────────────────────────────────────────────────────────
PRESENTATION_EVAL_PROMPT = """You are a drug discovery expert. Below is a hypothesis REPORT that was already written for {gene} in {disease}, followed by the source evidence.

Your job: produce a Target Validity scorecard that is CONSISTENT with the report's conclusions. The scorecard MUST be derived from the report's "## 1. Target Validity" section. Read each subsection and translate its written conclusion (including its stated Confidence) into a rating:
- genetic_association     <- from "1a. Genetic Evidence"
- functional_association  <- from "1b. Functional Evidence"
- clinical_relevance      <- from "1c. Clinical Relevance"
- network_context         <- from "1d. Expression & Network"
- target_validity_overall <- from "1e. Overall Validity Score"
Do not contradict the report. If a subsection says evidence is weak/absent, the rating must be Low; do not inflate it.

--- REPORT ---
{hypothesis}
--- END REPORT ---

--- SOURCE EVIDENCE ---
{context}
--- END SOURCE EVIDENCE ---

RATING SCALE — choose exactly one per field:
- Very High : gene-disease association is clear AND the gene is linked to disease severity/progression
- High      : gene-disease association is clear; gene function/pathway directly related to disease mechanism
- Middle    : no direct gene-disease data, but association is explainable via PPI partners or pathway overlap
- Low       : no relevant information found

Output ONLY this JSON (no markdown, no extra text). Fill every field with real data:

{{
  "genetic_association":    {{"rating": "Very High / High / Middle / Low", "finding": "specific finding in ≤20 words"}},
  "functional_association": {{"rating": "Very High / High / Middle / Low", "finding": "specific finding in ≤20 words"}},
  "clinical_relevance":     {{"rating": "Very High / High / Middle / Low", "finding": "specific finding in ≤20 words"}},
  "network_context":        {{"rating": "Very High / High / Middle / Low", "finding": "specific finding in ≤20 words"}},
  "target_validity_overall":{{"rating": "Very High / High / Middle / Low", "finding": "one-sentence overall summary"}}
}}
"""

PRESENTATION_EVAL_PROMPT_JA = """あなたは創薬の専門家です。以下は{gene}の{disease}に対してすでに作成された仮説レポートと、その根拠エビデンスです。

あなたの仕事：レポートの結論と【矛盾しない】Target Validity スコアカードを作成すること。スコアカードは必ずレポートの「## 1. ターゲット妥当性評価」セクションから導出してください。各小項目の記述（明記された信頼度を含む）を読み取り、評価に変換します：
- genetic_association     ← 「1a. 遺伝的エビデンス」
- functional_association  ← 「1b. 機能的エビデンス」
- clinical_relevance      ← 「1c. 臨床的関連性」
- network_context         ← 「1d. 発現プロファイル・ネットワーク」
- target_validity_overall ← 「1e. ターゲット妥当性総合スコア」
レポートと食い違う評価を出してはいけません。小項目が「エビデンスが弱い/ない」と述べている場合は必ず Low とし、過大評価しないこと。

--- レポート ---
{hypothesis}
--- レポートここまで ---

--- 根拠エビデンス ---
{context}
--- 根拠エビデンスここまで ---

評価基準（各フィールドに必ずいずれか1つを選択）:
- Very High : 遺伝子と疾患の関連性が明らかで、かつ疾患の悪性度・進行と関連している
- High      : 遺伝子と疾患の関連性が明確。遺伝子の機能・パスウェイが疾患メカニズムと直接関連
- Middle    : 直接的な関連情報はないが、PPIパートナーやパスウェイの重複で説明可能
- Low       : 関連情報が全くない

以下のJSONのみを出力（マークダウン・余分なテキスト不要）。実際のデータを使うこと：

{{
  "genetic_association":    {{"rating": "Very High / High / Middle / Low", "finding": "具体的な知見を30字以内で"}},
  "functional_association": {{"rating": "Very High / High / Middle / Low", "finding": "具体的な知見を30字以内で"}},
  "clinical_relevance":     {{"rating": "Very High / High / Middle / Low", "finding": "具体的な知見を30字以内で"}},
  "network_context":        {{"rating": "Very High / High / Middle / Low", "finding": "具体的な知見を30字以内で"}},
  "target_validity_overall":{{"rating": "Very High / High / Middle / Low", "finding": "総合評価を一文で"}}
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


# 評価カードの正規キー
EVAL_FIELDS = [
    "genetic_association",
    "functional_association",
    "clinical_relevance",
    "network_context",
    "target_validity_overall",
]

# rating のゆれを正規ラベルへ寄せる
_RATING_CANON = [
    (r"very[\s_-]*high|非常に高", "Very High"),
    (r"\bhigh\b|^高",            "High"),
    (r"middle|moderate|中",      "Middle"),
    (r"\blow\b|^低|weak|弱|none|no[\s_-]*data|なし", "Low"),
]


def _normalize_rating(value: str) -> str:
    """rating 文字列を Very High / High / Middle / Low のいずれかに正規化。"""
    import re
    s = (value or "").strip()
    # プレースホルダ（"Very High / High / ..." など）はそのまま返してきたら未判定扱い
    if not s or "/" in s or s.upper() == "RATING":
        return ""
    low = s.lower()
    for pat, canon in _RATING_CANON:
        if re.search(pat, low):
            return canon
    return ""


def _extract_json(raw: str) -> dict:
    """LLM出力から JSON オブジェクトを頑健に抽出する。"""
    import json, re
    if not raw:
        return {}
    text = raw.strip()
    # ```json ... ``` フェンスを除去
    text = re.sub(r"^```(?:json)?\s*|\s*```$", "", text, flags=re.IGNORECASE).strip()
    # 最初の { から対応する } までをバランス抽出
    start = text.find("{")
    if start == -1:
        return {}
    depth, end = 0, -1
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    candidate = text[start:end] if end != -1 else text[start:]
    # 末尾カンマを除去
    candidate = re.sub(r",(\s*[}\]])", r"\1", candidate)
    try:
        return json.loads(candidate)
    except json.JSONDecodeError:
        return {}


def generate_presentation_eval(
    gene: str,
    disease: str,
    context: str,
    llm_client,
    lang: str = "en",
    hypothesis: str = "",
) -> dict:
    """Target Validity 評価カード（5項目）を生成し、正規化して返す。

    hypothesis を渡すと、その本文の結論と整合する評価を生成する（表と文章を一致させる）。
    Ollama の format="json" で JSON を強制し、抽出・正規化で表記ゆれを吸収する。
    """
    template = PRESENTATION_EVAL_PROMPT_JA if lang == "ja" else PRESENTATION_EVAL_PROMPT
    prompt = template.format(
        gene=gene, disease=disease, context=context,
        hypothesis=hypothesis or "(report not provided — judge from evidence only)",
    )
    try:
        raw = llm_client.generate(prompt, temperature=0.1, max_tokens=1000, format="json")
    except TypeError:
        # format 未対応クライアントへのフォールバック
        raw = llm_client.generate(prompt, temperature=0.1, max_tokens=1000)

    parsed = _extract_json(raw)
    if not parsed:
        return {}

    result = {}
    for key in EVAL_FIELDS:
        item = parsed.get(key)
        if not isinstance(item, dict):
            continue
        rating  = _normalize_rating(str(item.get("rating", "")))
        finding = str(item.get("finding", "")).strip()
        # プレースホルダ的な finding は捨てる
        if finding in {"...", "FINDING", "finding", "RATING"}:
            finding = ""
        result[key] = {"rating": rating, "finding": finding}
    return result
