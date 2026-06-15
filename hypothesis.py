"""Hypothesis generation — builds prompt and calls LLM."""

HYPOTHESIS_PROMPT_TEMPLATE = """You are an expert medicinal chemist and computational biologist specializing in drug target identification and drug discovery hypothesis generation.

Based on the evidence below, generate a detailed drug discovery hypothesis for targeting the gene {gene} in the context of {disease}.

{context}

---

Please provide a structured hypothesis report with the following sections:

## 1. Target Validity Assessment
- Summarize the strength of evidence linking {gene} to {disease}
- Assign a confidence score (Low / Medium / High / Very High) with justification
- Key supporting evidence points (genetic, functional, clinical)

## 2. Proposed Disease Mechanism
- How does dysregulation/mutation of {gene} contribute to {disease} pathophysiology?
- Relevant biological pathways and network context
- Upstream/downstream molecular players

## 3. Therapeutic Hypothesis
State a clear, testable hypothesis in the format:
"If [mechanism of intervention on {gene}] then [expected therapeutic outcome] because [mechanistic rationale]"

## 4. Modality Recommendation
Recommend the most suitable drug modality (small molecule, antibody, PROTAC, ASO, gene therapy, etc.) and explain:
- Why this modality fits the target biology
- Key technical considerations
- Potential advantages over alternative modalities

## 5. Existing Drug Landscape & Repositioning Opportunities
- Summarize existing drugs/compounds targeting {gene}
- Identify repositioning candidates or combination strategies
- Competitive/freedom-to-operate considerations

## 6. Safety & Toxicity Risk Assessment
- Key on-target safety concerns based on gene function
- Off-target/mechanism-based toxicity risks
- Patient population considerations

## 7. Recommended Next Experiments
List 3–5 priority experiments to validate this hypothesis (in vitro, in vivo, translational):
- Experiment type | Endpoint | Expected result

## 8. Key Uncertainties & Limitations
- What evidence is missing?
- Alternative interpretations of the data
- Major risks to the hypothesis

Be specific and maintain scientific rigor. Do not fabricate data points not present in the evidence.

**Citation rules:**
- Cite evidence inline using the [Ref N] tags provided in the context (e.g., "strong GWAS association [Ref 3]").
- Use the abstract summaries in the "Recent Literature" section to support mechanistic reasoning — quote or paraphrase key findings from abstracts where relevant.
- At the end of the report, add a "## References" section that lists every [Ref N] you cited, copied verbatim from the References block in the context.
- Do not invent reference numbers not present in the context.
"""

HYPOTHESIS_PROMPT_JA = """あなたは創薬ターゲット同定と創薬仮説生成を専門とする薬化学者・計算生物学者です。

以下のエビデンスに基づき、{disease}に対する{gene}をターゲットとした創薬仮説を日本語で詳述してください。

{context}

---

以下の構成でレポートを作成してください：

## 1. ターゲット妥当性評価
- {gene}と{disease}を結びつけるエビデンスの強度を要約
- 確信度スコア（低 / 中 / 高 / 非常に高）を根拠とともに提示
- 主要な支持エビデンス（遺伝的・機能的・臨床的）

## 2. 疾患メカニズムの考察
- {gene}の機能異常・変異が{disease}の病態生理にどう寄与するか
- 関連する生物学的経路とネットワーク
- 上流・下流の分子メカニズム

## 3. 治療仮説
以下の形式で明確かつ検証可能な仮説を述べてください：
「{gene}に対して〔介入機序〕を行うと、〔期待される治療効果〕が得られる。これは〔メカニズムの根拠〕による。」

## 4. モダリティ提案
最適な創薬モダリティ（低分子・抗体・PROTAC・ASO・遺伝子療法など）を推奨し、以下を説明：
- このモダリティがターゲット生物学に適合する理由
- 主要な技術的考慮事項
- 他のモダリティに対する優位性

## 5. 既存薬景観とリポジショニング機会
- {gene}をターゲットとする既存薬・化合物の概要
- リポジショニング候補または併用戦略
- 競合状況・知的財産上の考慮

## 6. 安全性・毒性リスク評価
- 遺伝子機能に基づくオンターゲット安全性懸念
- オフターゲット・機序由来毒性リスク
- 患者集団上の考慮事項

## 7. 推奨次期実験
仮説を検証するための優先実験を3〜5件列挙（in vitro・in vivo・トランスレーショナル）：
- 実験種別 | エンドポイント | 期待される結果

## 8. 主要な不確実性・限界
- 不足しているエビデンスは何か
- データの代替解釈
- 仮説に対する主要なリスク

提供されたデータに基づき具体的に記述し、科学的厳密性を保ってください。データにない情報を捏造しないでください。

**引用のルール:**
- コンテキスト内の [Ref N] タグを使って根拠を本文中に引用してください（例: 「強いGWAS関連性が確認されている [Ref 3]」）。
- 「Recent Literature」セクションのアブストラクト要約を積極的に活用し、メカニズムの考察や治療仮説の根拠として引用・要約してください。
- レポートの末尾に「## 参考文献」セクションを追加し、引用した [Ref N] をコンテキストの References ブロックからそのまま転記してください。
- コンテキストに存在しない参照番号を創作しないでください。
"""

PRESENTATION_EVAL_PROMPT = """You are a drug discovery expert preparing a concise slide-ready evaluation for a scientific presentation.

Gene: {gene}
Disease: {disease}

Evidence data:
{context}

---

For each evidence category below, provide a ONE-LINE evaluation with:
- A rating: ✅ Strong / 🟡 Moderate / 🔴 Weak / ⬜ No data
- One sentence of key finding (max 20 words)

Output ONLY the following JSON structure (no markdown, no extra text):

{{
  "target_validity": {{
    "rating": "✅ Strong",
    "finding": "one sentence"
  }},
  "genetic_evidence": {{
    "rating": "🟡 Moderate",
    "finding": "one sentence"
  }},
  "disease_mechanism": {{
    "rating": "✅ Strong",
    "finding": "one sentence"
  }},
  "existing_drugs": {{
    "rating": "🟡 Moderate",
    "finding": "one sentence"
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


def generate_hypothesis(
    gene: str,
    disease: str,
    context: str,
    llm_client,
    temperature: float = 0.3,
    lang: str = "en",
) -> str:
    template = HYPOTHESIS_PROMPT_JA if lang == "ja" else HYPOTHESIS_PROMPT_TEMPLATE
    prompt = template.format(gene=gene, disease=disease, context=context)
    return llm_client.generate(prompt, temperature=temperature, max_tokens=4096)


PRESENTATION_EVAL_PROMPT_JA = """あなたは創薬の専門家です。科学プレゼンテーション用に、各エビデンスカテゴリを簡潔に評価してください。

遺伝子: {gene}
疾患: {disease}

エビデンスデータ:
{context}

---

以下の各カテゴリについて、評価と一文の所見（最大30字）を日本語で提供してください。
評価: ✅ 強い / 🟡 中程度 / 🔴 弱い / ⬜ データなし

以下のJSON構造のみを出力（マークダウン・余分なテキスト不要）:

{{
  "target_validity": {{"rating": "✅ 強い", "finding": "一文"}},
  "genetic_evidence": {{"rating": "🟡 中程度", "finding": "一文"}},
  "disease_mechanism": {{"rating": "✅ 強い", "finding": "一文"}},
  "existing_drugs": {{"rating": "🟡 中程度", "finding": "一文"}},
  "repositioning_potential": {{"rating": "⬜ データなし", "finding": "一文"}},
  "safety_risk": {{"rating": "🟡 中程度", "finding": "一文"}},
  "modality_fit": {{"rating": "✅ 強い", "finding": "一文"}},
  "overall_confidence": {{"rating": "🟡 中程度", "finding": "仮説の一文要約"}}
}}
"""


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
    raw = llm_client.generate(prompt, temperature=0.1, max_tokens=800)
    m = re.search(r"\{.*\}", raw, re.DOTALL)
    if m:
        return json.loads(m.group())
    return {}
