# 創薬仮説生成システム 仕様書

ローカル LLM（Ollama）と公開バイオデータベースを組み合わせ、
**「遺伝子 × 疾患」ごとに創薬仮説レポートを自動生成する**システム。

---

## 1. 概要

- 入力: 疾患名（1件）と対象遺伝子リスト（複数可）
- 処理: 15以上の公開DBからエビデンスを収集 → PPIネットワーク・機能解析 → LLMで仮説生成・評価
- 出力: 遺伝子ごとの Markdown レポート（評価表・PPI図・エンリッチメント・仮説本文）＋ バッチサマリー
- 特徴:
  - 完全ローカル LLM（API課金なし、データ外部送信なし）
  - 使用データは原則 **CC BY 4.0（商用利用可）**。KEGG / TRANSFAC / BioGRID 等の商用・非商用限定ソースは既定で不使用

---

## 2. アーキテクチャ

```
batch_hypothesis_generator.ipynb   ← UI（疾患選択・遺伝子入力・実行）
        │
        ▼
   pipeline.py                     ← 処理フロー全体の制御
        │
   ┌────┼─────────────┬───────────────┬──────────────┐
   ▼    ▼             ▼               ▼              ▼
aggregator   network        hypothesis        report      llm/
 (収集+      (PPI+          (プロンプト+      (Markdown/   ollama_client
  整形)       enrichment)    LLM生成+評価)     HTML整形)    (LLM呼び出し)
   │
   ▼
collectors/ × 18                   ← 各DBのAPIラッパー（単一責務）
```

### モジュール責務

| モジュール | 役割 |
|---|---|
| `collectors/*.py` | 各DB APIを叩き、辞書/リストで返す薄いラッパー |
| `aggregator.py` | 全コレクターを並列実行（`collect_all`）／ LLM用コンテキスト整形（`build_llm_context`） |
| `network.py` | PPIネットワーク構築・パートナー順位付け・エンリッチメント・PNG描画 |
| `hypothesis.py` | プロンプト定義・仮説本文生成・Target Validity 評価生成 |
| `report.py` | 評価表・エンリッチメント表・サマリーの Markdown / HTML 生成（文字列のみ） |
| `pipeline.py` | 1遺伝子の処理（`process_gene`）とバッチ実行（`run_batch`）・保存 |
| `llm/ollama_client.py` | Ollama への HTTP 呼び出し（ストリーミング・JSON強制に対応） |

---

## 3. データソース一覧

すべて API キー不要（BioGRID を除く）。括弧内はライセンス。

| ソース | 取得内容 | ライセンス / 商用 |
|---|---|---|
| PubMed | 論文・アブストラクト | パブリックドメイン |
| OpenTargets | 遺伝子×疾患 関連スコア・既知薬 | CC0 / 可 |
| UniProt | タンパク質機能・局在・GO | CC BY 4.0 / 可 |
| GWAS Catalog | 遺伝的関連研究 | EMBL-EBI / 可 |
| ClinVar | 病的バリアント | パブリックドメイン |
| ChEMBL | 既存薬・フェーズ・作用機序 | CC BY-SA 3.0 / 可 |
| IntAct | タンパク質相互作用 (PPI) | CC BY 4.0 / 可 |
| SIGNOR | シグナル伝達相互作用 (PPI) | CC BY 4.0 / 可 |
| Reactome | パスウェイ・相互作用 (PPI) | CC BY 4.0 / 可 |
| gnomAD | 制約スコア pLI / LOEUF | 可 |
| GTEx | 組織別発現 TPM | 可 |
| Human Protein Atlas | 発現・細胞内局在 | CC BY-SA / 可 |
| DGIdb | 薬剤–遺伝子相互作用 | 可 |
| ClinicalTrials.gov | 臨床試験 | パブリックドメイン |
| AlphaFold DB | 構造信頼度 pLDDT | CC BY 4.0 / 可 |
| PubChem / openFDA | 毒性アッセイ・副作用 | パブリックドメイン |
| g:Profiler | 機能エンリッチメント | BSD 2-Clause / 可 |

### 不使用（ライセンス上の理由）

- **KEGG** … 商用利用に別途ライセンスが必要 → g:Profiler のソースから除外
- **TRANSFAC (TF)** … 商用ライセンス → 除外
- **BioGRID** … 非商用・学術限定 → 既定で無効（`use_biogrid=True` で明示的に有効化は可能）

---

## 4. 処理フロー（遺伝子1件あたり）

`pipeline.process_gene()` が以下を順に実行する。

1. **エビデンス収集** `aggregator.collect_all`
   - 全コレクターを最大10並列で実行（各コレクターは最大3回リトライ）
2. **PPIネットワーク構築** `network.build_ppi_network`
   - IntAct + SIGNOR + Reactome の相互作用を統合
3. **エンリッチメント解析** `network.run_network_enrichment`
   - ネットワーク内遺伝子を g:Profiler で解析（FDR < 0.05）
4. **コンテキスト整形** `aggregator.build_llm_context`
   - 収集結果を引用タグ付きの構造化テキストに変換
5. **仮説生成** `hypothesis.generate_hypothesis`
   - レポート本文（8セクション）をストリーミング生成
6. **Target Validity 評価** `hypothesis.generate_presentation_eval`
   - **生成した本文の結論に基づき**評価（表と本文を一致させる）
   - Ollama の `format="json"` で JSON を強制取得し、表記ゆれを正規化
7. **レポート保存** `pipeline._save_report`
   - Markdown 本体 + PPI画像(PNG) + eval/raw JSON を保存

---

## 5. PPI パートナーの選定ロジック

`network.rank_partners()` が次の優先順位（降順）で並べ替え、上位30件を採用する。

1. **複数DBで共通する遺伝子**（裏付けDB数が多いほど上位）
2. **スコアが高いもの**（IntAct intactScore / SIGNOR score）
3. エッジの重み（観測された相互作用の回数）

---

## 6. 出力仕様

### 個別レポート `reports/{GENE}_{DISEASE}/{timestamp}_{EN|JA}.md`

```
# Drug Discovery Hypothesis: GENE × DISEASE
---
## Target Validity              ← 評価表（VH/H/M/L、本文の結論と一致）
---
## PPI Network                  ← PNG（★=対象遺伝子が上部、パートナーが下部）
## Functional Enrichment        ← ソース別テーブル（GO/Reactome/WikiPathways…）
---
## 1〜8（仮説本文）             ← 妥当性/機序/治療仮説/モダリティ/既存薬/安全性/実験/不確実性
---
## Evidence Context             ← LLMに渡した根拠（引用タグ付き）
```

同フォルダに以下も保存:
- `{timestamp}_ppi.png` … PPIネットワーク画像
- `{timestamp}_eval.json` … 評価結果
- `{timestamp}_raw.json` … 収集した生データ

### Target Validity 評価基準（4段階）

| 略号 | 基準 |
|---|---|
| **VH** (Very High) | 遺伝子と疾患の関連が明確、かつ疾患の悪性度・進行と関連 |
| **H** (High) | 関連が明確。遺伝子の機能・パスウェイが疾患機序と直接関連 |
| **M** (Middle) | 直接の関連情報はないが、PPI/パスウェイ重複で説明可能 |
| **L** (Low) | 関連情報なし |

評価項目: Genetic Association / Functional Association / Clinical Relevance /
Expression・Network / Overall（本文 Section 1 の 1a〜1e に対応）

### バッチサマリー `reports/{DISEASE}_summary_{timestamp}.md`

行 = 評価項目、列 = 遺伝子 の一覧表（HTMLでも画面表示）。

---

## 7. 設定項目（ノートブック Step 1）

```python
MODEL = 'qwen2.5:14b'   # Ollama モデル名
LANG  = 'en'            # 'ja' or 'en'

CONTEXT_CONFIG = dict(
    max_papers=5, abstract_chars=600,   # PubMed
    max_drugs=8, max_gwas=5, max_clinvar=5,
    max_interactions=10, max_trials=6, max_reactome=10,
    gtex_top_n=5, hpa_top_n=8, max_dgidb=8,
    uniprot_chars=500, uniprot_keywords=10, uniprot_go_terms=8,
)
```
値を大きくすると根拠が豊富になるが、生成時間も増える。

---

## 8. 使い方

### 前提

```bash
pip install -r requirements.txt
ollama serve            # 別ターミナルで起動
ollama pull qwen2.5:14b
```

### 実行（`batch_hypothesis_generator.ipynb`）

1. **Step 1**: モデル・言語・コンテキスト設定を実行
2. **Step 2**: 疾患名を検索しリストから選択（OpenTargets）
3. **Step 3**: `GENE_LIST_INPUT` に遺伝子を入力して実行（HGNCで検証）
4. **Step 4**: バッチ実行 → `reports/` にレポートが保存され、サマリーが表示される

---

## 9. 拡張方法

### データソースを追加する

1. `collectors/xxx.py` に `get_xxx(...)` を実装（dict/list を返す）
2. `aggregator.collect_all` の `tasks` に1行追加
3. 必要なら `aggregator.build_llm_context` に整形ブロックを追加

### PPIソースを追加する

1. `collectors/xxx.py` に `get_interactions()` を実装
   （`{"source","target","score"}` または `{"partners":[...]}` 形式で返す）
2. `network.build_ppi_network` に `add_edges(...)` 呼び出しを追加
3. `network._db_color` に色を追加

---

## 10. 依存関係

- Python 3.10+
- `requests`, `networkx`, `matplotlib`, `ipywidgets`
- Ollama（ローカル LLM ランタイム）
